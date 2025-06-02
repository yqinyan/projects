program compute_phi
    use mpi_f08
    use decomp_2d
    implicit none

    integer :: Nx, Ny, Nz
    integer :: i, j, k, ierr, nproc, rank
    integer :: p_row, p_col
    real(kind=8), allocatable :: phi(:, :, :), result(:, :, :)
    real(kind=8) :: lambda, dx, dy, dz
    real(kind=8) :: start_time, end_time, comp_time, comm_time
    real(kind=8) :: total_comp_time, total_comm_time
    integer :: local_size
    real(kind=8), allocatable :: tmp(:)

    ! Taille globale du domaine
    Nx = 128
    Ny = 128
    Nz = 128
    ! Ici, on laisse le decomp2d en mode auto-tuning pour le processeur grid
    p_row = 0
    p_col = 0

    ! Domaine normalise [0,1]
    dx = 1.0 / Nx
    dy = 1.0 / Ny
    dz = 1.0 / Nz
    lambda = 1.0

    ! Initialisation MPI et affichage du nombre de processus
    call MPI_Init(ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, nproc, ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
    if (rank == 0) print *, "MPI_COMM_WORLD size = ", nproc

    ! Initialisation decomp_2d (la fonction auto-tune détermine p_row et p_col)
    call decomp_2d_init(Nx, Ny, Nz, p_row, p_col)

    ! Allocation des tableaux locaux
    allocate(phi(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)))
    allocate(result(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)))

    ! Mesure du temps de calcul
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    start_time = MPI_Wtime()

    ! Definition de phi
    do k = xstart(3), xend(3)
        do j = xstart(2), xend(2)
            do i = xstart(1), xend(1)
                phi(i, j, k) = sin(i * dx) * sin(j * dy) * sin(k * dz)
            end do
        end do
    end do

    ! Calcul du probleme
    result = phi - lambda**2 * phi

    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    end_time = MPI_Wtime()
    comp_time = end_time - start_time

    ! Mesure du temps de communication 
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    start_time = MPI_Wtime()

    ! Copie le tableau local 'result' dans un buffer contigu 'tmp'
    local_size = size(result)
    allocate(tmp(local_size))
    tmp = reshape(result, [local_size])

    ! Réduction MPI : utilisation de MPI_IN_PLACE pour le processus 0
    if (rank == 0) then
        call MPI_Reduce(MPI_IN_PLACE, tmp, local_size, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    else
        call MPI_Reduce(tmp, tmp, local_size, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    endif

    ! Copier le resultat reduit dans 'result' sur le processus 0
    if (rank == 0) then
        result = reshape(tmp, shape(result))
    end if
    deallocate(tmp)

    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    end_time = MPI_Wtime()
    comm_time = end_time - start_time

    ! Collecte des temps de calcul et de communication
    call MPI_Reduce(comp_time, total_comp_time, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    call MPI_Reduce(comm_time, total_comm_time, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

    if (rank == 0) then
        total_comp_time = total_comp_time / nproc
        total_comm_time = total_comm_time / nproc

        print *, "Temps de calcul moyen : ", total_comp_time, " secondes"
        print *, "Temps de communication moyen : ", total_comm_time, " secondes"
    end if

    ! Liberation des ressources
    deallocate(phi)
    deallocate(result)

    ! Finalisation decomp_2d et MPI
    call decomp_2d_finalize()
    call MPI_Finalize(ierr)
end program compute_phi