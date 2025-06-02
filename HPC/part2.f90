program cg_solver
    use mpi_f08
    use decomp_2d
    implicit none

    integer :: nx, ny, nz
    integer :: ierr, nproc, rank
    integer :: p_row, p_col
    real(kind=8), allocatable :: phi(:,:,:), S(:,:,:), r(:,:,:), p(:,:,:), Ap(:,:,:)
    real(kind=8) :: alpha, beta, r_norm, r_norm_new, pAp
    integer :: i, j, k, iter
    real(kind=8) :: t_start, t_end, comp_time, start_time, end_tiem, comm_time
    integer :: source, dest, tag

    ! Taille globale de la grille
    nx = 128
    ny = 128
    nz = 128

    p_row = 0
    p_col = 0

    ! Initialisation de MPI
    call MPI_Init(ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, nproc, ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)

    if (rank == 0) print *, "Taille de MPI_COMM_WORLD = ", nproc

    ! Initialisation de decomp_2d
    call decomp_2d_init(nx, ny, nz, p_row, p_col)

    ! Allocation des tableaux
    allocate(phi(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)))
    allocate(S(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)))
    allocate(r(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)))
    allocate(p(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)))
    allocate(Ap(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)))

    ! Calcul des conditions initiales S = X^3
    do k = xstart(3), xend(3)
        do j = xstart(2), xend(2)
            do i = xstart(1), xend(1)
                S(i,j,k) = real(i, kind=8)**3  
                phi(i,j,k) = 0.0               
            end do
        end do
    end do

    ! Calcul du résidu initial r = S - A * phi
    call apply_A(phi, Ap)
    r = S - Ap
    p = r
    r_norm = sum(r * r)


    do iter = 1, 1000  
        t_start = MPI_Wtime()
        ! Calcul de Ap = A * p
        call apply_A(p, Ap)

        ! Calcul de alpha_k = (r_k, r_k) / (p_k, Ap_k)
        pAp = sum(p * Ap)
        alpha = r_norm / pAp

        ! Mise à jour de phi_k+1 = phi_k + alpha_k * p_k
        phi = phi + alpha * p

        ! Mise à jour de r_k+1 = r_k - alpha_k * Ap_k
        r = r - alpha * Ap

        r_norm_new = sum(r * r)

        ! Vérification de la convergence
        if (sqrt(r_norm_new) < 1.0e-6) exit

        ! Calcul de beta_k = (r_k+1, r_k+1) / (r_k, r_k)
        beta = r_norm_new / r_norm

        ! Mise à jour de p_k+1 = r_k+1 + beta_k * p_k
        p = r + beta * p

        r_norm = r_norm_new
        t_end = MPI_Wtime()
        comp_time = t_end - t_start

        ! P2P Communication
        dest = mod(rank + 1, nproc)
        source = mod(rank - 1, nproc)

        ! Send r to the next process
        call MPI_Send(r, size(r), MPI_REAL8, dest, 0, MPI_COMM_WORLD, ierr)
        start_time = MPI_Wtime()
        ! Receive r from the previous process
        call MPI_Recv(r, size(r), MPI_REAL8, source, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
        t_end = MPI_Wtime()
        ! Calcul du temps de communication
        comm_time = t_end - t_start

    end do


    if (rank == 0) then
        print *, "Méthode du gradient conjugué terminée, nombre d'itérations :", iter
        print *, "Temps de calcul :", comp_time, "secondes"
        print *, "Temps de communication :", comm_time, "secondes"
    end if

    ! Libération des ressources
    deallocate(phi, S, r, p, Ap)
    call decomp_2d_finalize()
    call MPI_Finalize(ierr)

contains

    ! Calcul de A * phi, c'est-à-dire A(phi) = phi - λ² phi
    subroutine apply_A(input_phi, output_Ap)
        real(kind=8), intent(in) :: input_phi(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3))
        real(kind=8), intent(out) :: output_Ap(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3))
        integer :: i, j, k
        real(kind=8) :: lambda2
        lambda2 = 1.0

        do k = xstart(3), xend(3)
            do j = xstart(2), xend(2)
                do i = xstart(1), xend(1)
                    output_Ap(i, j, k) = input_phi(i, j, k) - lambda2 * input_phi(i, j, k)
                end do
            end do
        end do
    end subroutine apply_A

end program cg_solver
