# -*- coding: utf-8 -*-
"""
Created on Tue Jan 16 14:33:08 2024

@author: Yang
"""

import numpy as np
import matplotlib.pyplot as plt
import time


# Domaine de calcul
Lx, Ly = 2., 1.
Nx, Ny = 4, 4  

hx, hy = Lx / Nx, Ly / Ny

# Maillage face et centre en x
Xf = np.linspace(0, Lx, Nx, endpoint=False)
Xc = np.linspace(hx / 2., Lx - hx / 2., Nx)

# Maillage face et centre en y
Yf = np.linspace(0, Ly, Ny, endpoint=False)
Yc = np.linspace(hy / 2., Ly - hy / 2., Ny)

# Déclaration des variables comme champs de données indicés comme (i, j)
size = (Nx, Ny)
u, v, pres, fi = np.zeros(size), np.zeros(size), np.zeros(size), np.zeros(size)
u_star, v_star = np.zeros(size), np.zeros(size)
NLu, NLv, Lu, Lv = np.zeros(size), np.zeros(size), np.zeros(size), np.zeros(size)
sfi = np.zeros(size)

# Déclaration des variables comme vecteur (p=i+j*nx) pour résoudre la correction de pression
b, x = np.zeros(Nx * Ny), np.zeros(Nx * Ny)

# Pour traiter les conditions aux limites périodiques
II = np.arange(Nx)
IP = np.arange(Nx) + 1
IM = np.arange(Nx) - 1
IP[Nx - 1], IM[0] = 0, Nx - 1

JJ = np.arange(Ny)
JP = np.arange(Ny) + 1
JM = np.arange(Ny) - 1
JP[Ny - 1], JM[0] = 0, Ny - 1
        
        
def periodic_index(i, I):
    if i < 0:
        return I + i
    elif i >= I:
        return i - I
    else:
        return i

def init_field():
    u, v = np.zeros(size), np.zeros(size)
    Pj, Rj, Ax = 20, 0.25, 0.5
    lbd_x = 0.5 * Lx
    for i in II:
        for j in JJ:
            u[i, j] = 0.5 * (1 + np.tanh(0.5 * Pj * (1 - np.abs(Ly / 2. - Yc[j]) / Rj)))
            u[i, j] = u[i, j] * (1 + Ax * np.sin(2 * np.pi * Xf[i] / lbd_x))
    return [u, v]

# Calcul des termes non-linéaires NLu & NLv
def get_NonLinear_terms(u, v):
    NLu, NLv = np.zeros(np.shape(u)), np.zeros(np.shape(v))
    for i in II:
        for j in JJ:
            im, ip, jm, jp = IM[i], IP[i], JM[j], JP[j]
            ul = 0.5 * (u[im, j] + u[i, j])
            ur = 0.5 * (u[ip, j] + u[i, j])

            ub = 0.5 * (u[i, jm] + u[i, j])
            ut = 0.5 * (u[i, jp] + u[i, j])

            up = u[i, j]
            vp = 0.25 * (v[i, j] + v[i, jp] + v[im, jp] + v[im, j])

            NLu[i, j] = up * (ur - ul) / hx + vp * (ut - ub) / hy

            vb = 0.5 * (v[i, jm] + v[i, j])
            vt = 0.5 * (v[i, jp] + v[i, j])

            vl = 0.5 * (v[im, j] + v[i, j])
            vr = 0.5 * (v[ip, j] + v[i, j])

            vp = v[i, j]
            up = 0.25 * (u[i, j] + u[ip, j] + u[i, jm] + u[ip, jm])

            NLv[i, j] = up * (vr - vl) / hx + vp * (vt - vb) / hy

    return [NLu, NLv]

# Calcul des termes non-linéaires Lu & Lv
def get_Linear_terms(u, v):
    Lu, Lv = np.zeros(np.shape(u)), np.zeros(np.shape(v))
    for i in II:
        for j in JJ:
            im, ip, jm, jp = IM[i], IP[i], JM[j], JP[j]
            Lu[i, j] = (u[ip, j] - 2 * u[i, j] + u[im, j]) / hx**2
            Lu[i, j] = Lu[i, j] + (u[i, jp] - 2 * u[i, j] + u[i, jm]) / hy**2

            Lv[i, j] = (v[ip, j] - 2 * v[i, j] + v[im, j]) / hx**2
            Lv[i, j] = Lv[i, j] + (v[i, jp] - 2 * v[i, j] + v[i, jm]) / hy**2

    return [Lu, Lv]

# Calcul de la divergence
def get_divergence(u, v):
    Dx, Dy = np.zeros(np.shape(u)), np.zeros(np.shape(v))
    for i in II:
        for j in JJ:
            im, ip, jm, jp = IM[i], IP[i], JM[j], JP[j]
            Dx[i, j] = (u[ip, j] - u[i, j]) / hx
            Dy[i, j] = (v[i, jp] - v[i, j]) / hy
    div = Dx + Dy
    return div

def get_gradient(fi):
    Gx, Gy = np.zeros(np.shape(fi)), np.zeros(np.shape(fi))
    for i in II:
        for j in JJ:
            im, ip, jm, jp = IM[i], IP[i], JM[j], JP[j]
            Gx[i, j] = (fi[i, j] - fi[im, j]) / hx
            Gy[i, j] = (fi[i, j] - fi[i, jm]) / hy
    return [Gx, Gy]

def index(i, j):
    p = i + j * Nx
    return p

def BuildPoisson(hx, hy, nx, ny):
    M = np.zeros((nx * ny, nx * ny))
    for i in range(nx):
        for j in range(ny):
            im, ip, ii = IM[i], IP[i], II[i]
            jm, jp, jj = JM[j], JP[j], JJ[j]
            row = index(ii, jj)
            M[row, index(ii, jj)] = -2 / hx**2 - 2 / hy**2
            M[row, index(im, jj)] = 1 / hx**2
            M[row, index(ip, jj)] = 1 / hx**2
            M[row, index(ii, jm)] = 1 / hy**2
            M[row, index(ii, jp)] = 1 / hy**2

    M[0, :] = M[0, :] * 0
    M[0, 0] = 1

    return M

# En précalcul, j'assemble le système linéaire M
M = BuildPoisson(hx, hy, Nx, Ny)
# En précalcul, on calcule l'inverse de M
Minv = np.linalg.inv(M)

nu = 1e-3
dt = 0.25 * hx**2  # Stabilité

xc, yc = np.meshgrid(Xc, Yc, indexing='ij')
xu, yu = np.meshgrid(Xf, Yc, indexing='ij')
xv, yv = np.meshgrid(Xc, Yf, indexing='ij')

u, v = init_field()

t, Tmax = 0, 1.055
plot_interval = 0.1  # Intervalle pour tracer les résultats

temps1 = [0]

while t < Tmax:
    
    start_time = time.perf_counter()
    
    t = t + dt
    NLu, NLv = get_NonLinear_terms(u, v)
    Lu, Lv = get_Linear_terms(u, v)
    dpdx, dpdy = get_gradient(pres)

    u_star = u + dt * (-dpdx - NLu + nu * Lu)
    v_star = v + dt * (-dpdy - NLv + nu * Lv)

    div = get_divergence(u_star, v_star)
    div = +div / dt

    for i in II:
        for j in JJ:
            row = index(i, j)
            b[row] = div[i, j]

    
    b[0] = 0
    x = Minv @ b
    

    for i in II:
        for j in JJ:
            row = index(i, j)
            fi[i, j] = x[row]

    Gx, Gy = get_gradient(fi)
    u = u_star - dt * Gx
    v = v_star - dt * Gy
    pres = pres + fi
    div = get_divergence(u, v)
    
    #print("%f %e" % (t, np.max(np.abs(div)))) 

    # Tracer les résultats à intervalles réguliers
    if t % plot_interval < dt:
        omega = np.copy(u)
        for i in II:
            for j in JJ:
                omega[i, j] = (v[i, j] - v[IM[i], j]) / hx - (u[i, j] - u[i, JM[j]]) / hy
        

        # Afficher les résultats
        counter = 0
        
        plt.figure(figsize=(12, 10))

        plt.subplot(2, 2, 1)
        plt.contourf(Xc, Yc, u, cmap='viridis', levels=20)
        plt.colorbar()
        plt.title('Field of u')

        plt.subplot(2, 2, 2)
        plt.contourf(Xc, Yc, v, cmap='viridis', levels=20)
        plt.colorbar()
        plt.title('Field of v')

        plt.subplot(2, 2, 3)
        plt.contourf(Xc, Yc, pres, cmap='viridis', levels=20)
        plt.colorbar()
        plt.title('Pressure Field')

        plt.subplot(2, 2, 4)
        plt.contourf(Xc, Yc, omega, cmap='viridis', levels=20)
        plt.colorbar()
        plt.title('Vorticity Field')

        plt.tight_layout()
        plt.show()
        
        counter += 1
    end_time = time.perf_counter()
    elapsed_time = end_time - start_time
    temps1.append(temps1[-1] + elapsed_time)  
    
execution_time1 = np.zeros(10)
execution_time1[0] = temps1[0]
execution_time1[-1] = temps1[-1]
delta = round(len(temps1)/10)
for i in range(1,9):
    execution_time1[i] = temps1[i*delta]        
print(execution_time1) 





#%% Domaine de calcul
Lx, Ly = 2., 1.
Nx, Ny = 18, 18  

hx, hy = Lx / Nx, Ly / Ny

# Maillage face et centre en x
Xf = np.linspace(0, Lx, Nx, endpoint=False)
Xc = np.linspace(hx / 2., Lx - hx / 2., Nx)

# Maillage face et centre en y
Yf = np.linspace(0, Ly, Ny, endpoint=False)
Yc = np.linspace(hy / 2., Ly - hy / 2., Ny)

# Déclaration des variables comme champs de données indicés comme (i, j)
size = (Nx, Ny)
u, v, pres, fi = np.zeros(size), np.zeros(size), np.zeros(size), np.zeros(size)
u_star, v_star = np.zeros(size), np.zeros(size)
NLu, NLv, Lu, Lv = np.zeros(size), np.zeros(size), np.zeros(size), np.zeros(size)
sfi = np.zeros(size)

# Déclaration des variables comme vecteur (p=i+j*nx) pour résoudre la correction de pression
b, x = np.zeros(Nx * Ny), np.zeros(Nx * Ny)

# Pour traiter les conditions aux limites périodiques
II = np.arange(Nx)
IP = np.arange(Nx) + 1
IM = np.arange(Nx) - 1
IP[Nx - 1], IM[0] = 0, Nx - 1

JJ = np.arange(Ny)
JP = np.arange(Ny) + 1
JM = np.arange(Ny) - 1
JP[Ny - 1], JM[0] = 0, Ny - 1
        
        
def periodic_index(i, I):
    if i < 0:
        return I + i
    elif i >= I:
        return i - I
    else:
        return i

def init_field():
    u, v = np.zeros(size), np.zeros(size)
    Pj, Rj, Ax = 20, 0.25, 0.5
    lbd_x = 0.5 * Lx
    for i in II:
        for j in JJ:
            u[i, j] = 0.5 * (1 + np.tanh(0.5 * Pj * (1 - np.abs(Ly / 2. - Yc[j]) / Rj)))
            u[i, j] = u[i, j] * (1 + Ax * np.sin(2 * np.pi * Xf[i] / lbd_x))
    return [u, v]

# Calcul des termes non-linéaires NLu & NLv
def get_NonLinear_terms(u, v):
    NLu, NLv = np.zeros(np.shape(u)), np.zeros(np.shape(v))
    for i in II:
        for j in JJ:
            im, ip, jm, jp = IM[i], IP[i], JM[j], JP[j]
            ul = 0.5 * (u[im, j] + u[i, j])
            ur = 0.5 * (u[ip, j] + u[i, j])

            ub = 0.5 * (u[i, jm] + u[i, j])
            ut = 0.5 * (u[i, jp] + u[i, j])

            up = u[i, j]
            vp = 0.25 * (v[i, j] + v[i, jp] + v[im, jp] + v[im, j])

            NLu[i, j] = up * (ur - ul) / hx + vp * (ut - ub) / hy

            vb = 0.5 * (v[i, jm] + v[i, j])
            vt = 0.5 * (v[i, jp] + v[i, j])

            vl = 0.5 * (v[im, j] + v[i, j])
            vr = 0.5 * (v[ip, j] + v[i, j])

            vp = v[i, j]
            up = 0.25 * (u[i, j] + u[ip, j] + u[i, jm] + u[ip, jm])

            NLv[i, j] = up * (vr - vl) / hx + vp * (vt - vb) / hy

    return [NLu, NLv]

# Calcul des termes non-linéaires Lu & Lv
def get_Linear_terms(u, v):
    Lu, Lv = np.zeros(np.shape(u)), np.zeros(np.shape(v))
    for i in II:
        for j in JJ:
            im, ip, jm, jp = IM[i], IP[i], JM[j], JP[j]
            Lu[i, j] = (u[ip, j] - 2 * u[i, j] + u[im, j]) / hx**2
            Lu[i, j] = Lu[i, j] + (u[i, jp] - 2 * u[i, j] + u[i, jm]) / hy**2

            Lv[i, j] = (v[ip, j] - 2 * v[i, j] + v[im, j]) / hx**2
            Lv[i, j] = Lv[i, j] + (v[i, jp] - 2 * v[i, j] + v[i, jm]) / hy**2

    return [Lu, Lv]

# Calcul de la divergence
def get_divergence(u, v):
    Dx, Dy = np.zeros(np.shape(u)), np.zeros(np.shape(v))
    for i in II:
        for j in JJ:
            im, ip, jm, jp = IM[i], IP[i], JM[j], JP[j]
            Dx[i, j] = (u[ip, j] - u[i, j]) / hx
            Dy[i, j] = (v[i, jp] - v[i, j]) / hy
    div = Dx + Dy
    return div

def get_gradient(fi):
    Gx, Gy = np.zeros(np.shape(fi)), np.zeros(np.shape(fi))
    for i in II:
        for j in JJ:
            im, ip, jm, jp = IM[i], IP[i], JM[j], JP[j]
            Gx[i, j] = (fi[i, j] - fi[im, j]) / hx
            Gy[i, j] = (fi[i, j] - fi[i, jm]) / hy
    return [Gx, Gy]

def index(i, j):
    p = i + j * Nx
    return p

def BuildPoisson(hx, hy, nx, ny):
    M = np.zeros((nx * ny, nx * ny))
    for i in range(nx):
        for j in range(ny):
            im, ip, ii = IM[i], IP[i], II[i]
            jm, jp, jj = JM[j], JP[j], JJ[j]
            row = index(ii, jj)
            M[row, index(ii, jj)] = -2 / hx**2 - 2 / hy**2
            M[row, index(im, jj)] = 1 / hx**2
            M[row, index(ip, jj)] = 1 / hx**2
            M[row, index(ii, jm)] = 1 / hy**2
            M[row, index(ii, jp)] = 1 / hy**2

    M[0, :] = M[0, :] * 0
    M[0, 0] = 1

    return M

# En précalcul, j'assemble le système linéaire M
M = BuildPoisson(hx, hy, Nx, Ny)
# En précalcul, on calcule l'inverse de M
Minv = np.linalg.inv(M)

nu = 1e-3
dt = 0.25 * hx**2  # Stabilité

xc, yc = np.meshgrid(Xc, Yc, indexing='ij')
xu, yu = np.meshgrid(Xf, Yc, indexing='ij')
xv, yv = np.meshgrid(Xc, Yf, indexing='ij')

u, v = init_field()

t, Tmax = 0, 1.055
plot_interval = 0.1  # Intervalle pour tracer les résultats

temps2 = [execution_time1[-1]]


while t < Tmax:
    
    start_time = time.perf_counter()
    
    t = t + dt
    NLu, NLv = get_NonLinear_terms(u, v)
    Lu, Lv = get_Linear_terms(u, v)
    dpdx, dpdy = get_gradient(pres)

    u_star = u + dt * (-dpdx - NLu + nu * Lu)
    v_star = v + dt * (-dpdy - NLv + nu * Lv)

    div = get_divergence(u_star, v_star)
    div = +div / dt

    for i in II:
        for j in JJ:
            row = index(i, j)
            b[row] = div[i, j]

    b[0] = 0
    x = Minv @ b

    for i in II:
        for j in JJ:
            row = index(i, j)
            fi[i, j] = x[row]

    Gx, Gy = get_gradient(fi)
    u = u_star - dt * Gx
    v = v_star - dt * Gy
    pres = pres + fi
    div = get_divergence(u, v)
    
    #print("%f %e" % (t, np.max(np.abs(div)))) 

    # Tracer les résultats à intervalles réguliers
    if t % plot_interval < dt:
        omega = np.copy(u)
        for i in II:
            for j in JJ:
                omega[i, j] = (v[i, j] - v[IM[i], j]) / hx - (u[i, j] - u[i, JM[j]]) / hy
        

        # Afficher les résultats
        
        counter = counter
        
        plt.figure(figsize=(12, 10))

        plt.subplot(2, 2, 1)
        plt.contourf(Xc, Yc, u, cmap='viridis', levels=20)
        plt.colorbar()
        plt.title('Field of u')

        plt.subplot(2, 2, 2)
        plt.contourf(Xc, Yc, v, cmap='viridis', levels=20)
        plt.colorbar()
        plt.title('Field of v')

        plt.subplot(2, 2, 3)
        plt.contourf(Xc, Yc, pres, cmap='viridis', levels=20)
        plt.colorbar()
        plt.title('Pressure Field')

        plt.subplot(2, 2, 4)
        plt.contourf(Xc, Yc, omega, cmap='viridis', levels=20)
        plt.colorbar()
        plt.title('Vorticity Field')

        plt.tight_layout()
        plt.show()
        
        counter += 1
    
    end_time = time.perf_counter()
    elapsed_time = end_time - start_time
    temps2.append(temps2[-1] + elapsed_time)
execution_time2 = np.zeros(10)
execution_time2[0] = temps2[0]
execution_time2[-1] = temps2[-1]
delta2 = round(len(temps2)/10)
print(delta2)
for i in range(1,9):
    execution_time2[i] = temps2[i*delta2] 

         
print(execution_time2)


#%% Domaine de calcul
Lx, Ly = 2., 1.
Nx, Ny = 32, 32  

hx, hy = Lx / Nx, Ly / Ny

# Maillage face et centre en x
Xf = np.linspace(0, Lx, Nx, endpoint=False)
Xc = np.linspace(hx / 2., Lx - hx / 2., Nx)

# Maillage face et centre en y
Yf = np.linspace(0, Ly, Ny, endpoint=False)
Yc = np.linspace(hy / 2., Ly - hy / 2., Ny)

# Déclaration des variables comme champs de données indicés comme (i, j)
size = (Nx, Ny)
u, v, pres, fi = np.zeros(size), np.zeros(size), np.zeros(size), np.zeros(size)
u_star, v_star = np.zeros(size), np.zeros(size)
NLu, NLv, Lu, Lv = np.zeros(size), np.zeros(size), np.zeros(size), np.zeros(size)
sfi = np.zeros(size)

# Déclaration des variables comme vecteur (p=i+j*nx) pour résoudre la correction de pression
b, x = np.zeros(Nx * Ny), np.zeros(Nx * Ny)

# Pour traiter les conditions aux limites périodiques
II = np.arange(Nx)
IP = np.arange(Nx) + 1
IM = np.arange(Nx) - 1
IP[Nx - 1], IM[0] = 0, Nx - 1

JJ = np.arange(Ny)
JP = np.arange(Ny) + 1
JM = np.arange(Ny) - 1
JP[Ny - 1], JM[0] = 0, Ny - 1
        
        
def periodic_index(i, I):
    if i < 0:
        return I + i
    elif i >= I:
        return i - I
    else:
        return i

def init_field():
    u, v = np.zeros(size), np.zeros(size)
    Pj, Rj, Ax = 20, 0.25, 0.5
    lbd_x = 0.5 * Lx
    for i in II:
        for j in JJ:
            u[i, j] = 0.5 * (1 + np.tanh(0.5 * Pj * (1 - np.abs(Ly / 2. - Yc[j]) / Rj)))
            u[i, j] = u[i, j] * (1 + Ax * np.sin(2 * np.pi * Xf[i] / lbd_x))
    return [u, v]

# Calcul des termes non-linéaires NLu & NLv
def get_NonLinear_terms(u, v):
    NLu, NLv = np.zeros(np.shape(u)), np.zeros(np.shape(v))
    for i in II:
        for j in JJ:
            im, ip, jm, jp = IM[i], IP[i], JM[j], JP[j]
            ul = 0.5 * (u[im, j] + u[i, j])
            ur = 0.5 * (u[ip, j] + u[i, j])

            ub = 0.5 * (u[i, jm] + u[i, j])
            ut = 0.5 * (u[i, jp] + u[i, j])

            up = u[i, j]
            vp = 0.25 * (v[i, j] + v[i, jp] + v[im, jp] + v[im, j])

            NLu[i, j] = up * (ur - ul) / hx + vp * (ut - ub) / hy

            vb = 0.5 * (v[i, jm] + v[i, j])
            vt = 0.5 * (v[i, jp] + v[i, j])

            vl = 0.5 * (v[im, j] + v[i, j])
            vr = 0.5 * (v[ip, j] + v[i, j])

            vp = v[i, j]
            up = 0.25 * (u[i, j] + u[ip, j] + u[i, jm] + u[ip, jm])

            NLv[i, j] = up * (vr - vl) / hx + vp * (vt - vb) / hy

    return [NLu, NLv]

# Calcul des termes non-linéaires Lu & Lv
def get_Linear_terms(u, v):
    Lu, Lv = np.zeros(np.shape(u)), np.zeros(np.shape(v))
    for i in II:
        for j in JJ:
            im, ip, jm, jp = IM[i], IP[i], JM[j], JP[j]
            Lu[i, j] = (u[ip, j] - 2 * u[i, j] + u[im, j]) / hx**2
            Lu[i, j] = Lu[i, j] + (u[i, jp] - 2 * u[i, j] + u[i, jm]) / hy**2

            Lv[i, j] = (v[ip, j] - 2 * v[i, j] + v[im, j]) / hx**2
            Lv[i, j] = Lv[i, j] + (v[i, jp] - 2 * v[i, j] + v[i, jm]) / hy**2

    return [Lu, Lv]

# Calcul de la divergence
def get_divergence(u, v):
    Dx, Dy = np.zeros(np.shape(u)), np.zeros(np.shape(v))
    for i in II:
        for j in JJ:
            im, ip, jm, jp = IM[i], IP[i], JM[j], JP[j]
            Dx[i, j] = (u[ip, j] - u[i, j]) / hx
            Dy[i, j] = (v[i, jp] - v[i, j]) / hy
    div = Dx + Dy
    return div

def get_gradient(fi):
    Gx, Gy = np.zeros(np.shape(fi)), np.zeros(np.shape(fi))
    for i in II:
        for j in JJ:
            im, ip, jm, jp = IM[i], IP[i], JM[j], JP[j]
            Gx[i, j] = (fi[i, j] - fi[im, j]) / hx
            Gy[i, j] = (fi[i, j] - fi[i, jm]) / hy
    return [Gx, Gy]

def index(i, j):
    p = i + j * Nx
    return p

def BuildPoisson(hx, hy, nx, ny):
    M = np.zeros((nx * ny, nx * ny))
    for i in range(nx):
        for j in range(ny):
            im, ip, ii = IM[i], IP[i], II[i]
            jm, jp, jj = JM[j], JP[j], JJ[j]
            row = index(ii, jj)
            M[row, index(ii, jj)] = -2 / hx**2 - 2 / hy**2
            M[row, index(im, jj)] = 1 / hx**2
            M[row, index(ip, jj)] = 1 / hx**2
            M[row, index(ii, jm)] = 1 / hy**2
            M[row, index(ii, jp)] = 1 / hy**2

    M[0, :] = M[0, :] * 0
    M[0, 0] = 1

    return M

# En précalcul, j'assemble le système linéaire M
M = BuildPoisson(hx, hy, Nx, Ny)
# En précalcul, on calcule l'inverse de M
Minv = np.linalg.inv(M)

nu = 1e-3
dt = 0.25 * hx**2  # Stabilité

xc, yc = np.meshgrid(Xc, Yc, indexing='ij')
xu, yu = np.meshgrid(Xf, Yc, indexing='ij')
xv, yv = np.meshgrid(Xc, Yf, indexing='ij')

u, v = init_field()

t, Tmax = 0, 1.055
plot_interval = 0.1  # Intervalle pour tracer les résultats

temps3 = [execution_time2[-1]]


while t < Tmax:
    
    start_time = time.perf_counter()
    
    t = t + dt
    NLu, NLv = get_NonLinear_terms(u, v)
    Lu, Lv = get_Linear_terms(u, v)
    dpdx, dpdy = get_gradient(pres)

    u_star = u + dt * (-dpdx - NLu + nu * Lu)
    v_star = v + dt * (-dpdy - NLv + nu * Lv)

    div = get_divergence(u_star, v_star)
    div = +div / dt

    for i in II:
        for j in JJ:
            row = index(i, j)
            b[row] = div[i, j]

    b[0] = 0
    x = Minv @ b

    for i in II:
        for j in JJ:
            row = index(i, j)
            fi[i, j] = x[row]

    Gx, Gy = get_gradient(fi)
    u = u_star - dt * Gx
    v = v_star - dt * Gy
    pres = pres + fi
    div = get_divergence(u, v)
    
    #print("%f %e" % (t, np.max(np.abs(div)))) 

    # Tracer les résultats à intervalles réguliers
    if t % plot_interval < dt:
        omega = np.copy(u)
        for i in II:
            for j in JJ:
                omega[i, j] = (v[i, j] - v[IM[i], j]) / hx - (u[i, j] - u[i, JM[j]]) / hy
        

        # Afficher les résultats
        
        counter = counter
        
        plt.figure(figsize=(12, 10))

        plt.subplot(2, 2, 1)
        plt.contourf(Xc, Yc, u, cmap='viridis', levels=20)
        plt.colorbar()
        plt.title('Field of u')

        plt.subplot(2, 2, 2)
        plt.contourf(Xc, Yc, v, cmap='viridis', levels=20)
        plt.colorbar()
        plt.title('Field of v')

        plt.subplot(2, 2, 3)
        plt.contourf(Xc, Yc, pres, cmap='viridis', levels=20)
        plt.colorbar()
        plt.title('Pressure Field')

        plt.subplot(2, 2, 4)
        plt.contourf(Xc, Yc, omega, cmap='viridis', levels=20)
        plt.colorbar()
        plt.title('Vorticity Field')

        plt.tight_layout()
        plt.show()
        
        counter += 1
    
    end_time = time.perf_counter()
    elapsed_time = end_time - start_time
    temps3.append(temps3[-1] + elapsed_time)
execution_time3 = np.zeros(10)
execution_time3[0] = temps3[0]
execution_time3[-1] = temps3[-1]
delta3 = round(len(temps3)/10)
print(delta3)
for i in range(1,9):
    execution_time3[i] = temps3[i*delta3] 

         
print(execution_time3)


#%% Domaine de calcul
Lx, Ly = 2., 1.
Nx, Ny = 46, 46  

hx, hy = Lx / Nx, Ly / Ny

# Maillage face et centre en x
Xf = np.linspace(0, Lx, Nx, endpoint=False)
Xc = np.linspace(hx / 2., Lx - hx / 2., Nx)

# Maillage face et centre en y
Yf = np.linspace(0, Ly, Ny, endpoint=False)
Yc = np.linspace(hy / 2., Ly - hy / 2., Ny)

# Déclaration des variables comme champs de données indicés comme (i, j)
size = (Nx, Ny)
u, v, pres, fi = np.zeros(size), np.zeros(size), np.zeros(size), np.zeros(size)
u_star, v_star = np.zeros(size), np.zeros(size)
NLu, NLv, Lu, Lv = np.zeros(size), np.zeros(size), np.zeros(size), np.zeros(size)
sfi = np.zeros(size)

# Déclaration des variables comme vecteur (p=i+j*nx) pour résoudre la correction de pression
b, x = np.zeros(Nx * Ny), np.zeros(Nx * Ny)

# Pour traiter les conditions aux limites périodiques
II = np.arange(Nx)
IP = np.arange(Nx) + 1
IM = np.arange(Nx) - 1
IP[Nx - 1], IM[0] = 0, Nx - 1

JJ = np.arange(Ny)
JP = np.arange(Ny) + 1
JM = np.arange(Ny) - 1
JP[Ny - 1], JM[0] = 0, Ny - 1
        
        
def periodic_index(i, I):
    if i < 0:
        return I + i
    elif i >= I:
        return i - I
    else:
        return i

def init_field():
    u, v = np.zeros(size), np.zeros(size)
    Pj, Rj, Ax = 20, 0.25, 0.5
    lbd_x = 0.5 * Lx
    for i in II:
        for j in JJ:
            u[i, j] = 0.5 * (1 + np.tanh(0.5 * Pj * (1 - np.abs(Ly / 2. - Yc[j]) / Rj)))
            u[i, j] = u[i, j] * (1 + Ax * np.sin(2 * np.pi * Xf[i] / lbd_x))
    return [u, v]

# Calcul des termes non-linéaires NLu & NLv
def get_NonLinear_terms(u, v):
    NLu, NLv = np.zeros(np.shape(u)), np.zeros(np.shape(v))
    for i in II:
        for j in JJ:
            im, ip, jm, jp = IM[i], IP[i], JM[j], JP[j]
            ul = 0.5 * (u[im, j] + u[i, j])
            ur = 0.5 * (u[ip, j] + u[i, j])

            ub = 0.5 * (u[i, jm] + u[i, j])
            ut = 0.5 * (u[i, jp] + u[i, j])

            up = u[i, j]
            vp = 0.25 * (v[i, j] + v[i, jp] + v[im, jp] + v[im, j])

            NLu[i, j] = up * (ur - ul) / hx + vp * (ut - ub) / hy

            vb = 0.5 * (v[i, jm] + v[i, j])
            vt = 0.5 * (v[i, jp] + v[i, j])

            vl = 0.5 * (v[im, j] + v[i, j])
            vr = 0.5 * (v[ip, j] + v[i, j])

            vp = v[i, j]
            up = 0.25 * (u[i, j] + u[ip, j] + u[i, jm] + u[ip, jm])

            NLv[i, j] = up * (vr - vl) / hx + vp * (vt - vb) / hy

    return [NLu, NLv]

# Calcul des termes non-linéaires Lu & Lv
def get_Linear_terms(u, v):
    Lu, Lv = np.zeros(np.shape(u)), np.zeros(np.shape(v))
    for i in II:
        for j in JJ:
            im, ip, jm, jp = IM[i], IP[i], JM[j], JP[j]
            Lu[i, j] = (u[ip, j] - 2 * u[i, j] + u[im, j]) / hx**2
            Lu[i, j] = Lu[i, j] + (u[i, jp] - 2 * u[i, j] + u[i, jm]) / hy**2

            Lv[i, j] = (v[ip, j] - 2 * v[i, j] + v[im, j]) / hx**2
            Lv[i, j] = Lv[i, j] + (v[i, jp] - 2 * v[i, j] + v[i, jm]) / hy**2

    return [Lu, Lv]

# Calcul de la divergence
def get_divergence(u, v):
    Dx, Dy = np.zeros(np.shape(u)), np.zeros(np.shape(v))
    for i in II:
        for j in JJ:
            im, ip, jm, jp = IM[i], IP[i], JM[j], JP[j]
            Dx[i, j] = (u[ip, j] - u[i, j]) / hx
            Dy[i, j] = (v[i, jp] - v[i, j]) / hy
    div = Dx + Dy
    return div

def get_gradient(fi):
    Gx, Gy = np.zeros(np.shape(fi)), np.zeros(np.shape(fi))
    for i in II:
        for j in JJ:
            im, ip, jm, jp = IM[i], IP[i], JM[j], JP[j]
            Gx[i, j] = (fi[i, j] - fi[im, j]) / hx
            Gy[i, j] = (fi[i, j] - fi[i, jm]) / hy
    return [Gx, Gy]

def index(i, j):
    p = i + j * Nx
    return p

def BuildPoisson(hx, hy, nx, ny):
    M = np.zeros((nx * ny, nx * ny))
    for i in range(nx):
        for j in range(ny):
            im, ip, ii = IM[i], IP[i], II[i]
            jm, jp, jj = JM[j], JP[j], JJ[j]
            row = index(ii, jj)
            M[row, index(ii, jj)] = -2 / hx**2 - 2 / hy**2
            M[row, index(im, jj)] = 1 / hx**2
            M[row, index(ip, jj)] = 1 / hx**2
            M[row, index(ii, jm)] = 1 / hy**2
            M[row, index(ii, jp)] = 1 / hy**2

    M[0, :] = M[0, :] * 0
    M[0, 0] = 1

    return M

# En précalcul, j'assemble le système linéaire M
M = BuildPoisson(hx, hy, Nx, Ny)
# En précalcul, on calcule l'inverse de M
Minv = np.linalg.inv(M)

nu = 1e-3
dt = 0.25 * hx**2  # Stabilité

xc, yc = np.meshgrid(Xc, Yc, indexing='ij')
xu, yu = np.meshgrid(Xf, Yc, indexing='ij')
xv, yv = np.meshgrid(Xc, Yf, indexing='ij')

u, v = init_field()

t, Tmax = 0, 1.055
plot_interval = 0.1  # Intervalle pour tracer les résultats

temps4 = [execution_time3[-1]]


while t < Tmax:
    
    start_time = time.perf_counter()
    
    t = t + dt
    NLu, NLv = get_NonLinear_terms(u, v)
    Lu, Lv = get_Linear_terms(u, v)
    dpdx, dpdy = get_gradient(pres)

    u_star = u + dt * (-dpdx - NLu + nu * Lu)
    v_star = v + dt * (-dpdy - NLv + nu * Lv)

    div = get_divergence(u_star, v_star)
    div = +div / dt

    for i in II:
        for j in JJ:
            row = index(i, j)
            b[row] = div[i, j]

    b[0] = 0
    x = Minv @ b

    for i in II:
        for j in JJ:
            row = index(i, j)
            fi[i, j] = x[row]

    Gx, Gy = get_gradient(fi)
    u = u_star - dt * Gx
    v = v_star - dt * Gy
    pres = pres + fi
    div = get_divergence(u, v)
    
    #print("%f %e" % (t, np.max(np.abs(div)))) 

    # Tracer les résultats à intervalles réguliers
    if t % plot_interval < dt:
        omega = np.copy(u)
        for i in II:
            for j in JJ:
                omega[i, j] = (v[i, j] - v[IM[i], j]) / hx - (u[i, j] - u[i, JM[j]]) / hy
        

        # Afficher les résultats
        
        counter = counter
        
        plt.figure(figsize=(12, 10))

        plt.subplot(2, 2, 1)
        plt.contourf(Xc, Yc, u, cmap='viridis', levels=20)
        plt.colorbar()
        plt.title('Field of u')

        plt.subplot(2, 2, 2)
        plt.contourf(Xc, Yc, v, cmap='viridis', levels=20)
        plt.colorbar()
        plt.title('Field of v')

        plt.subplot(2, 2, 3)
        plt.contourf(Xc, Yc, pres, cmap='viridis', levels=20)
        plt.colorbar()
        plt.title('Pressure Field')

        plt.subplot(2, 2, 4)
        plt.contourf(Xc, Yc, omega, cmap='viridis', levels=20)
        plt.colorbar()
        plt.title('Vorticity Field')

        plt.tight_layout()
        plt.show()
        
        counter += 1
    
    end_time = time.perf_counter()
    elapsed_time = end_time - start_time
    temps4.append(temps4[-1] + elapsed_time)
execution_time4 = np.zeros(10)
execution_time4[0] = temps4[0]
execution_time4[-1] = temps4[-1]
delta4 = round(len(temps4)/10)
print(delta4)
for i in range(1,9):
    execution_time4[i] = temps4[i*delta4] 

         
print(execution_time4)


#%% Domaine de calcul
Lx, Ly = 2., 1.
Nx, Ny = 60, 60  

hx, hy = Lx / Nx, Ly / Ny

# Maillage face et centre en x
Xf = np.linspace(0, Lx, Nx, endpoint=False)
Xc = np.linspace(hx / 2., Lx - hx / 2., Nx)

# Maillage face et centre en y
Yf = np.linspace(0, Ly, Ny, endpoint=False)
Yc = np.linspace(hy / 2., Ly - hy / 2., Ny)

# Déclaration des variables comme champs de données indicés comme (i, j)
size = (Nx, Ny)
u, v, pres, fi = np.zeros(size), np.zeros(size), np.zeros(size), np.zeros(size)
u_star, v_star = np.zeros(size), np.zeros(size)
NLu, NLv, Lu, Lv = np.zeros(size), np.zeros(size), np.zeros(size), np.zeros(size)
sfi = np.zeros(size)

# Déclaration des variables comme vecteur (p=i+j*nx) pour résoudre la correction de pression
b, x = np.zeros(Nx * Ny), np.zeros(Nx * Ny)

# Pour traiter les conditions aux limites périodiques
II = np.arange(Nx)
IP = np.arange(Nx) + 1
IM = np.arange(Nx) - 1
IP[Nx - 1], IM[0] = 0, Nx - 1

JJ = np.arange(Ny)
JP = np.arange(Ny) + 1
JM = np.arange(Ny) - 1
JP[Ny - 1], JM[0] = 0, Ny - 1
        
        
def periodic_index(i, I):
    if i < 0:
        return I + i
    elif i >= I:
        return i - I
    else:
        return i

def init_field():
    u, v = np.zeros(size), np.zeros(size)
    Pj, Rj, Ax = 20, 0.25, 0.5
    lbd_x = 0.5 * Lx
    for i in II:
        for j in JJ:
            u[i, j] = 0.5 * (1 + np.tanh(0.5 * Pj * (1 - np.abs(Ly / 2. - Yc[j]) / Rj)))
            u[i, j] = u[i, j] * (1 + Ax * np.sin(2 * np.pi * Xf[i] / lbd_x))
    return [u, v]

# Calcul des termes non-linéaires NLu & NLv
def get_NonLinear_terms(u, v):
    NLu, NLv = np.zeros(np.shape(u)), np.zeros(np.shape(v))
    for i in II:
        for j in JJ:
            im, ip, jm, jp = IM[i], IP[i], JM[j], JP[j]
            ul = 0.5 * (u[im, j] + u[i, j])
            ur = 0.5 * (u[ip, j] + u[i, j])

            ub = 0.5 * (u[i, jm] + u[i, j])
            ut = 0.5 * (u[i, jp] + u[i, j])

            up = u[i, j]
            vp = 0.25 * (v[i, j] + v[i, jp] + v[im, jp] + v[im, j])

            NLu[i, j] = up * (ur - ul) / hx + vp * (ut - ub) / hy

            vb = 0.5 * (v[i, jm] + v[i, j])
            vt = 0.5 * (v[i, jp] + v[i, j])

            vl = 0.5 * (v[im, j] + v[i, j])
            vr = 0.5 * (v[ip, j] + v[i, j])

            vp = v[i, j]
            up = 0.25 * (u[i, j] + u[ip, j] + u[i, jm] + u[ip, jm])

            NLv[i, j] = up * (vr - vl) / hx + vp * (vt - vb) / hy

    return [NLu, NLv]

# Calcul des termes non-linéaires Lu & Lv
def get_Linear_terms(u, v):
    Lu, Lv = np.zeros(np.shape(u)), np.zeros(np.shape(v))
    for i in II:
        for j in JJ:
            im, ip, jm, jp = IM[i], IP[i], JM[j], JP[j]
            Lu[i, j] = (u[ip, j] - 2 * u[i, j] + u[im, j]) / hx**2
            Lu[i, j] = Lu[i, j] + (u[i, jp] - 2 * u[i, j] + u[i, jm]) / hy**2

            Lv[i, j] = (v[ip, j] - 2 * v[i, j] + v[im, j]) / hx**2
            Lv[i, j] = Lv[i, j] + (v[i, jp] - 2 * v[i, j] + v[i, jm]) / hy**2

    return [Lu, Lv]

# Calcul de la divergence
def get_divergence(u, v):
    Dx, Dy = np.zeros(np.shape(u)), np.zeros(np.shape(v))
    for i in II:
        for j in JJ:
            im, ip, jm, jp = IM[i], IP[i], JM[j], JP[j]
            Dx[i, j] = (u[ip, j] - u[i, j]) / hx
            Dy[i, j] = (v[i, jp] - v[i, j]) / hy
    div = Dx + Dy
    return div

def get_gradient(fi):
    Gx, Gy = np.zeros(np.shape(fi)), np.zeros(np.shape(fi))
    for i in II:
        for j in JJ:
            im, ip, jm, jp = IM[i], IP[i], JM[j], JP[j]
            Gx[i, j] = (fi[i, j] - fi[im, j]) / hx
            Gy[i, j] = (fi[i, j] - fi[i, jm]) / hy
    return [Gx, Gy]

def index(i, j):
    p = i + j * Nx
    return p

def BuildPoisson(hx, hy, nx, ny):
    M = np.zeros((nx * ny, nx * ny))
    for i in range(nx):
        for j in range(ny):
            im, ip, ii = IM[i], IP[i], II[i]
            jm, jp, jj = JM[j], JP[j], JJ[j]
            row = index(ii, jj)
            M[row, index(ii, jj)] = -2 / hx**2 - 2 / hy**2
            M[row, index(im, jj)] = 1 / hx**2
            M[row, index(ip, jj)] = 1 / hx**2
            M[row, index(ii, jm)] = 1 / hy**2
            M[row, index(ii, jp)] = 1 / hy**2

    M[0, :] = M[0, :] * 0
    M[0, 0] = 1

    return M

# En précalcul, j'assemble le système linéaire M
M = BuildPoisson(hx, hy, Nx, Ny)
# En précalcul, on calcule l'inverse de M
Minv = np.linalg.inv(M)

nu = 1e-3
dt = 0.25 * hx**2  # Stabilité

xc, yc = np.meshgrid(Xc, Yc, indexing='ij')
xu, yu = np.meshgrid(Xf, Yc, indexing='ij')
xv, yv = np.meshgrid(Xc, Yf, indexing='ij')

u, v = init_field()

t, Tmax = 0, 1.055
plot_interval = 0.1  # Intervalle pour tracer les résultats

temps5 = [execution_time4[-1]]


while t < Tmax:
    
    start_time = time.perf_counter()
    
    t = t + dt
    NLu, NLv = get_NonLinear_terms(u, v)
    Lu, Lv = get_Linear_terms(u, v)
    dpdx, dpdy = get_gradient(pres)

    u_star = u + dt * (-dpdx - NLu + nu * Lu)
    v_star = v + dt * (-dpdy - NLv + nu * Lv)

    div = get_divergence(u_star, v_star)
    div = +div / dt

    for i in II:
        for j in JJ:
            row = index(i, j)
            b[row] = div[i, j]

    b[0] = 0
    x = Minv @ b

    for i in II:
        for j in JJ:
            row = index(i, j)
            fi[i, j] = x[row]

    Gx, Gy = get_gradient(fi)
    u = u_star - dt * Gx
    v = v_star - dt * Gy
    pres = pres + fi
    div = get_divergence(u, v)
    
    #print("%f %e" % (t, np.max(np.abs(div)))) 

    # Tracer les résultats à intervalles réguliers
    if t % plot_interval < dt:
        omega = np.copy(u)
        for i in II:
            for j in JJ:
                omega[i, j] = (v[i, j] - v[IM[i], j]) / hx - (u[i, j] - u[i, JM[j]]) / hy
        

        # Afficher les résultats
        
        counter = counter
        
        plt.figure(figsize=(12, 10))

        plt.subplot(2, 2, 1)
        plt.contourf(Xc, Yc, u, cmap='viridis', levels=20)
        plt.colorbar()
        plt.title('Field of u')

        plt.subplot(2, 2, 2)
        plt.contourf(Xc, Yc, v, cmap='viridis', levels=20)
        plt.colorbar()
        plt.title('Field of v')

        plt.subplot(2, 2, 3)
        plt.contourf(Xc, Yc, pres, cmap='viridis', levels=20)
        plt.colorbar()
        plt.title('Pressure Field')

        plt.subplot(2, 2, 4)
        plt.contourf(Xc, Yc, omega, cmap='viridis', levels=20)
        plt.colorbar()
        plt.title('Vorticity Field')

        plt.tight_layout()
        plt.show()
        
        counter += 1
    
    end_time = time.perf_counter()
    elapsed_time = end_time - start_time
    temps5.append(temps5[-1] + elapsed_time)
execution_time5 = np.zeros(10)
execution_time5[0] = temps5[0]
execution_time5[-1] = temps5[-1]
delta5 = round(len(temps5)/10)
print(delta5)
for i in range(1,9):
    execution_time5[i] = temps5[i*delta5] 

         
print(execution_time5)
#%%
execution_time = np.zeros(50)
execution_time[0:10] = execution_time1
execution_time[10:20] = execution_time2
execution_time[20:30] = execution_time3
execution_time[30:40] = execution_time4
execution_time[40:50] = execution_time5

plt.figure()
num_graphe = np.linspace(0.0, 50.5, 50)
plt.plot(num_graphe,execution_time)
plt.title('Temps de calcul',size=13)
plt.xlabel("Num_Graphe", size=10)
plt.ylabel("temps", size=10)
plt.legend()
plt.show()