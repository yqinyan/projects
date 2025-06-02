# -*- coding: utf-8 -*-
"""
Created on Sat Jun 29 19:54:35 2024

@author: YY
"""

import numpy as np
import matplotlib.pyplot as plt


def sol_exacte(t, x, y, omega, mu, e): 
# mu est la permittivite magnetique, e est la permitivite electrique
    
    c = 1/np.sqrt(mu*e)
    
    E = np.zeros((len(t), len(x)))
    H = np.zeros((len(t), len(y)))
    
    for n in range(len(t)):
        E[n, :] = np.sin(omega * x[:]) * np.sin(c * omega * t[n])
        
        H[n, :] = - 1/(mu*c) * np.cos(c * omega * (delta_t/2+t[n])) * np.cos(omega * y[:]) 
    return H, E



def Maxwell1d(N, I, T, l, L, omega, mu, e):
    
    c = 1/np.sqrt(mu*e)
    
    delta_t = T/N  #le pas de discritisation temporelle
    delta_x = (L-l)/I  #le pas de discritisation spatiale
    coeff = delta_t/delta_x
    
    # I+1 points de discritisation sur pour Ez sur [0,L]
    z = np.arange(l, L + 0.1 * delta_x, delta_x) 
    # I+1 points de discritisation pour Hy sur [delta_x/2,L+delta_x/2]
    y = np.arange(l + delta_x/2, L + 0.1 * delta_x, delta_x)
    # I+1 points de discritisation pour Hy sur [delta_x/2,L+delta_x/2]
    t = np.arange(0, T + 0.1 * delta_t, delta_t) 
    
    # forme des champs magnetique et eletrique
    
    Hy = np.zeros((N+1,I))
    Ez = np.zeros((N+1, I+1))
    
    # conditions initiales 
    
    
    Hy[0, :] = - 1/(mu*c) * np.cos(c * omega * (delta_t/2+t[0])) * np.cos(omega * y[:])
    Ez[0, :] = np.sin(omega * z[:]) * np.sin(c * omega * t[0])
    
    Ez[1, 1:-1] = Ez[0, 1:-1] + (coeff/e) * (Hy[0, 1:] - Hy[0, :-1])
    # conditions de bords
    Ez[1, 0] = 0
    Ez[1, -1] = 0
    
    for n in range(1,N):
        
        Hy[n, :-1] = Hy[n-1, :-1] + (coeff/mu) * (Ez[n, 1:-1] - Ez[n, :-2])
        Hy[n, -1] = Hy[n-1, -1] + (coeff/mu) * (Ez[n, -1] - Ez[n, -2])  
        
        Ez[n + 1, 1:-1] = Ez[n, 1:-1] + (coeff/e) * (Hy[n, 1:] - Hy[n, :-1])
        
        # conditions de bords
        Ez[n + 1, 0] = 0
        Ez[n + 1, -1] = 0
        
    Hy[N, :-1] = Hy[N-1, :-1] + (coeff/mu) * (Ez[N, 1:-1] - Ez[N, :-2])
    Hy[N, -1] = Hy[N-1, -1] + (coeff/mu) * (Ez[N, 1] - Ez[N, -1]) 
    return Hy, Ez, t, z, y

N = 1500
k = 6
omega = k * np.pi / 3  # pour assurer que E et H sont 1-periodique
T = 1.5

delta_t = T /N

# entre 0 et L/3
I1 = 200
l1= 0
L1 = 1

mu1 = 2
e1 = 1/3

Hy1, Ez1, t1, z1, y1 = Maxwell1d(N, I1, T, l1, L1, omega, mu1, e1)
H1, E1 = sol_exacte(t1, z1, y1, omega, mu1, e1)

delta_x1 = L1 / I1
CFL1 = delta_t / (mu1* e1 * delta_x1)
print(f"CFL1 = {CFL1}")

# entre L/3 et 2L/3
I2 = 300
l2 =1
L2 = 2

mu2 = 1
e2 = 1/2

Hy2, Ez2, t2, z2, y2 = Maxwell1d(N, I2, T, l2, L2, omega, mu2, e2)
H2, E2 = sol_exacte(t2, z2, y2, omega, mu2, e2)

delta_x2 = L2 / I2
CFL2 = delta_t / (mu2 * e2 * delta_x2)
print(f"CFL2 = {CFL2}")

# entre 2L/3 et L
I3 = 400
l3 =2
L3 = 3

mu3 = 2
e3 = 1/3

Hy3, Ez3, t3, z3, y3 = Maxwell1d(N, I3, T, l3, L3, omega, mu3, e3)
H3, E3 = sol_exacte(t3, z3, y3, omega, mu3, e3)

delta_x3 = L3 / I3
CFL3 = delta_t / (mu3 * e3 * delta_x3)
print(f"CFL3 = {CFL3}")

#graphe du champ electrique et du champ magnetique
plt.figure(num=1, figsize=(12,8))
t1 = 1500
plt.subplot(211)
plt.title(f"champ electrique pour t={t1}")

# entre 0 et L/3
plt.plot(z1, Ez1[t1,:], label='Fluide Electrique Yee (Ez)', color='g', lw = 2, linestyle = '--')#, marker='+')
plt.plot(z1, E1[t1,:], label='Fluide Electrique exact (E)', color='r', lw = 1)
# entre L/3 et 2L/3
plt.plot(z2, Ez2[t1,:], color='g', lw = 2, linestyle = '--')
plt.plot(z2, E2[t1,:], color='r', lw = 1)
# entre 2L/3 et L
plt.plot(z3, Ez3[t1,:], color='g', lw = 2, linestyle = '--')
plt.plot(z3, E3[t1,:], color='r', lw = 1)

plt.xlabel('z')
plt.ylabel('Ez')
plt.legend()


plt.subplot(212)
plt.title(f"champ magnetique t={t1}")

# entre 0 et L/3
plt.plot(y1, Hy1[t1,:], label='Fluide Magnetique Yee (Hy)', color='g', linestyle = '--', lw =2)#, marker='+')
plt.plot(y1, H1[t1,:], label='Fluide Magnetique exact (H)', color='r', lw=1)
# entre L/3 et 2L/3
plt.plot(y2, Hy2[t1,:], color='g', linestyle = '--', lw =2)
plt.plot(y2, H2[t1,:], color='r', lw=1)
# entre 2L/3 et L
plt.plot(y3, Hy3[t1,:], color='g', linestyle = '--', lw =2)
plt.plot(y3, H3[t1,:], color='r', lw=1)

plt.xlabel('y')
plt.ylabel('Hy')
plt.legend()

plt.tight_layout()
plt.show()

#%%  Animation de la variation du champ electrique et du champ magnetique

from matplotlib.animation import FuncAnimation

fig, ((ax1, ax3), (ax2, ax4)) = plt.subplots(2, 2, figsize=(12, 8))

def animate(n):
    ax1.clear()
    ax2.clear()
    ax3.clear()
    ax4.clear()

    ax1.set_title('Fluide Electrique par Yee (Ez) at t = {}*dt'.format(n))
    ax1.set_xlabel('En axe Z')
    ax1.set_ylabel('Ez')
    ax1.plot(z1, Ez1[n, :], color='g', label='Ez par Yee', lw=2.5, ls='--')
    ax1.plot(z1, E1[n, :], color='r', label='E exacte', lw =1) 
    ax1.plot(z2, Ez2[n, :], color='g', lw=2.5, ls='--')
    ax1.plot(z2, E2[n, :], color='r', lw =1)
    ax1.plot(z3, Ez3[n, :], color='g', lw=2.5, ls='--')
    ax1.plot(z3, E3[n, :], color='r', lw =1)
    ax1.set_ylim(-1.1, 1.1)
    ax1.legend()

    ax2.set_title('Fluide Magnetique par Yee (Hy) at t = {}*dt'.format(n))
    ax2.set_xlabel('En axe Y')
    ax2.set_ylabel('Hy')
    ax2.plot(y1, Hy1[n, :], color='b', label='Hy par Yee', lw=2.5, ls='--')
    ax2.plot(y1, H1[n, :], color='r', label='H exacte', lw=1)
    ax2.plot(y2, Hy2[n, :], color='b', lw=2.5, ls='--')
    ax2.plot(y2, H2[n, :], color='r', lw=1)
    ax2.plot(y3, Hy3[n, :], color='b', lw=2.5, ls='--')
    ax2.plot(y3, H3[n, :], color='r', lw=1)
    ax2.set_ylim(-1.1, 1.1)
    ax2.legend()
    
    ax3.set_title('Solution Electrique (E) at t = {}*dt'.format(n))
    ax3.set_xlabel('En axe Z')
    ax3.set_ylabel('E')
    ax3.plot(z1, E1[n, :], color='r', label='E')
    ax3.plot(z2, E2[n, :], color='r')
    ax3.plot(z3, E3[n, :], color='r')
    ax3.set_ylim(-1.1, 1.1)
    ax3.legend()

    
    ax4.set_title('Solution Magnetique (H) at t = {}*dt'.format(n))
    ax4.set_xlabel('En axe Y')
    ax4.set_ylabel('H')
    ax4.plot(y1, H1[n, :], color='c', label='H')
    ax4.plot(y2, H2[n, :], color='c')
    ax4.plot(y3, H3[n, :], color='c')
    ax4.set_ylim(-1.1, 1.1)
    ax4.legend()

    plt.tight_layout()

ani = FuncAnimation(fig, animate, frames = range(0, N, 8), interval=50)

plt.show()