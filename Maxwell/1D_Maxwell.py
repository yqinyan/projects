# -*- coding: utf-8 -*-
"""
Created on Wed Jun 12 18:41:51 2024

@author: YY
"""

import numpy as np
import matplotlib.pyplot as plt

def sol_exacte(t, x, y, omega, mu, e):
    
    c = 1/np.sqrt(mu*e)    
    E = np.zeros((len(t), len(x)))
    H = np.zeros((len(t), len(x)))
    
    for n in range(len(t)):
        E[n, :] = np.sin(omega * x[:]) * np.sin(c * omega * t[n])
        
        H[n, :] = - 1/(mu*c) * np.cos(c * omega * (delta_t/2+t[n])) * np.cos(omega * y[:]) 
    return H, E


def sol_exacte(t, x, y, omega):
    
    delta_t = t[1] - t[0]
    
    E = np.zeros((len(t), len(x)))
    H = np.zeros((len(t), len(x)))
    
    for n in range(len(t)):
        E[n, :] = np.sin(omega * x[:]) * np.sin(omega * t[n])
        # pour H on a un decalage en temps de delta_t/2
        H[n, :] = - np.cos(omega * (delta_t/2+t[n])) * np.cos(omega * y[:]) 
    return H, E



def Maxwell1d(N, I, T, L, omega):
    
    delta_t = T/N  #le pas de discritisation temporelle
    delta_x = L/I  #le pas de discritisation spatiale
    coeff = delta_t/delta_x
    
    z = np.arange(0, L + 0.1 * delta_x, delta_x) # I+1 points de discritisation sur pour Ez sur [0,L]
    y = np.arange(delta_x/2, L+delta_x/2 + 0.1 * delta_x, delta_x) # I+1 points de discritisation pour Hy sur [delta_x/2,L+delta_x/2]
    t = np.arange(0, T + 0.1 * delta_t, delta_t) # N+1 points de discritisation temporelle sur [0,L]
    
    # forme des champs magnetique et eletrique
    
    Hy = np.zeros((N+1,I+1))
    Ez = np.zeros((N+1, I+1))
    
    # conditions initiales 
    
    Hy[0, :] = - np.cos(omega * (delta_t/2+t[0])) * np.cos(omega * y[:])
    Ez[0, :] = np.sin(omega * z[:]) * np.sin(omega * t[0])
    
    # le schema: Ez[1, i] = Ez[0, i] + coeff * (Hy[0, i] - Hy[0, i-1]), pour i dans {0,...I-1}
    Ez[1, 1:-1] = Ez[0, 1:-1] + coeff * (Hy[0, 1:-1] - Hy[0, :-2])
    
    # condtions de bords pour temps indice=1
    Ez[1, 0] = 0
    Ez[1, -1] = 0  
    
    # mise a jour des champs electrique et magnetique
    for n in range(1,N):
        
        # le schema: Hy[n, i] = Hy[n-1, i] + coeff * (Ez[n, i+1] - Ez[n, i]), pour i dans {0,...I-1}
        Hy[n, :-1] = Hy[n-1, :-1] + coeff * (Ez[n, 1:] - Ez[n, :-1])
        Hy[n, -1] = Hy[n-1, -1] + coeff * (Ez[n, 1] - Ez[n, -1]) 
        
        # le schema: Ez[n + 1, i] = Ez[n, i] + coeff * (Hy[n, i] - Hy[n, i-1]), pour i dans {0,...I-1}
        Ez[n + 1, 1:-1] = Ez[n, 1:-1] + coeff * (Hy[n, 1:-1] - Hy[n, :-2])
        
        # condtions de bords 
        Ez[n + 1, 0] = 0
        Ez[n + 1, -1] = 0
    
    # mise a jour de Hy pour le dernier indice temporelle
    Hy[N, :-1] = Hy[N-1, :-1] + coeff * (Ez[N, 1:] - Ez[N, :-1])
    Hy[N, -1] = Hy[N-1, -1] + coeff * (Ez[N, 1] - Ez[N, -1]) 
        
        
    return Hy, Ez, t, z, y


k = 2
L = 5
T = 4

omega = k * np.pi / L

N = 1500 # on a N+1=1501 points de discritisation temporelle
I = 600 ## on a I+1=601 points de discritisation spatiale

delta_t = T /N
delta_x = L / I
# dans ce cas pour verifier la condition CFL, il suffit de prendre que delta_t / delta_x < 1
CFL = delta_t / delta_x  
print(f"CFL = {delta_t/delta_x}")

Hy, Ez, t, x, y = Maxwell1d(N, I, T, L, omega)

H, E = sol_exacte(t, x, y, omega)




#graphe du champ electrique et du champ magnetique
plt.figure(num=1, figsize=(12,8))
t1 = 600
plt.subplot(211)
plt.title(f"champ electrique pour t={t1}")
plt.plot(x, Ez[t1,:], label='Fluide Electrique Yee (Ez)', color='g', linestyle = '--', lw = 2)
plt.plot(x, E[t1,:], label='Fluide Electrique Exact (E)', color='r', lw = 1)
plt.xlabel('z')
plt.ylabel('Ez')
plt.legend()


plt.subplot(212)
plt.title(f"champ magnetique t={t1}")
plt.plot(y, Hy[t1,:], label='Fluide Magnetique Yee (Hy)', color='g', linestyle = '--', lw =2)
plt.plot(y, H[t1,:], label='Fluide Magnetique Exact (H)', color='r', lw=1.5)
plt.xlabel('y')
plt.ylabel('Hy')
plt.legend()

plt.tight_layout()
plt.show()


plt.figure(num=2, figsize=(12,8))
t2 = 1500
plt.subplot(211)
plt.title(f"champ electrique pour t={t2}")
plt.plot(x, Ez[t2,:], label='Fluide Electrique Yee (Ez)', color='g', linestyle = '--', lw = 2)
plt.plot(x, E[t2,:], label='Fluide Electrique Exact (E)', color='r', lw = 1)
plt.xlabel('z')
plt.ylabel('Ez')
plt.legend()


plt.subplot(212)
plt.title(f"champ magnetique pour t={t2}")
plt.plot(y, Hy[t2,:], label='Fluide Magnetique Yee (Hy)', color='g', lw = 2, linestyle = '--')
plt.plot(y, H[t2,:], label='Fluide Magnetique Exact (H)', color='r', lw = 1)
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
    ax1.plot(Ez[n, :], color='g', label='Ez par Yee', lw=2.5, ls='--')
    ax1.plot(E[n, :], color='r', label='E exacte', lw =1)
    ax1.set_ylim([-1.1, 1.1])
    ax1.legend()

    ax2.set_title('Fluide Magnetique par Yee (Hy) at t = {}*dt'.format(n))
    ax2.set_xlabel('En axe Y')
    ax2.set_ylabel('Hy')
    ax2.plot(Hy[n, :], color='g', label='Hy par Yee', lw=2.5, ls='--')
    ax2.plot(H[n, :], color='r', label='H exacte', lw=1)
    ax2.set_ylim([-1.2, 1.2])
    #ax2.set_ylim([-0.1, 0.1])
    ax2.legend()
    
    ax3.set_title('Solution Electrique (E) at t = {}*dt'.format(n))
    ax3.set_xlabel('En axe Z')
    ax3.set_ylabel('E')
    ax3.plot(E[n, :], color='r', label='E')
    ax3.set_ylim([-1.1, 1.1])
    ax3.legend()

    
    ax4.set_title('Solution Magnetique (H) at t = {}*dt'.format(n))
    ax4.set_xlabel('En axe Y')
    ax4.set_ylabel('H')
    ax4.plot(H[n, :], color='c', label='H')
    ax4.set_ylim([-1.2, 1.2])
    #ax4.set_ylim([-0.1, 0.1])
    ax4.legend()

    plt.tight_layout()

ani = FuncAnimation(fig, animate, frames = range(0, N, 25), interval=50)

plt.show()


#%%   Calcul de l'ordre de convergence
k = 4
T = 3
L = 3
omega = k * np.pi / L

I_list  = np.arange(800, 2002, 200) # un ensembre de nombre de discritisation spatiale
N_list = 3 * (T * I_list // L) # # un ensembre de nombre de discritisation temporelle, pour la condition de CFL: delta_t/delta_x= 1/3 


err_H = np.zeros(len(I_list))
err_E = np.zeros(len(I_list))

for i, N in enumerate(N_list):
    
    err_HH = np.zeros(N)
    err_EE = np.zeros(N)
    
    # mise a jour du pas spatial et du pas temporel
    delta_x = L / I_list[i]
    delta_t = T / N
    
    # mise a jour de la solution numerique et la solution exacte
    Hy, Ez, t, z, y = Maxwell1d(N, I_list[i], T, L, omega)
    H, E = sol_exacte(t, z, y, omega)
    
    for j in range(N):
        
        # L'effet du temps sur l'erreur 
        err_HH[j] = np.linalg.norm(H[j, :] - Hy[j, :], ord=2) *delta_x**(1/2)
        #err_HH[j] = np.linalg.norm(H[j, :] - Hy[j, :], ord=np.inf)
        err_EE[j] = np.linalg.norm(E[j, :] - Ez[j, :], ord=2) *delta_x**(1/2)
        #err_EE[j] = np.linalg.norm(E[j, :] - Ez[j, :], ord=np.inf)
    
    #L'effet de l'axe x du temps et de l'espace sur l'erreur 
    err_H[i] = np.linalg.norm(err_HH, ord=2) *delta_t**(1/2)
    err_E[i] = np.linalg.norm(err_EE, ord=2) *delta_t**(1/2)
    #err_H[i] = np.linalg.norm(err_HH, ord=np.inf)
    #err_E[i] = np.linalg.norm(err_EE, ord=np.inf)

#la norme 2 marche mieux que la norme infinie

#print(err_E)
#print(err_H)

h = L / I_list  # une liste du pas de discritisation spatialle 
print(h)
H_ordre = np.polyfit(np.log(h), np.log(err_H), 1)[0] # calculer la pente de la courbe de log(h) et log(err_H)
E_ordre = np.polyfit(np.log(h), np.log(err_E), 1)[0] # calculer la pente de la courbe de log(h) et log(err_E)
print(f"Ordre de convergence spatial pour H: {H_ordre}")
print(f"Ordre de convergence spatial pour E: {E_ordre}")

print(f'max err_H = {max(err_H)}, max err_E = {max(err_E)}')

# figure d'erreur
plt.figure()
plt.plot(h, err_E, label="erreur de E", color='r')
plt.plot(h, err_H, label="erreur de H", color='b', ls='--')
plt.title("L'erreur du schema")
plt.xlabel('delta_x')
plt.ylabel('erreur')
plt.legend()
plt.show()


# figure d'ordre de convergence
plt.figure()
plt.plot(np.log(h), np.log(err_E), label=f"ordre d'erreur de E = {round(E_ordre,2)}", color='r')
plt.plot(np.log(h), np.log(err_H), label=f"ordre d'erreur de H = {round(H_ordre,2)}", color='b', ls='--')
plt.title("L'ordre de convergence")
plt.xlabel('log(delta_x)')
plt.ylabel('log(erreur)')
plt.legend()
plt.show()

cfl = (T/N_list[:]) / (L/I_list[:])
print("cfl = " + str(cfl))

