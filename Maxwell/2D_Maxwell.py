# -*- coding: utf-8 -*-
"""
Created on Wed Jun 19 16:26:39 2024

@author: YY
"""

import numpy as np
from numpy import cos, sin, pi
import matplotlib.pyplot as plt



def Maxwell_2d(N, I, T, L, omega):
    
    delta_t = T/N
    delta_x = L/I
    delta_y = L/I
    coeff_x = delta_t / delta_x
    coeff_y = delta_t / delta_y
    
    x = np.arange(0, L + 0.1 * delta_x, delta_x)
    y = np.arange(0, L + 0.1 * delta_x, delta_x)
    t = np.arange(0, T + 0.1 * delta_t, delta_t)
    
    # forme des champs magnetique et eletrique
    Hx = np.zeros((N+1, I+1, I+1))
    Hy = np.zeros((N+1, I+1, I+1))
    Ez = np.zeros((N+1, I+1, I+1))
    
    # conditons initiales
    for i in range(I+1):
        Ez[0, i, :] = sin(omega*x[i]) * sin(omega*y[:]) * sin(np.sqrt(2)*omega* t[0])  
        Hx[0, i, :] = 1/np.sqrt(2) * sin(omega* x[i]) * cos(omega* (y[:] + delta_y/2)) * cos(np.sqrt(2)*omega*(t[0]+ delta_t/2)) 
        Hy[0, i, :] = -1/np.sqrt(2) * cos(omega* (x[i] + delta_x/2)) * sin(omega* y[:]) * cos(np.sqrt(2)*omega*(t[0]+ delta_t/2)) 
    
    
    for n in range(N):
        
        Ez[n+1, 1:-1, 1:-1] = Ez[n, 1:-1, 1:-1] + coeff_x * (Hy[n, 1:-1 , 1:-1] - Hy[n, 0:-2, 1:-1]) - coeff_y * (Hx[n, 1:-1, 1:-1] - Hx[n, 1:-1, 0:-2])
        Hx[n+1, 0:-1, 0:-1] = Hx[n, 0:-1, 0:-1] - coeff_y * (Ez[n+1, 0:-1, 1: ] - Ez[n+1, 0:-1, 0:-1])
        Hy[n+1, 0:-1, 0:-1] = Hy[n, 0:-1, 0:-1] + coeff_x * (Ez[n+1, 1: , 0:-1] - Ez[n+1, 0:-1, 0:-1]) 
        
        
        Ez[n+1, 0, :] = 0
        Ez[n+1, -1, :] = 0
        Ez[n+1, :, 0] = 0
        Ez[n+1, :, -1] = 0
        
        Hx[n+1, -1, :] = 0 
        Hx[n+1, :, -1] = 1/np.sqrt(2) * sin(omega* (x[:] +delta_x/2)) * cos(omega*(y[-1] + delta_y/2)) * cos(np.sqrt(2)*omega* (t[n+1]+ delta_t/2)) 
        
        Hy[n+1, :, -1] = 0
        Hy[n+1, -1, :] = -1/np.sqrt(2) * cos(omega* (x[-1] +delta_x/2)) * sin(omega*(y[:] + delta_y/2)) * cos(np.sqrt(2)*omega* (t[n+1]+ delta_t/2)) 
        
    return t, x, y, Ez, Hx, Hy


def solution(t, x, y, omega):
    delta_t = t[1]-t[0]
    Ez_s = np.zeros((len(t), len(x), len(y)))
    Hx_s = np.zeros((len(t), len(x), len(y)))
    Hy_s = np.zeros((len(t), len(x), len(y)))
    
    for n in range(len(t)):
        for i in range(len(x)):
            Ez_s[n,i,:] = sin(omega*x[i]) * sin(omega*y[:]) * sin(np.sqrt(2)*omega* t[n]) 
            Hx_s[n,i,:] = 1/np.sqrt(2) * sin(omega* x[i]) * cos(omega* (y[:] + delta_y/2)) * cos(np.sqrt(2)*omega* (t[n]+ delta_t/2)) 
            Hy_s[n,i,:] = -1/np.sqrt(2) * cos(omega* (x[i] +delta_x/2)) * sin(omega* y[:]) * cos(np.sqrt(2)*omega* (t[n]+ delta_t/2))
    
    return Ez_s, Hx_s, Hy_s

k = 4  # k pair pour que la solution soit periodique
L = 4
T = 3

N = 700
I = 250
delta_t = T/N
delta_x = L/I
delta_y = L/I
omega = k * pi / L

CFL = delta_t/delta_x + delta_t/delta_y

t, x, y, Ez, Hx, Hy = Maxwell_2d(N, I, T, L, omega)
Ez_s, Hx_s, Hy_s = solution(t, x, y, omega)


print(f"condition CFL = {CFL}")

#%%  

temps_indice = 700

# variation du champ electrique et du champ magnetique sur l'axe y
x_indice = 100

plt.figure(figsize=(18, 12))

plt.subplot(2, 3, 1)
plt.plot(y, Ez[temps_indice, x_indice, :], label='Yee', color = 'g', lw =2, ls= '--')#, marker='+')
plt.plot(y, Ez_s[temps_indice, x_indice, :], label='sol', color = 'r', lw =0.8)
plt.title('Ez')
plt.xlabel('Y axis')
plt.ylabel('fluide electrique')
plt.grid(True)
plt.legend()
#plt.ylim([-0.3, 0.3])

plt.subplot(2, 3, 2)
plt.plot(y, Hx[temps_indice, x_indice, :], label='Yee', color = 'g', lw =2, ls= '--')#, marker='+')
plt.plot(y, Hx_s[temps_indice, x_indice, :], label='sol', color = 'r', lw =0.8)
plt.title('Hx')
plt.xlabel('Y axis')
plt.ylabel('fuide magnetique')
plt.grid(True)
plt.legend()
#plt.ylim([-0.45, 0.45])

plt.subplot(2, 3, 3)
plt.plot(y, Hy[temps_indice, x_indice, :], label='Yee', color = 'g', lw =2, ls= '--')#, marker='+')
plt.plot(y, Hy_s[temps_indice, x_indice, :], label='sol', color = 'r', lw =0.8)
plt.title('Hy')
plt.xlabel('Y axis')
plt.ylabel('fuide magnetique')
plt.grid(True)
plt.legend()

plt.subplot(2, 3, 4)
plt.plot(y, Ez_s[temps_indice, x_indice, :], color = 'r', label='sol')
plt.title('Ez_s')
plt.xlabel('Y axis')
plt.ylabel('fuide electrique')
plt.grid(True)
plt.legend()

plt.subplot(2, 3, 5)
plt.plot(y, Hx_s[temps_indice, x_indice, :], color = 'r', label='sol')
plt.title('Hx_s')
plt.xlabel('Y axis')
plt.ylabel('fuide magnetique')
plt.grid(True)
plt.legend()

plt.subplot(2, 3, 6)
plt.plot(y, Hy_s[temps_indice, x_indice, :], color = 'r', label='sol')
plt.title('Hy_s')
plt.xlabel('Y axis')
plt.ylabel('fuide magnetique')
plt.grid(True)
plt.legend()

plt.tight_layout()
plt.show()


# variation du champ electrique et du champ magnetique sur l'axe x
y_indice = 100

plt.figure(figsize=(18, 12))

plt.subplot(2, 3, 1)
plt.plot(x, Ez[temps_indice, :, y_indice], label='Yee', color = 'g', lw =2, ls= '--')#, marker='+')
plt.plot(x, Ez_s[temps_indice, :, y_indice], label='sol', color = 'r', lw =0.8)
plt.title('Ez')
plt.xlabel('X axis')
plt.ylabel('fuide electrique')
plt.grid(True)
plt.legend()
#plt.ylim([-0.3, 0.3])

plt.subplot(2, 3, 2)
plt.plot(x, Hx[temps_indice, :, y_indice], label='Yee', color = 'g', lw =2, ls= '--')#, marker='+')
plt.plot(x, Hx_s[temps_indice, :, y_indice], label='sol', color = 'r', lw =0.8)
plt.title('Hx')
plt.xlabel('X axis')
plt.ylabel('fuide magnetique')
plt.grid(True)
plt.legend()

plt.subplot(2, 3, 3)
plt.plot(x, Hy[temps_indice, :, y_indice], label='Yee', color = 'g', lw =2, ls= '--')#, marker='+')
plt.plot(x, Hy_s[temps_indice, :, y_indice], label='sol', color = 'r', lw =0.8)
plt.title('Hy')
plt.xlabel('X axis')
plt.ylabel('fuide magnetique')
plt.grid(True)
plt.legend()

plt.subplot(2, 3, 4)
plt.plot(x, Ez_s[temps_indice, :, y_indice], color = 'r', label='sol')
plt.title('Ez_s')
plt.xlabel('X axis')
plt.ylabel('fuide electrique')
plt.grid(True)
plt.legend()

plt.subplot(2, 3, 5)
plt.plot(x, Hx_s[temps_indice, :, y_indice], color = 'r', label='sol')
plt.title('Hx_s')
plt.xlabel('X axis')
plt.ylabel('fuide magnetique')
plt.grid(True)
plt.legend()

plt.subplot(2, 3, 6)
plt.plot(x, Hy_s[temps_indice, :, y_indice], color = 'r', label='sol')
plt.title('Hy_s')
plt.xlabel('X axis')
plt.ylabel('fuide magnetique')
plt.grid(True)
plt.legend()

plt.tight_layout()
plt.show()

#%%  graphes du champ electrique et du champ magnetique

temps_indice = 700

plt.figure(figsize=(18, 12))

plt.subplot(2, 3, 1)
plt.contourf(x, y, Ez[temps_indice, :, :], cmap='magma')
plt.title('Ez')
plt.colorbar(label='Ez')
plt.xlabel('X axis')
plt.ylabel('Y axis')

plt.subplot(2, 3, 2)
plt.contourf(x, y, Hx[temps_indice, :, :], cmap='magma')
plt.title('Hx')
plt.colorbar(label='Hx')
plt.xlabel('X axis')
plt.ylabel('Y axis')

plt.subplot(2, 3, 3)
plt.contourf(x, y, Hy[temps_indice, :, :], cmap='magma')
plt.title('Hy')
plt.colorbar(label='Hy')
plt.xlabel('X axis')
plt.ylabel('Y axis')

plt.subplot(2, 3, 4)
plt.contourf(x, y, Ez_s[temps_indice, :, :], cmap='magma')
plt.title('Ez_s')
plt.colorbar(label='Ez_s')
plt.xlabel('X axis')
plt.ylabel('Y axis')

plt.subplot(2, 3, 5)
plt.contourf(x, y, Hx_s[temps_indice, :, :], cmap='magma')
plt.title('Hx_s')
plt.colorbar(label='Hx_s')
plt.xlabel('X axis')
plt.ylabel('Y axis')

plt.subplot(2, 3, 6)
plt.contourf(x, y, Hy_s[temps_indice, :, :], cmap='magma')
plt.title('Hy_s')
plt.colorbar(label='Hy_s')
plt.xlabel('X axis')
plt.ylabel('Y axis')

plt.tight_layout()
plt.show()

#%%   Animation de la variation du champ electrique et du champ magnetique
from matplotlib.animation import FuncAnimation

fig, axs = plt.subplots(2, 3, figsize=(18, 12))

def update(frame):
    axs[0, 0].clear()
    axs[0, 0].contourf(x, y, Ez[frame, :, :], cmap='magma')
    axs[0, 0].set_title('Ez at time = {:.2f}'.format(t[frame]))
    axs[0, 0].set_xlabel('X axis')
    axs[0, 0].set_ylabel('Y axis')
    
    axs[0, 1].clear()
    axs[0, 1].contourf(x, y, Hx[frame, :, :], cmap='magma')
    axs[0, 1].set_title('Hx at time = {:.2f}'.format(t[frame]))
    axs[0, 1].set_xlabel('X axis')
    axs[0, 1].set_ylabel('Y axis')
    
    axs[0, 2].clear()
    axs[0, 2].contourf(x, y, Hy[frame, :, :], cmap='magma')
    axs[0, 2].set_title('Hy at time = {:.2f}'.format(t[frame]))
    axs[0, 2].set_xlabel('X axis')
    axs[0, 2].set_ylabel('Y axis')
    
    axs[1, 0].clear()
    axs[1, 0].contourf(x, y, Ez_s[frame, :, :], cmap='magma')
    axs[1, 0].set_title('Ez_s at time = {:.2f}'.format(t[frame]))
    axs[1, 0].set_xlabel('X axis')
    axs[1, 0].set_ylabel('Y axis')
    
    axs[1, 1].clear()
    axs[1, 1].contourf(x, y, Hx_s[frame, :, :], cmap='magma')
    axs[1, 1].set_title('Hx_s at time = {:.2f}'.format(t[frame]))
    axs[1, 1].set_xlabel('X axis')
    axs[1, 1].set_ylabel('Y axis')
    
    
    axs[1, 2].clear()
    axs[1, 2].contourf(x, y, Hy_s[frame, :, :], cmap='magma')
    axs[1, 2].set_title('Hy_s at time = {:.2f}'.format(t[frame]))
    axs[1, 2].set_xlabel('X axis')
    axs[1, 2].set_ylabel('Y axis')
print(len(t)) 

ani = FuncAnimation(fig, update, frames=50, interval=100)
plt.show()


#%%  Calcul de l'ordre de convergence

# un ensembre de nombre de discritisation spatiale
I_list  = np.arange(250, 280, 10) 
# un ensembre de nombre de discritisation temporelle, pour la condition de CFL: delta_t/delta_x <= 1
N_list = 3 * (T * I_list // L) 


err_Hx = np.zeros(len(I_list))
err_Hy = np.zeros(len(I_list))
err_Ez = np.zeros(len(I_list))

for i, I in enumerate(I_list):
    
    err_xxx = np.zeros((N, I)) # erreur de Hx sur l'axe y 
    err_yyy = np.zeros((N, I)) # erreur de Hy sur l'axe y 
    err_zzz = np.zeros((N, I)) # erreur de Ez sur l'axe y 
    
    err_xx = np.zeros(N) # erreur de Hx sur l'axe x
    err_yy = np.zeros(N) # erreur de Hy sur l'axe x
    err_zz = np.zeros(N) # erreur de Ez sur l'axe x
    
    delta_x = L / I
    delta_t = T / N_list[i]
    print("CFL =",delta_t/delta_x + delta_t/delta_y)
    
    t, x, y, Ez, Hx, Hy = Maxwell_2d(N_list[i], I, T, L, omega)
    Ez_s, Hx_s, Hy_s = solution(t, x, y, omega)
    
    for n in range(N_list[i]):
        for j in range(I): 
            # L'effet de l'axe y sur l'erreur 
            err_xxx[n, j] = np.linalg.norm(Hx_s[n, j, :] - Hx[n, j, :], ord=2) *delta_x**(1/2)
            err_yyy[n, j] = np.linalg.norm(Hy_s[n, j, :] - Hy[n, j, :], ord=2) *delta_x**(1/2)
            err_zzz[n, j] = np.linalg.norm(Ez_s[n, j, :] - Ez[n, j, :], ord=2) *delta_x**(1/2)
        
        # L'effet de l'axe y du temps et de l'espace sur l'erreur
        err_xx[n] = np.linalg.norm(err_xxx, ord=2) *delta_t**(1/2)
        err_yy[n] = np.linalg.norm(err_yyy, ord=2) *delta_t**(1/2)
        err_zz[n] = np.linalg.norm(err_zzz, ord=2) *delta_t**(1/2)
    
    #L'effet de l'axe y et l'axe x du temps et de l'espace sur l'erreur
    err_Hx[i] = np.linalg.norm(err_xx, ord=2) *delta_t**(1/2)
    err_Hy[i] = np.linalg.norm(err_yy, ord=2) *delta_t**(1/2)
    err_Ez[i] = np.linalg.norm(err_zz, ord=2) *delta_t**(1/2)

h = L / I_list 

Hx_ordre = np.polyfit(np.log(h), np.log(err_Hx), 1)[0]
Hy_ordre = np.polyfit(np.log(h), np.log(err_Hy), 1)[0]
Ez_ordre = np.polyfit(np.log(h), np.log(err_Ez), 1)[0]

print(Hx_ordre)
print(Hy_ordre)
print(Ez_ordre)

#%%
# figure d'erreur
plt.figure()
plt.plot(h, err_Ez, label="erreur de Ez", color='r')
plt.plot(h, err_Hx, label="erreur de Hx", color='b', ls='--')
plt.plot(h, err_Hy, label="erreur de Hy", color='g', ls='-.')
plt.title("L'erreur du schema")
plt.xlabel('delta_x')
plt.ylabel('erreur')
plt.legend()
plt.show()


# figure d'ordre de convergence
plt.figure()
plt.plot(np.log(h), np.log(err_Ez), label=f"ordre d'erreur de Ez = {round(Ez_ordre,2)}", color='r')
plt.plot(np.log(h), np.log(err_Hx), label=f"ordre d'erreur de Hx = {round(Hx_ordre,2)}", color='b', ls='--')
plt.plot(np.log(h), np.log(err_Hy), label=f"ordre d'erreur de Hy = {round(Hy_ordre,2)}", color='g', ls='-.')
plt.title("L'ordre de convergence")
plt.xlabel('log(delta_x)')
plt.ylabel('log(erreur)')
plt.legend()
plt.show()
