# -*- coding: utf-8 -*-
"""
Created on Fri Oct  4 10:55:46 2024

@author: YY
"""

import numpy as np
import matplotlib.pyplot as plt



def f(x):
    return x[0]**2 + x[1]**2 - 4*x[0] - x[1] - x[0]*x[1]


def Nelder_Mead(x, f, maxit, rho=1, gamma=2, alpha=1/2, sigma=1/2, tolerance=0.1):
    
    n1,n2 = x.shape
    simplex = [x[:,i].T for i in range(n2)] 
    simplex = np.array(simplex)  
    nit = 0
    all_simplex = [simplex.copy()]
    
    while nit < maxit:
        simplex_sorted = sorted(simplex, key = f, reverse=False)
        all_simplex.append(simplex_sorted)
        
        m = np.mean(simplex_sorted[0:-1], axis=0)
        
        r = m + rho* (m- simplex_sorted[-1])  #reflection point 
        
        if f(simplex_sorted[0]) <= f(r) < f(simplex_sorted[-2]):
            simplex_sorted[-1] = r
            
        
        elif f(r) < f(simplex_sorted[0]):
            e = m + gamma* (r - m) #extension point
            if f(e) < f(r):
                simplex_sorted[-1] = e
            else:
                simplex_sorted[-1] = r
            
        
        elif f(r) < f(simplex_sorted[-1]):
            co = m + alpha* (r - m) #contraction point outside
            if f(co) < f(simplex_sorted[-1]):
                simplex_sorted[-1] = co
            
        
        else:
            ci = m + alpha* (simplex_sorted[-1] - m) #contractoion point inside
            if f(ci) < f(simplex_sorted[-1]):
                simplex_sorted[-1] = ci
            else:
                simplex_sorted[1:] = simplex_sorted[0] + sigma* (simplex_sorted[1:] - simplex_sorted[0])
                
        simplex = simplex_sorted
        
        if np.linalg.norm(simplex_sorted[-1] - simplex_sorted[0]) < tolerance:
            break
        
        nit +=1
    
    best_point = simplex_sorted[0]  
    best_value = f(best_point)       
    return all_simplex, best_point, best_value, nit  
        

def draw_figure_2d(V, color):
    X = [v[0] for v in V] + [V[0][0]]  
    Y = [v[1] for v in V] + [V[0][1]]
    plt.plot(X, Y, '-', linewidth=0.7, color=color)
    plt.scatter(X, Y, color=color)  
    pass


def draw_sequel(V):
    colors = ['r', 'g', 'b', 'y', 'c']  # 使用不同颜色绘制每步的三角形
    for idx, v in enumerate(V):
        color = colors[idx % len(colors)]  # 循环使用颜色
        draw_figure_2d(v, color)


maxit = 10
x = np.array([[0, 0], [1.2, 0], [0, 0.8]]).T


all_simplex, best_point, best_value, iterations = Nelder_Mead(x, f, maxit)
print('ss',Nelder_Mead(x, f, maxit)[0][-1])
print(all_simplex[-1])
draw_sequel(all_simplex)


plt.title("Nelder-Mead Simplex Evolution")
plt.xlabel("x")
plt.ylabel("y")
plt.grid(True)
plt.show()

print(f"Best point: {best_point}")
print(f"Best value: {best_value}")
print(f"Iterations: {iterations}")
