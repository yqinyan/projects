# -*- coding: utf-8 -*-
"""
Created on Fri Oct  4 11:44:54 2024

@author: YY
"""


import matplotlib.pyplot as plt
import numpy as np

def draw_figure_2d(V) :
    X = [v[0] for v in V] + [V[0][0]]
    Y = [v[1] for v in V] + [V[0][1]]
    plt.plot(X,Y,'-',linewidth=0.5)
    pass

def draw_sequel(V) :
    for v in V :
        draw_figure_2d(v)
    pass

def reflection(W,M) :
    return 2*np.array(M) - np.array(W)

def extension(R,M,p=2) :
    return p*np.array(R) - np.array(M)

def f(x) :
    return x[0]**2 + x[1]**2 - 4*x[0] - x[1] - x[0]*x[1]

def NelderMead_2d(f,V_init,eps=0.1) :
   
    V = [sorted(V_init,key=f)]
   
    while True :
        L = sorted(V[-1],key=f)
        B = L[0] #best
        G = L[1] #good
        W = L[2] #worst

        M = [(B[0]+G[0]) / 2, (B[1]+G[1]) / 2]
        R = reflection(W,M)
       
        if f(R) < f(G) : # Extend or reflect
            if f(B) < f(R) :
                W = R
            else :
                E = extension(R,M)
                if f(E) < f(B) :
                    W = E
                else :
                    W = R
               
        else : # contract or shrink
            if f(R) < f(W) :
                W = R
            C1 = [(W[0] + M[0]) / 2 , (W[1] + M[1]) / 2]
            C2 = [(R[0] + M[0]) / 2 , (R[1] + M[1]) / 2]
            C = min(C1,C2,key=f)
            if f(C) < f(W) :
                W = C
            else :
                S = [(W[0] + B[0]) / 2 , (W[1] + B[1]) / 2]
                W = S
                G = M
        V.append([B,G,W])

        if np.linalg.norm(np.array(B) - np.array(G)) < eps : #Precision
            break
    return V

V = NelderMead_2d(f,[[0,0],[1.2,0],[0,0.8]])
draw_sequel(V)

plt.show()
