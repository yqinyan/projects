import matplotlib.pyplot as plt
import numpy as np
from math import comb

## Plot ##

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

## Dimension 2 ##

def reflection2d(W,M) :
    return 2*np.array(M) - np.array(W)

def extension2d(R,M,p=2) :
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
        R = reflection2d(W,M)
       
        if f(R) < f(G) : # Extend or reflect
            if f(B) < f(R) :
                W = R
            else :
                E = extension2d(R,M)
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

## Dimension n ##

def isobarycenter(V) : 
    n = len(V[0])
    bary = [0] * n
    for v in V :
        for i in range(n) :
            bary[i] += v[i]
    return [x / len(V) for x in bary]

def reflection(x0,x,a=1) :
    return np.array(x0) + a*(np.array(x0) - np.array(x))

def homothetie(x0,xr,k=2) :
    return np.array(x0) + k*(np.array(xr)-np.array(x0))

def NelderMead(f,V_init,eps=0.01,e=2,k=2,s=2) : # e,k extension and contraction ratios
    V = [V_init]
    while True :
        L = sorted(V[-1], key=f)
        x0 = isobarycenter(L[:-1])
        xr = reflection(x0,L[-1])
        
        if f(L[0]) <= f(xr) and f(xr) < f(L[-2]) :
            L[-1] = xr
        elif f(xr) < f(L[0]) :
            xe = homothetie(x0,xr,e)
            if f(xe) <= f(xr) : 
                L[-1] = xe
            else :
                L[-1] = xr
        elif f(L[-2]) <= f(xr) :
            xc = homothetie(x0,L[-1],1/k)
            if f(xc) < f(L[-1]) :
                L[-1] = xc
            else :
                L = [homothetie(L[0],x,1/2) for x in L]
        V.append(L)
        if max(np.linalg.norm(np.array(L[i]) - np.array(L[j])) for i in range(len(L)) for j in range(i+1, len(L))) <= eps :
            break
    return V

#%%
###### Problem

# domain length
l = 1
# number of elements
n_el = 200
# element size
h = l/n_el
# number of nodes
n_no = n_el+1
# nodes coordinates
X = np.linspace(0,1,n_no)
# integration points coordinates
Xg = 0.5*( X[:-1] + X[1:] )
# target temperature
T_star = np.ones(X.shape) + (X/l)**2
# initial design vector that parameterizes heat sources
dim_opt = 6
Var_ini = np.full(dim_opt,0)
# conduction coefficient with possible noise
rng = np.random.default_rng(seed=2024)

def compute_conduction(noise_val):
    noise = rng.random(Xg.shape) * noise_val
    K = np.ones(Xg.shape) +  0.4*np.sin(4*np.pi*Xg/l) +  0.3*np.sin(12*np.pi*Xg/l) + noise
    return K

# without noise
K_ref = compute_conduction(0)

def bernstein_poly(i, n, x):
    return comb(n, i) * (x ** i) * ((1 - x) ** (n - i))

def compute_source(Var_opt):  

    S = np.ones(Xg.shape)
    for i in range(len(Var_opt)):
        #S += Var_opt[i] * np.sin((i+1) * np.pi * Xg)
        #S += Var_opt[i] * bernstein_poly(i, len(Var_opt) - 1, Xg)
        S += Var_opt[i] * Xg**i
    return S

def compute_matrix(K):   
    # Finite-Element matrix initialisation
    M = np.zeros((n_no,n_no))
    
    # Boundary conditions
    M[0][0] = 1
    M[n_no-1][n_no-1] = 1
    
    # Internal coeff
    for i in range(1,n_no-1):
        M[i][i] = (K[i-1]+K[i])/h
        M[i][i-1] = -K[i-1]/h
        M[i][i+1] = -K[i]/h    
    return M

def compute_rhs(S):
    # Finite-Element right-hand side initialisation
    Rhs = np.zeros((n_no,1))
    
    # Boundary conditions
    Rhs[0] = 1
    Rhs[n_no-1] = 2
    
    # internal coeff
    for i in range(1,n_no-1):
        Rhs[i] = (S[i-1]+S[i])*h/2
    return Rhs

def simulator(noise,Var):
    # conduction
    K = compute_conduction(noise)

    # matrix
    M = compute_matrix(K)
    
    # compute heat source
    Src = compute_source(Var)

    # right-hand side
    Rhs = compute_rhs(Src)

    # Finite-element solution
    T = np.matmul(np.linalg.inv(M),Rhs).reshape((n_no))
    return T

def J(S) :
    T = simulator(0,S)
    r = 0
    for i in range(len(T)):
        r += 0.5*((T[i]-T_star[i])**2)*h
    return r

if __name__ == '__main__' :
    #V = NelderMead(f,[[0,0],[1,0],[0,1]])
    #draw_sequel(V)
    #plt.plot()

    B = np.column_stack((np.array([0 for i in range(dim_opt)]),-3*np.eye(dim_opt)))
    G = NelderMead(J,B)
    
    S_final = G[-1][0]  # Extract the optimal source
    T_final = simulator(0, S_final)
    T_ini = simulator(0, Var_ini)
    
    #S_final = compute_source(G[-1])
    #T_final = simulator(0,S_final)
    #T_ini = simulator(0,Var_ini)
    
    plt.plot(X,T_ini,'b', label = 'initial')
    plt.plot(X,T_final,'g', label ='temperature_modifiee')
    plt.plot(X,T_star,'r--', label = 'cible')
    
    plt.title(f'dim_opt={dim_opt}')
    plt.xlabel('x')
    plt.ylabel('temperature')
    plt.legend()
    #plt.legend(["initial","", "cible"], loc="lower right")
    #plt.savefig('figure de {} dim avec Bernstein'.format(dim_opt))
    plt.show()
    
#%%
print(S_final)
print(J(S_final))

plt.figure()
plt.plot(X[:,-2],compute_source(S_final))
#%%  erreur sous differents dimension

def erreur(S, noise):
    return np.linalg.norm(simulator(noise, S) - T_star)

dim_opt_set = np.arange(2, 8, 1) #differentes dimensions
eps = 0.02
S = [None] * len(dim_opt_set)
err = [None] * len(dim_opt_set)

for i in range(len(dim_opt_set)):
    
    B = np.column_stack((np.array([0 for i in range(dim_opt_set[i])]),-3*np.eye(dim_opt_set[i])))
    G = NelderMead(J,B)
    S[i] = G[-1][0]
    err[i] = erreur(S[i], 0)
    if err[i] < eps:
        break
print(f'Arrêt à dimension = {dim_opt_set[i]}')
    
plt.figure()
plt.plot(dim_opt_set,err)
plt.title('erreur sous differentes dimensions')
plt.xlabel('dimension')
plt.ylabel('erreur')
plt.grid(True)
plt.legend()
plt.show()
    
#%% erreure sur positions differentes 
dim_opt_set = [6]
def erreur_x(S, noise):
    err_x = []
    for i in range(len(X)):
        #err_x.append(simulator(noise, S)[i] - T_star[i])
        err_x.append(np.linalg.norm(simulator(noise, S)[i] - T_star[i]))
    return err_x

for i in range(len(dim_opt_set)):
    
    B = np.column_stack((np.array([0 for i in range(dim_opt_set[i])]),-3*np.eye(dim_opt_set[i])))
    G = NelderMead(J,B)
    S[i] = G[-1][0]
    
    err_x = erreur_x(S[i], 0)
    plt.figure()
    plt.plot(X,err_x)
    plt.title(f'erreur le long de X pour dimension={dim_opt_set[i]}')
    plt.xlabel('x')
    plt.ylabel('erreur')
    plt.grid(True)
    plt.legend()
    plt.show()

