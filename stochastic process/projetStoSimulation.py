
#Libs

import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats


#%% Justification des paramètres par rapport au problème réel

Y0 = np.linspace(5, 25, 51) * 1e6 #capital entre 5 - 25 millions
alpha = 0.1 * 1e6 # le casino gagne 0.1 million per jour (100k/jour)
mu = alpha * (1 - 0.05)  #alpha (5% >) mu
mu = alpha * (1 + 0.05)  #alpha (5% <) mu
mu = alpha * (1 - 0.00)  #alpha = mu
annees = 6 #temps arbitraire que nous avons choisi pour visualiser la simulation

#%% modèle simulation individuel

### Paramètres du modèle

Y0 = 25 *1e6 # Capital initial du casino
alpha = 0.1 *1e6   # Taux de gain du casino par unité de temps
mu = alpha * (1 - 0.05)    # Moyenne des gains des joueurs (distribution exponentielle)
lambda_gain = 1 / mu  # Paramètre lambda pour la loi exponentielle des gains (variable gamma)
lambda_temps_gains = 1 # Paramètre lambda pour la loi exponentielle des temps de gains
annees = 6 #temps de simulation
#------------------------------------------------------------------------------------------

### Simuler un processus de Poisson pour les instants de gain

temps_total = annees * 365 # Durée totale de la simulation

#pour s'assurer que la somme cumulée est toujours supérieure au temps de simulation
size_gains = (lambda_temps_gains) * (temps_total + 100) 

#somme chi_i -> T_i
temps_gains = np.random.exponential(scale = 1/lambda_temps_gains, size = size_gains).cumsum() 
temps_gains = temps_gains[temps_gains <= temps_total]

#------------------------------------------------------------------------------------------

### Simuler les gains des joueurs
gains_joueurs = np.random.exponential(scale=1/lambda_gain, size=len(temps_gains)) #variable X_i

#------------------------------------------------------------------------------------------

### Calculer le capital du casino au fil du temps

temps = np.arange(temps_total + 1) # Durée totale de la simulation
capital_casino = np.zeros_like(temps) #vecteur composé de zéros pour stocker l'évolution du capital
index_gains = 0
somme_gains = 0

for i, t in enumerate(temps):
    
    while index_gains < len(temps_gains) and temps_gains[index_gains] <= t:
        # cond 1: Éviter un Débordement de l'Index
        # cond 2: Condition N_t definition
        somme_gains += gains_joueurs[index_gains]
        index_gains += 1
        
    #somme_gains accumulés jusqu'à présent t
    capital_casino[i] = Y0 + alpha * t - somme_gains
    
    if capital_casino[i] < 0:
        break

#------------------------------------------------------------------------------------------

### Visualiser la simulation

annees_temps = temps / 365
capital_casino_millions = capital_casino  #/ 1e6

plt.figure(figsize=(12, 6))
plt.plot(annees_temps, capital_casino_millions, label='Capital du Casino')
plt.axhline(Y0 / 1e6, color='grey', linestyle='--', label='Capital Initial')
plt.title('Simulation du Capital du Casino au Fil du Temps')
plt.xlabel('Années')
plt.ylabel('Capital casino (millions)')
plt.legend()
plt.grid(True)
plt.show()


#%% Fonctions
#------------------------------------------------------------------------------------------

### Probabilité

# Définition de la probabilité exacte si les gains des joueurs 
# sont distribués de manière exponentielle
def r(Y0,mu,alpha):
    gamma = 1/mu
    return np.exp(-Y0 * (gamma - (1/alpha)))/(alpha * gamma)


# essai de une simulation de l'évolution du capital du casino.
def test_simulation(Y0,mu,alpha,annees):
    
    lambda_gain = 1 / mu  # Paramètre lambda pour la loi exponentielle des gains (variable gamma)
    lambda_temps_gains = 1 # Paramètre lambda pour la loi exponentielle des temps de gains
    #------------------------------------------------------------------------------------------
    
    ### Simuler un processus de Poisson pour les instants de gain
    
    temps_total = annees * 365 # Durée totale de la simulation
    #pour s'assurer que la somme cumulée est toujours supérieure au temps de simulation
    size_gains = (lambda_temps_gains) * (temps_total + 100) 
    #somme chi_i -> T_i
    temps_gains = np.random.exponential(scale = 1/lambda_temps_gains, size = size_gains).cumsum() 
    temps_gains = temps_gains[temps_gains <= temps_total]
    
    #------------------------------------------------------------------------------------------
    
    ### Simuler les gains des joueurs
    gains_joueurs = np.random.exponential(scale=1/lambda_gain, size=len(temps_gains)) #variable X_i
    
    #------------------------------------------------------------------------------------------
    
    ### Calculer le capital du casino au fil du temps
    
    temps = np.arange(temps_total + 1)
    capital_casino = np.zeros_like(temps) #vecteur composé de zéros pour stocker l'évolution du capital
    index_gains = 0
    somme_gains = 0
    
    for i, t in enumerate(temps):
        
        while index_gains < len(temps_gains) and temps_gains[index_gains] <= t:
            # cond 1: Éviter un Débordement de l'Index
            # cond 2: Condition N_t definition
            somme_gains += gains_joueurs[index_gains]
            index_gains += 1
            
        #somme_gains accumulés jusqu'à présent t
        capital_casino[i] = Y0 + alpha * t - somme_gains
        
        if capital_casino[i] < 0:
            return 1
    
    return 0


# r(y) = P( \exists t \in R+ t.q. Y_t < 0 | Y_0 = y )
def r_simulation(Y0,mu,alpha,annees,nombre_tests):
    
    tests = np.zeros(nombre_tests) #vecteur composé de zéros pour stocker si Y_t < 0
    
    for i in range(nombre_tests):
        
        tests[i] = test_simulation(Y0, mu, alpha, annees) #stoker des tests
    
    return np.sum(tests)/nombre_tests #return probabilité





#%% probabilité theorique vs simulation


Y0 = np.linspace(5, 25, 51) * 1e6 #capital entre 5 - 25 millions
alpha = 0.1 * 1e6
 # le casino gagne 0.1 million per jour (100k)
mu = alpha * (1 - 0.05) #le casino perds 5% du alpha, ce qui gagnent des joueurs
annees = 5

proba_ruine = r(Y0,mu,alpha)
proba_ruine_simulation = np.zeros_like(Y0)

for i, y in enumerate(Y0):
    
    proba_ruine_simulation[i] = r_simulation(y, mu, alpha, annees, 2000)
    print(i)
    
#%% affichage probabilité theorique vs simulation

Y0_millions = Y0/1e6


plt.figure(figsize=(12, 6))
plt.plot(Y0_millions, proba_ruine, label='Théorique')
plt.plot(Y0_millions, proba_ruine_simulation, label='Simulation', marker='*', linestyle = '')
plt.title('r(y) vs y [alpha = 0.1 M et mu = 95%alpha]')
plt.xlabel('Capital Initial (y) - Millions')
plt.ylabel('r(y) - P(ruiné)')
plt.legend()
plt.grid(True)
plt.show()




#%% version 2 simulation pour montrer alpha>mu, alpha = mu et alpha<mu

def casino_simulation(Y0,mu,alpha,annees):
    
    lambda_gain = 1 / mu  # Paramètre lambda pour la loi exponentielle des gains (variable gamma)
    lambda_temps_gains = 1 # Paramètre lambda pour la loi exponentielle des temps de gains
    #------------------------------------------------------------------------------------------
    
    ### Simuler un processus de Poisson pour les instants de gain
    
    temps_total = annees * 365 # Durée totale de la simulation
    #pour s'assurer que la somme cumulée est toujours supérieure au temps de simulation
    size_gains = (lambda_temps_gains) * (temps_total + 100) 
    #somme chi_i -> T_i
    temps_gains = np.random.exponential(scale = 1/lambda_temps_gains, size = size_gains).cumsum() 
    temps_gains = temps_gains[temps_gains <= temps_total]
    
    #------------------------------------------------------------------------------------------
    
    ### Simuler les gains des joueurs
    gains_joueurs = np.random.exponential(scale=1/lambda_gain, size=len(temps_gains)) #variable X_i
    
    #------------------------------------------------------------------------------------------
    
    ### Calculer le capital du casino au fil du temps
    
    temps = np.arange(temps_total + 1)
    capital_casino = np.zeros_like(temps) #vecteur composé de zéros pour stocker l'évolution du capital
    index_gains = 0
    somme_gains = 0
    
    for i, t in enumerate(temps):
        
        while index_gains < len(temps_gains) and temps_gains[index_gains] <= t:
            # cond 1: Éviter un Débordement de l'Index
            # cond 2: Condition N_t definition
            somme_gains += gains_joueurs[index_gains]
            index_gains += 1
            
        #somme_gains accumulés jusqu'à présent t
        capital_casino[i] = Y0 + alpha * t - somme_gains
        
        if capital_casino[i] < 0:
            break
    
    return np.array(temps), np.array(capital_casino)

#%%

annees = 10
Y0 = 8 * 1e6
alpha = 0.1 * 1e6
mu = [ alpha * (1 + 0.05), alpha * (1 + 0.0), alpha * (1 - 0.05) ]

temps1, casino1 = casino_simulation(Y0, mu[0], alpha, annees) 
temps2, casino2 = casino_simulation(Y0, mu[1], alpha, annees)
temps3, casino3 = casino_simulation(Y0, mu[2], alpha, annees)

x_data = np.array([temps1, temps2, temps3])
y_data = np.array([casino1, casino2, casino3])
x_data = x_data/365
y_data = y_data/1e6
titles = ['Simulation du Capital du Casino - alpha (5% <) mu', 'Simulation du Capital du Casino - mu = alpha', 'Simulation du Capital du Casino - alpha (5% >) mu']  # Titres correspondants

# Nombre de lignes et de colonnes pour les sous-graphiques
n_rows = 3
n_cols = 1

# Création de la figure et des axes
fig, axs = plt.subplots(n_rows, n_cols, figsize=(8 * n_cols, 4 * n_rows))

if n_rows == 1 or n_cols == 1:
    axs = np.array(axs).reshape(n_rows, n_cols)

# Parcourir et créer chaque sous-graphique
for i in range(n_rows):
    for j in range(n_cols):
        index = i * n_cols + j
        if index < len(y_data):
            axs[i, j].plot(x_data[index], y_data[index])
            axs[i, j].set_title(titles[index])
            axs[i, j].set_xlabel('temps (années)')
            axs[i, j].set_ylabel('Capital Casino (millions)')

# Ajuster l'espacement
plt.tight_layout()
plt.show()



#%% Determination de C, generation de la distribution de Hn

def F_chapeau_empirique(echantillon, x):
    #Calculer la Fonction Répartition empirique à la valeur x.
    return np.mean(echantillon <= x)

def F_exponentielle(x, gamma):
    #Calculer la Fonction Répartition d'une distribution exponentielle à la valeur x
    return 1 - np.exp(-gamma * x)

def calculer_Hn(echantillon, gamma_chapeau, resolution):
    #Calculer la statistique H_n pour l'échantillon donné et le gamma estimé
    
    valeurs_y = np.linspace(0, max(echantillon), resolution)
    
    valeurs_Hn = np.zeros(resolution)
    
    for i, y in enumerate(valeurs_y):
        valeurs_Hn[i] = abs(F_chapeau_empirique(echantillon, y) - F_exponentielle(y, gamma_chapeau))
    
    return max(valeurs_Hn)

# Simulation d'un échantillon exponentiel
alpha = 0.1 * 1e6
mu = alpha * (1 - 0.05)
n = 1000  # taille de l'échantillon
gamma_vrai = 1/mu  # paramètre réel de la distribution exponentielle
echantillon = np.random.exponential(1/gamma_vrai, n)

# Estimation de gamma à partir de l'échantillon
gamma_chapeau = n / np.sum(echantillon)

# Calcul de H_n (pour un seul test)
resolution = 500
Hn = calculer_Hn(echantillon, gamma_chapeau, resolution)

# Affichage des résultats
print(f'Hn = {Hn}')
print(f'gamma chapeau = {gamma_chapeau}')
print(f'gamma_vrai = {gamma_vrai}')


#%%
#pour plusieurs tests

# Nombre de simulations
nombre_simulations = 5000
n = 1000  # taille de l'échantillon


# Simulation des valeurs de H_n sous l'hypothèse que la distribution est exponentielle
simulations_Hn = np.zeros(nombre_simulations)

for i in range(nombre_simulations):
    
    # Simulation d'un échantillon de la distribution exponentielle 
    echantillon_simule = np.random.exponential(scale = 1/gamma_vrai, size = n)
    
    # Calcul de la fonction répartition empirique de l'échantillon simulé
    gamma_chapeau_simule = n / np.sum(echantillon_simule)
    Hn_simule = calculer_Hn(echantillon_simule, gamma_chapeau_simule, resolution)
    simulations_Hn[i] = Hn_simule
    print(i)


#%%

# Définition du niveau de signification
alpha = 0.05

# Détermination du seuil C comme le quantile (1-alpha) des valeurs simulées de H_n 
C = np.quantile(simulations_Hn, 1 - alpha)
print(C)

shape, loc, scale = stats.gamma.fit(simulations_Hn, floc=0)
C_gamma = stats.gamma.ppf(1-alpha, shape, loc=0, scale=scale)
print(C_gamma)



#%% Affichage de la distribution de H_n

x_valeurs = np.linspace(0, max(simulations_Hn), 1000)
distribution_gamma = stats.gamma.pdf(x_valeurs, shape, loc=0, scale=scale)


plt.hist(simulations_Hn, density=True)
plt.plot(x_valeurs, distribution_gamma, label = f"loi Gamma k = {shape: .2f}, theta = {scale: .4f}")
plt.title("Distribution de H_n pour des échantillons exponentiels simulés")
plt.xlabel("Valeurs de H_n")
plt.ylabel("Fréquence")
plt.legend(loc='upper right')
plt.show()

#%% tendance vers zéro de Hn quand n vers le infini


# Paramètres de la simulation
tailles_echantillon = np.arange(100, 5001, 100)  # Tailles d'échantillon croissantes de 100 à 5000
# Simulation pour différentes tailles d'échantillon

valeurs_Hn = np.zeros(len(tailles_echantillon))
for i, n in enumerate(tailles_echantillon):
    
    # Simulation d'un échantillon de la distribution exponentielle 
    echantillon_simule = np.random.exponential(scale = 1/gamma_vrai, size = n)
    
    # Calcul de la fonction répartition empirique de l'échantillon simulé
    gamma_chapeau_simule = n / np.sum(echantillon_simule)
    
    Hn_simule = calculer_Hn(echantillon_simule, gamma_chapeau_simule, resolution)
    valeurs_Hn[i] = Hn_simule
    

# Affichage du graphique
plt.plot(tailles_echantillon, valeurs_Hn, marker='*', linestyle = '-', color = 'g')
plt.title("Convergence de H_n avec la taille de l'échantillon")
plt.xlabel("Taille de l'échantillon (n)")
plt.ylabel("Valeur de H_n")
plt.grid(True)
plt.show()
