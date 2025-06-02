import matplotlib.pyplot as plt
import numpy as np

# nb de processes
processes = np.array([2, 4, 8, 16, 32, 64, 128, 256])

# Partie1
# Temps de calcul
Calcul_temps1 = np.array([
    6.480969488620758E-003,
    4.164213314652443E-003,
    2.727994695305824E-003,
    1.397788058966398E-003,
    8.113912772387266E-004,
    6.913404213264585E-004,
    5.526696157176048E-003,
    1.255465584108606E-003
])

# Temps de communication
Comm_temps1 = np.array([
    1.303447782993317E-002,
    7.120549678802490E-003,
    4.384740255773067E-003,
    2.793956082314253E-003,
    1.818863442167640E-003,
    2.023664535954595E-003,
    7.846751541364938E-003,
    9.569121140521020E-002
])

# Temps total
temps_total1 = Comm_temps1 + Calcul_temps1

plt.figure()
plt.plot(processes, Calcul_temps1, marker='o', linestyle='-', label='Computing Time')
plt.plot(processes, Comm_temps1, marker='s', linestyle='--', label='Communication Time')
plt.plot(processes, temps_total1, marker='s', linestyle='--', label='Communication Time')
plt.xscale('log', base=2)
plt.xticks(processes, labels=[str(p) for p in processes])
plt.yscale('log')
plt.xlabel("Nombre de Processes")
plt.ylabel("Temps (secondes)")
plt.title("Scaling Analyse : Computation vs Communication")
plt.legend()
plt.grid(True, which="both", linestyle="--", linewidth=0.5)
plt.show()


# Partie2
# Temps de calcul
Calcul_temps2 = np.array([
    169.152249292587,
    85.2271906737005,
    43.2693037780700,
    22.2915858876659,
    11.4479497986613,
    6.19421366194729,
    3.47820109885652,
    1.99176577385515
])

# Temps de communication
Comm_temps2 = np.array([
    169.152249292587,
    85.2271906737005,
    43.2693037780700,
    22.2915858876659,
    11.4479497986613,
    6.19421366194729,
    3.47820109885652,
    1.99176577385515
])

# Temps total
temps_total2 = Comm_temps2 + Calcul_temps2

plt.figure()
plt.plot(processes, Calcul_temps2, marker='o', linestyle='-', lw = 2, label='Computing Time')
plt.plot(processes, Comm_temps2, marker='s', linestyle='--', label='Communication Time')
plt.plot(processes, temps_total2, marker='s', linestyle='-', label='Communication Time')
plt.xscale('log', base=2)
plt.xticks(processes, labels=[str(p) for p in processes])
plt.yscale('log')
plt.xlabel("Nombre de Processes")
plt.ylabel("Temps (secondes)")
plt.title("Scaling Analyse : Computation vs Communication")
plt.legend()
plt.grid(True, which="both", linestyle="--", linewidth=0.5)
plt.show()
