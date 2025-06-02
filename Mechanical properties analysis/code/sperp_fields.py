import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.interpolate import griddata

# load temporel interpolated reshaped exprimental data
def load_exp_data(i):
    file_path = f"E:/2023_2025/M2/Projet/Abaqus_simulation/interpolated_exp_data/image_{i:03}.csv"
        # check the existence of files
    if not os.path.exists(file_path):
        print(f"File not found: {file_path}")
        return None
    
    data = np.loadtxt(file_path, delimiter=",", skiprows=1)

    return data

def load_num_data(i):
    file_path =  f"E:/2023_2025/M2/Projet/Abaqus_simulation/numeric_data/champs_deformation/deformation_data_{i}.csv"
        # check the existence of files
    if not os.path.exists(file_path):
        print(f"File not found: {file_path}")
        return None
    
    data = np.loadtxt(file_path, delimiter=",", skiprows=1)

    return data

if __name__ == "__main__":    # exp to num
    # load initial geometric data
    #data_num = load_num_data(0)
    N = 285
    x_num, y_num = load_num_data(N)[:,5], load_num_data(N)[:,6]
    x_exp, y_exp = load_exp_data(N)[:,1], load_exp_data(N)[:,2]

    #print("x exp", np.shape(x_exp), "y exp", np.shape(y_exp))

    # interpolate exp data to num data
    # itp_y_exp = griddata((x_exp, y_exp), y_exp, (x_num, y_num), method='cubic')
    # itp_x_exp = griddata((y_exp, x_exp), x_exp, (y_num, x_num), method='cubic')

    # interpolate num data to exp data
    itp_y_num = griddata((x_num, y_num), y_num, (x_exp, y_exp), method='cubic')
    itp_x_num = griddata((y_num, x_num), x_num, (y_exp, x_exp), method='cubic')

    # plot the initial superposition
    plt.figure()
    #plt.scatter(data_num[:,5], data_num[:,6], color = 'g')
    #plt.scatter(x_num, y_num, color = 'red', alpha=1, label = 'num_data')

    # interpolated plot
    #plt.scatter(itp_x_exp, itp_y_exp, color = "black", alpha=0.3, label="interpolated_exp_data")
    plt.scatter(itp_x_num, itp_y_num, color = "black", alpha=1, label="interpolated_num_data")

    plt.scatter(x_exp, y_exp, color = "red", alpha=0.05, label = 'exp_data')

    plt.xlabel("X")
    plt.ylabel("Y")
    plt.title(f"Superposition at t= {N*144.04/290}s")
    plt.legend()
    plt.show()

# if __name__ == "__main__":    # pour tout le long du temps
     
#      for i in range(0, 289):
#         x_num, y_num = load_num_data(i)[:,5], load_num_data(i)[:,6]
#         x_exp, y_exp = load_exp_data(i)[:,1], load_exp_data(i)[:,2]

#         # interpolate exp data to num data
#         itp_y_exp = griddata((x_exp, y_exp), y_exp, (x_num, y_num), method='cubic')
#         itp_x_exp = griddata((y_exp, x_exp), x_exp, (y_num, x_num), method='cubic')

#         # interpolate num data to exp data
#         itp_y_num = griddata((x_num, y_num), y_num, (x_exp, y_exp), method='cubic')
#         itp_x_num = griddata((y_num, x_num), x_num, (y_exp, x_exp), method='cubic')

#         plt.figure()
#         plt.scatter(x_exp, y_exp, color = "black", alpha=0.9)
#         # interpolated plot
#         plt.scatter(itp_x_num, itp_y_num, color = "b", alpha=0.1)
#         plt.xlabel("X")
#         plt.ylabel("Y")
#         plt.title(f"Superposition at t= {i*144.04/290}s")
#         plt.legend()
#         plt.show()