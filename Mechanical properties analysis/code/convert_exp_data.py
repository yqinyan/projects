import numpy as np
import os
import matplotlib.pyplot as plt

# initial orgin and ratio number
ORIGIN_X = 354
ORIGIN_Y = 1522
PIXEL_TO_MM = 0.0694078947368421

def load_and_process_data(i):
    # file path 
    file_path = f"E:/2023_2025/M2/Projet/projet_data/Essai_traction/data_traction/alu_7075/image-00000{i:03}_0.csv"
    
    # check the existence of files
    if not os.path.exists(file_path):
        print(f"File not found: {file_path}")
        return None
    
    # read 8 and 9 colonnes
    data = np.loadtxt(file_path, delimiter=",", skiprows=1, usecols=(7, 8))
    
    # change coordinate system
    data[:, 0] -= ORIGIN_X  # X tanslation
    data[:, 1] -= ORIGIN_Y  # Y translation
    
    # pixels to mm
    data *= PIXEL_TO_MM
    
    coord_label = np.linspace(1, len(data[:,0]), len(data[:,0]))
    def_data = np.zeros((len(data[:,0]), 3))
    def_data[:,0] = coord_label
    def_data[:,1] = data[:,0]
    def_data[:,2] = data[:,1]
    
    return data

def save_data(data, i):
    # save new data
    output_path = f"E:/2023_2025/M2/Projet/Abaqus_simulation/experimental_data/exp_deformation_{i:03}.csv"
    np.savetxt(output_path, data, fmt="%.3f", delimiter=",", header="X, Y", comments="")
    print(f"Data saved to {output_path}")

# main program
if __name__ == "__main__":
    data_exp = []  # initialize the list
    for i in range(0, 289):
        data = load_and_process_data(i)
        if data is not None:
            data_exp.append(data)  # save data in the list
            save_data(data, i)


data_0 = load_and_process_data(0)
data_288 = load_and_process_data(288)
plt.figure()
plt.scatter(data_0[:,0], data_0[:,1], color = 'r')
plt.scatter(data_288[:,0], data_288[:,1], color = 'b')
plt.xlabel("X")
plt.ylabel("Y")
plt.legend()
plt.show()