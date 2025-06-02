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
    index_sorted = np.argsort(data[:, 6])
    data = data[index_sorted]

    index = np.arange(len(data))
    rdata = np.column_stack((index, data[:,5], data[:,6]))
    return rdata

# load initial geometric data
data_num_ini = load_num_data(0)
data_exp_ini = load_exp_data(0)

# center point of geometry (float)
y_cent_num = (max(data_num_ini[:,2]) + min(data_num_ini[:,2])) / 2
y_cent_exp = (max(data_exp_ini[:,2]) + min(data_exp_ini[:,2])) / 2

# difine y limits
y_max_exp,  y_min_exp = y_cent_exp + 15,  y_cent_exp -15 
y_max_num,  y_min_num = y_cent_num + 15,  y_cent_num -15

# define geometric domain for observation
y_exp_filter = (data_exp_ini[:,2] >= y_min_exp) & (data_exp_ini[:,2] <= y_max_exp)
filtered_data_exp = data_exp_ini[y_exp_filter]
x_exp = filtered_data_exp[:,1]
y_exp = filtered_data_exp[:,2] # updata y data

y_num_filter = (data_num_ini[:,2] >= y_min_num) & (data_num_ini[:,2] <= y_max_num)
filtered_data_num = data_num_ini[y_num_filter]
x_num = filtered_data_num[:,1]
y_num = filtered_data_num[:,2] # updata y data

# get the coordinate index
###########################################################
index_num = data_num_ini[:,0]
index_exp = data_exp_ini[:,0] - 1 

filtered_index_num = index_num[y_num_filter] # Filtered indices
filtered_index_exp = index_exp[y_exp_filter] # Filtered indices
###########################################################
i_min_num, i_max_num = int(min(filtered_index_num)), int(max(filtered_index_num))     #  minimum index and maximum index for useful numeric zone
i_min_exp, i_max_exp = int(min(filtered_index_exp)), int(max(filtered_index_exp))     #  minimum index and maximum index for useful experimental zone
###########################################################


def superp(i): # num to exp
    # just use the useful zone xy coord. 
    data_num = load_num_data(i)
    data_exp = load_exp_data(i)

    # Use boolean masks to select data
    num_mask = (index_num >= i_min_num) & (index_num <= i_max_num)
    exp_mask = (index_exp >= i_min_exp) & (index_exp <= i_max_exp)

    x_num = data_num[num_mask, 1]
    y_num = data_num[num_mask, 2]
    x_exp = data_exp[exp_mask, 1]
    y_exp = data_exp[exp_mask, 2]

    points_num = np.column_stack((x_num, y_num))  
    points_exp = np.column_stack((x_exp, y_exp)) 

    #interpolate num data to exp data
    itp_x = griddata(points_num, x_num, points_exp, method='cubic')
    itp_y = griddata(points_num, y_num, points_exp, method='cubic')

    itp_data = np.column_stack((itp_x, itp_y))
    data_exp = np.column_stack((x_exp, y_exp))

    return itp_data, data_exp

if __name__ == '__main__':  # num to exp
    i = 289
    itp_data, exp_data = superp(i)
    itp_x_num, itp_y_num = itp_data[:, 0], itp_data[:, 1]
    x_exp, y_exp = exp_data[:, 0], exp_data[:, 1]
    
    data_num = load_num_data(i)
    num_mask = (index_num >= i_min_num) & (index_num <= i_max_num)
    x_num = data_num[num_mask, 1]
    y_num = data_num[num_mask, 2]

    # plot the i superposition
    plt.figure()
    plt.scatter(x_num, y_num, color = 'black', label='num field')
    plt.scatter(itp_x_num, itp_y_num, color = 'r', label= 'interpolated num field')
    plt.scatter(x_exp, y_exp, color = "b", alpha=0.35, label='exp field')
    plt.xlabel("X")
    plt.ylabel("Y")
    plt.title(f"Superposition at t= {i*144.04/290}s")
    plt.legend()
    plt.show()

