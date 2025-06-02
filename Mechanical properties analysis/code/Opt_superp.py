from abaqus import *
from abaqusConstants import *
from odbAccess import openOdb
import job

import numpy as np
#import matplotlib.pyplot as plt
import os
import gc
from scipy.interpolate import griddata
from scipy.optimize import minimize


# load temporel interpolated reshaped exprimental data
def load_exp_data(i):
    file_path = f"E:/2023_2025/M2/Projet/Abaqus_simulation/interpolated_exp_data/image_{i:03}.csv"
        # check the existence of files
    if not os.path.exists(file_path):
        print(f"File not found: {file_path}")
        return None
    
    data = np.loadtxt(file_path, delimiter=",", skiprows=1)

    return data

# load abaqus deformation data
def load_num_data(i):
    odb_path = 'E:/2023_2025/M2/Projet/Abaqus_simulation/traction.odb'
    odb = openOdb(path=odb_path)
    print(f"ODB opened: {odb}")

    # get the step1 and the initial time 
    last_step = odb.steps.keys()[-1] # setp1
    frames = odb.steps[last_step].frames[i] # for each time step

    # get force and displacement data
    displacement_field = frames.fieldOutputs['U']

    # get instance
    instance_name = list(odb.rootAssembly.instances.keys())[0]
    instance = odb.rootAssembly.instances[instance_name]

    # difine initial lists to save displacement field
    node_labels = []
    x_coords = []
    y_coords = []
    xdisplacements = []
    ydisplacements = []

    # append datas
    for value in displacement_field.values:

        node_label = value.nodeLabel  # node number
        node_labels.append(node_label) 
        node = instance.nodes[node_label - 1] # get instance nodes

        x_coords.append(node.coordinates[0])  # X coord
        y_coords.append(node.coordinates[1])  # Y coord
        xdisplacements.append(value.data[0])  # X displacement
        ydisplacements.append(value.data[1])  # Y displacement

    # close ODB field
    odb.close()

    # Create a numpy array to store the displacement field
    disp_field = np.zeros((len(node_labels), 7))
    disp_field[:,0] = node_labels
    disp_field[:,1] = x_coords
    disp_field[:,2] = y_coords
    disp_field[:,3] = xdisplacements 
    disp_field[:,4] = ydisplacements
    # displacement field
    disp_field[:,5] = disp_field[:,1] + disp_field[:,3] 
    disp_field[:,6] = disp_field[:,2] + disp_field[:,4]

    index_sorted = np.argsort(disp_field[:, 6])
    filtered_data = disp_field[index_sorted]
    index = np.arange(len(filtered_data))
    data = np.column_stack((index, filtered_data[:,5], filtered_data[:,6]))

    return data

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

# interpolate num data to exp data
itp_y_num_ini = griddata((x_num, y_num), y_num, (x_exp, y_exp), method='cubic')
itp_x_num_ini = griddata((y_num, x_num), x_num, (y_exp, x_exp), method='cubic')

# img_path = r"E:\\2023_2025\\M2\\Projet\\Abaqus_simulation\\images\\superposition_initial.png"
# # plot the initial superposition
# plt.figure()
# plt.scatter(itp_x_num_ini, itp_y_num_ini, color = 'r')
# plt.scatter(x_exp, y_exp, color = "y")
# plt.xlabel("X")
# plt.ylabel("Y")
# plt.title("Superposition at t=0")
# plt.legend()
# plt.savefig(img_path)
# plt.close()
# print(f"plot image is saved at: {img_path}")

root_path = r"E:\\2023_2025\\M2\\Projet"
def abaqus_sim(i, JC_params):

    # CAE file path
    cae_file_path = os.path.join(root_path, "Abaqus_simulation","simulation_traction.cae")
    # open CAE file
    openMdb(pathName = cae_file_path)

    # get elastic and Johnson-Cook paramamters
  
    #E, v, A, B, n = JC_params
    A, B, n = JC_params

    # updata paramaters
    model_name = "Model-1"  
    material_name = "Alu-7075"  

    # open Modele and Materiele
    model = mdb.models[model_name]
    material = model.materials[material_name]
    #material.Elastic(table=((E, v),))
    material.Plastic(hardening=JOHNSON_COOK, table=((A, B, n),))

    # save ODB file and submit the job
    os.chdir(r"E:\\2023_2025\\M2\\Projet\\Abaqus_simulation")
    job_name = 'traction'
    mdb.Job(name=job_name, model=model_name)
    mdb.jobs[job_name].submit(consistencyChecking=OFF)
    mdb.jobs[job_name].waitForCompletion()

    cae_name = "simulation_traction"
    # savgarder du modele
    mdb.saveAs(os.path.join(root_path, "Abaqus_simulation", cae_name))

    data = load_exp_data(i)

    return data

def superp(i, JC_params):
    # just use the useful zone xy coord. 
    data_num = load_num_data(i)
    data_exp = abaqus_sim(i,JC_params)

    # filteration of data
    num_mask = (index_num >= i_min_num) & (index_num <= i_max_num)
    exp_mask = (index_exp >= i_min_exp) & (index_exp <= i_max_exp)
    x_num = data_num[num_mask, 1]
    y_num = data_num[num_mask, 2]
    x_exp = data_exp[exp_mask, 1]
    y_exp = data_exp[exp_mask, 2]

    points_num = np.column_stack((x_num, y_num))  
    points_exp = np.column_stack((x_exp, y_exp)) 

    data_exp = np.column_stack((x_exp, y_exp))
    data_num = np.column_stack((x_num, y_num))

    #interpolate num data to exp data
    # itp_x = griddata(points_num, x_num, points_exp, method='cubic')
    # itp_y = griddata(points_num, y_num, points_exp, method='cubic')

    #interpolate exp data to num data
    itp_x = griddata(points_exp, x_exp, points_num, method='linear')
    itp_y = griddata(points_exp, y_exp, points_num, method='linear')

    # nitp_data = np.column_stack((itp_x, itp_y))
    eitp_data = np.column_stack((itp_x, itp_y))

    # return nitp_data, data_exp
    return eitp_data, data_num

def cost_fct(JC_params, i):

    itp_data, data_exp = superp(i, JC_params)

    denominator = np.sum(data_exp ** 2)
    if denominator == 0:
        raise ValueError("Experimental force data sum is zero, normalization is invalid.")
    error = np.sqrt(np.sum((data_exp - itp_data)** 2)) / np.sqrt(denominator)

    return error


if __name__ == "__main__":

    JCparam_err = []
    JC_params = [450, 270, 0.3]  # initialize the list

    for i in range(122, 299):  # time of platic deformation
        try:
            result = minimize(cost_fct, JC_params, method='Nelder-Mead', options={'disp': False, 'maxiter': 50}, args=(i,))
            print("optimized Johnson Cook parameters:", result.x) 

            if result.success:

                # updata parameters
                JCparams = result.x
                JCparam_err.append(list(JCparams) + [result.fun])
                print(f"Johnson Cook parameters and errors{JCparam_err}")
        
            # Explicit memory cleanup
            del result
            gc.collect()
        
        except Exception as e:
            print(f"Optimization failed at step {i}: {str(e)}")
            continue