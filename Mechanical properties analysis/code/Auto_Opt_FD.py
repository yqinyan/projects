from abaqus import *
from abaqusConstants import *
from odbAccess import openOdb
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from scipy.optimize import differential_evolution
import numpy as np
import os


root_path = r"E:\\2023_2025\\M2\\Projet"
# CAE file path
cae_file_path = os.path.join(root_path, "Abaqus_simulation","simulation_traction.cae")
# open CAE file
openMdb(pathName = cae_file_path)


def abaqus_FD(JC_params):
    # ODB file path
    odb_path = os.path.join(root_path, "Abaqus_simulation", "traction.odb")

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

    # open ODB file
    odb = openOdb(path=odb_path)
    print(f"ODB_traction_yy opened: {odb}")

    # get the step1 and its data
    last_step = list(odb.steps.keys())[-1] # setp1
    step1_frames = odb.steps[last_step].frames # data of step1

    # difine initial lists to save displacement field
    U2 = []
    RF2 = []

    assembly = odb.rootAssembly
    region_name = 'SET-RP'  # replace set name
    region = assembly.nodeSets[region_name]


    for frame in step1_frames:
        # get reaction force RF and displacement U data
        force = 0.0
        disp = 0.0

        # get RF2
        RF2_field = frame.fieldOutputs['RF'] # get reaction force RF data
        RF2_values = RF2_field.getSubset(region=region).values # get reaction force RF values
        for val in RF2_values:
            force += val.data[1]  # index of RF: [0]=X, [1]=Y, [2]=Z

        # get U2
        U2_field = frame.fieldOutputs['U']  # get displacement U data
        U2_values = U2_field.getSubset(region=region).values # get displacement U values
        for val in U2_values:
            disp += val.data[1]  # index of U same like before

        # save present Frame data
        RF2.append(force)
        U2.append(disp)

    # close ODB field..
    odb.close()

    # FD = np.zeros((len(U2),2))
    # FD[:,0] = U2
    # FD[:,1] = RF2
    return np.array(RF2), np.array(U2)#, FD


###############################################################################################

###############################################################################################


# load experimental data
exp_path = os.path.join(root_path, "FD_exp.csv")
data_exp = np.genfromtxt(exp_path, delimiter=",", skip_header=1) 

# define numeric and experimental force and displacement

force_exp, displacement_exp = data_exp[:, 1], data_exp[:, 0] 

# define interpolate domain 
N = len(displacement_exp)*1 # number of usual points
common_displacement = np.linspace(min(displacement_exp), max(displacement_exp), N)


# difine cost function
def cost_fct(JC_params):

    try:
        force_num, displacement_num = abaqus_FD(JC_params)
    except Exception as e:
        print(f"Error during abaqus_FD: {e}")
        return float('inf')  # return infinit value

    # check NaN value
    if np.any(np.isnan(force_num)) or np.any(np.isnan(displacement_num)):
        print("NaN encountered in numerical results")
        return float('inf')
    
    #force_num, displacement_num = abaqus_FD(JC_params)

    exp_force_interp = np.interp(common_displacement, displacement_exp, force_exp)
    num_force_interp = np.interp(common_displacement, displacement_num, force_num) 

    # normalized error 
    denominator = np.sum(exp_force_interp ** 2)
    if denominator == 0:
        raise ValueError("Experimental force data sum is zero, normalization is invalid.")
    err = np.sum((exp_force_interp - num_force_interp) ** 2) / denominator
    return err 


JCparameters = []

JCparams = [450, 270, 0.3]
#JCparams = [71000, 0.3, 450, 270, 0.3]

# optimize iterations
nb_iter = 1
for i in range(1, nb_iter+1):
    print(f" A = {JCparams[0]}, B = {JCparams[1]}, n = {JCparams[2]}")
    #print(f"E = {JCparams[0]}, v = {JCparams[1]},  A = {JCparams[2]}, B = {JCparams[3]}, n = {JCparams[4]}")

    # optimization
    
    #bound = [(50000,71000), (0.15, 0.39), (300, 650), (200, 600), (0.1, 0.7)]
    bound = [(300, 650), (200, 600), (0.1, 0.7)]
    result = minimize(cost_fct, JCparams, method='Nelder-Mead', bounds=bound, options={'maxiter': 1})#, 'disp': True})
    #result = differential_evolution(cost_fct,  bounds=bound, maxiter=1)
    print("optimized Johnson Cook model:", result) 

    # updata parameters
    JCparams = result.x

    JCparameters.append(list(JCparams) + [result.fun])

print(f"Johnson Cook parameters and errors{JCparameters}")

params_path = os.path.join(root_path, "Abaqus_simulation", "numeric_data", "JC_parameters.csv")

# save as a csv file
np.savetxt(params_path, JCparameters, fmt="%.3f", delimiter=",", 
           header="Young's Module, Poisson coeff, A,B,n, error", comments="")
# FX and FY represent the deformation field
print(f"Data saved to {params_path}")

###############################################################################################

###############################################################################################

# force_num, displacement_num = abaqus_FD(JCparams)
# num_force_interp = np.interp(common_displacement, displacement_num, force_num)

# plt.figure(figsize=(10, 6))
# plt.plot(displacement_exp, force_exp, label="Experimental", color="blue", linestyle='--')
# plt.plot(num_force_interp, force_num, label="Numerical", color="red")
# plt.xlabel("Displacement")
# plt.ylabel("Force")
# plt.legend()
# plt.title(f"Comparison of Experimental and Numerical Results (Iter {i})")
# plt.grid()
# plt.savefig(os.path.join(root_path, "Abaqus_simulation", "plots", f"comparison_iter_{i}.png"))
# plt.close()