from abaqus import *
from abaqusConstants import *
from odbAccess import openOdb
#import matplotlib.pyplot as plt
import numpy as np


# open ODB file
odb_path = r'E:\\2023_2025\\M2\\Projet\\Abaqus_simulation\\traction.odb'
odb = openOdb(path=odb_path)
print(f"ODB opened: {odb}")

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

FD = np.zeros((len(U2),2))
FD[:,0] = U2
FD[:,1] = RF2
# plt.figure()
# plt.plot(U2, RF2, label = 'FD')
# plt.title('force-displacemnt')
# plt.xlabel('displacemnt')
# plt.ylabel('force')
# plt.legend()
# plt.show()


output_path = r"E:\\2023_2025\\M2\\Projet\\Abaqus_simulation\\numeric_data\\force_displacement\\force_displacement.csv"

# save as a csv file
np.savetxt(output_path, FD, fmt="%.3f", delimiter=",", 
           header="displacement, force", comments="")
# FX and FY represent the deformation field
print(f"Data saved to {output_path}")


# from abaqus import *
# from abaqusConstants import *
# from odbAccess import openOdb

# odb_path = r'E:\\2023_2025\\M2\\Projet\Abaqus_simulation\\traction.odb'
# odb = openOdb(path=odb_path)
# print(f"ODB opened: {odb}")
# assembly = odb.rootAssembly

# node_set_names = [node_set.name for node_set in assembly.nodeSets.values()]
# print("Available node set names:", node_set_names)

# odb.close()