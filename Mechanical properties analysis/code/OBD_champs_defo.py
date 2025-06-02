from abaqus import *
from abaqusConstants import *
from odbAccess import openOdb
import job
import os
import shutil
import numpy as np

# extraire les champs numeric pour tout les neuds au cours du temps

# open ODB file
odb_path = 'E:/2023_2025/M2/Projet/Abaqus_simulation/traction.odb'
odb = openOdb(path=odb_path)
print(f"ODB opened: {odb}")

# get the step1 and the initial time 
last_step = odb.steps.keys()[-1] # setp1
first_frame = odb.steps[last_step].frames[0] # initial time

num_frames = len(odb.steps[last_step].frames)
print(f"Total frames: {num_frames}")


for i in range(num_frames):

    frame = odb.steps[last_step].frames[i]
    # get force and displacement data
    displacement_field = frame.fieldOutputs['U']

    # # get instance
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
        ydisplacements.append(value.data[1])

    # Create a numpy array to store the displacement field
    disp_field = np.zeros((len(node_labels), 7))
    disp_field[:,0] = node_labels
    disp_field[:,1] = x_coords
    disp_field[:,2] = y_coords
    disp_field[:,3] = xdisplacements 
    disp_field[:,4] = ydisplacements
    # defomed coord.
    disp_field[:,5] = disp_field[:,1] + disp_field[:,3] 
    disp_field[:,6] = disp_field[:,2] + disp_field[:,4]



    # if the file already exists, automatically back it up
    output_path = f"E:/2023_2025/M2/Projet/Abaqus_simulation/numeric_data/champs_deformation/deformation_data_{i}.csv"

    if not os.path.exists(output_path):
        base_dir = "E:/2023_2025/M2/Projet/Abaqus_simulation/numeric_data/champs_deformation"
        os.makedirs(base_dir, exist_ok=True)

    # save as a csv file
    np.savetxt(output_path, disp_field, fmt="%.3f", delimiter=",", 
            header="NodeLabel,X,Y,UX,UY, defmX, defmY", comments="")
    # FX and FY represent the deformation field
    print(f"Data saved to {output_path}")

# # close ODB field
odb.close()