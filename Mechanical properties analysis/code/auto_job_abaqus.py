import os
from abaqus import *
from abaqusConstants import *

N_folder = 2
#  A_fit, B_fit, n_fit are optimized
A_fit = 500  
B_fit = 300  
n_fit = 0.3  

# set CAE file path
cae_file_path = 'E:/2023_2025/M2/Projet/Abaqus_simulation/simulation_traction.cae'
# open CAE file
openMdb(pathName = cae_file_path)

#  Modifier les parametres du material
# choose the modele
model_name = "Model-1"  
material_name = "Alu-7075"  

# open Modele and Materiele
model = mdb.models[model_name]
material = model.materials[material_name]


# update Johnson-Cook parameters
material.Plastic(hardening=JOHNSON_COOK, table=((A_fit, B_fit, n_fit),))


folder_path = []
for i in range(1, N_folder+1):
    folder_name = f'Job_opt{i}'
    folder_path = f'E:/2023_2025/M2/Projet/Abaqus_simulation/Jobs_opt/{folder_name}'

    # creat a folder if it doesn't exist
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)
        print(f"folder '{folder_path}' created")
    else:
        print(f"folder '{folder_path}' exists before")
    

    # creer ou mise a jour le job
    job_name = f"Optimized_Job{i}" 

    # si job optimise existe, remplace le
    if job_name in mdb.jobs.keys():
        del mdb.jobs[job_name]


    # set Job path
    os.chdir(folder_path)

    # creer nouveau job
    mdb.Job(name=job_name, model=model_name)

    # submit job
    mdb.jobs[job_name].writeInput()
    mdb.jobs[job_name].submit(consistencyChecking=OFF)
    
    # wait for the job complet
    mdb.jobs[job_name].waitForCompletion()

    cae_name = f"simulation_traction{i}"
    # savgarder du modele
    mdb.saveAs(os.path.join(folder_path, cae_name))