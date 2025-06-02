from abaqus import *
from abaqusConstants import *
from odbAccess import openOdb
import numpy as np


def get_abaqus_time(odb_name):
    # open ODB file
    odb = openOdb(path=odb_name)

    # get the step1 and its data
    last_step = list(odb.steps.keys())[-1] # setp1
    step1_frames = odb.steps[last_step].frames # data of step1

    abqtime = []

    for frame in step1_frames:
        abqtime.append(frame.frameValue)

    return abqtime

if __name__ == '__main__':
    time_abaqus = get_abaqus_time('E:/2023_2025/M2/Projet/Abaqus_simulation/traction.odb')
    output_path = "E:/2023_2025/M2/Projet/Abaqus_simulation/time_abaqus.csv"
    np.savetxt(output_path, time_abaqus, fmt="%.3f", delimiter=",", header="time", comments="")
    print(time_abaqus)