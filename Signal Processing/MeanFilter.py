# -*- coding: utf-8 -*-
"""
Created on Mon Apr 29 22:06:20 2024

@author: qinyan
"""

from pydub import AudioSegment
import numpy as np

def moyen_filtre(s):
    N = len(s)
    if N <= 2:
        return s 
    result = np.zeros_like(s)
    
    result[1:-1] = (s[:-2] + s[1:-1] + s[2:]) / 3
   
    result[0] = (s[0] + s[1]) / 2
    result[-1] = (s[-1] + s[-2]) / 2
    return result

def apply_filter_to_audio(file_path, output_path):
    
    audio = AudioSegment.from_file(file_path)
    
    
    frame_rate = audio.frame_rate
    channels = audio.channels
    sample_width = audio.sample_width
    
    
    samples = np.array(audio.get_array_of_samples())
    if channels > 1:
        
        samples = samples.reshape(-1, channels)
        filtered_samples = np.vstack([moyen_filtre(samples[:, ch]) for ch in range(channels)]).T.flatten()
    else:
        
        filtered_samples = moyen_filtre(samples.astype(float))
    
    
    filtered_audio = audio._spawn(filtered_samples.astype(audio.array_type), overrides={
        "frame_rate": frame_rate,
        "sample_width": sample_width,
        "channels": channels
    })

    
    filtered_audio.export(output_path, format="mp3")


input_path = "E:\\2023_2024\\M1\\S2\\UE\\Traitement_du_signal\\presentation\\input_signal.mp3"
output_path = "E:\\2023_2024\\M1\\S2\\UE\\Traitement_du_signal\\presentation\\output_signal.mp3"


apply_filter_to_audio(input_path, output_path)


