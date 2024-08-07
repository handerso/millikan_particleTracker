#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 04/25/2023 for The Ohio State University, Department of Physics
Particle Tracker for Millikan Oil Drop experiment.
Please refer to Millikan_Experimental_Report.pdf
Credit to ricktjwong
@author: Hyun Wallace Anderson
"""

# Imports
import matplotlib as mpl # V
import numpy as np # V
import pandas as pd # V
import pims # V
import trackpy as tp # V
from slicerator import pipeline # V
import os
    
# Video preprocessing (convert to to grayscale)
@pipeline
def as_grey(frame):
    red = frame[:, :, 0]
    green = frame[:, :, 1]
    blue = frame[:, :, 2]
    return 0.2125 * red + 0.7154 * green + 0.0721 * blue

mpl.rc('figure',  figsize=(10, 5))
mpl.rc('image', cmap='gray')

def compute_traj(vid, filename, csv_name):
    
    # Open video for preprocessing
    frames = as_grey(vid)
    
    # Declare frames to analyze (+/- midpoint of video)
    midpoint = len(frames)/2
    start = 0
    stop = len(frames)
    
    # Segmentation parameters (tweak for finer particle identification)
    particle_diameter = 23
    min_brightness = 100
    gyration_brightness = 8.0
    py_engine = 'numba'
    windows_os = 1 # If on windows OS set to 1 to disable multiprocessing
    survival_threshold = 120 # number of frames particle needs to survive
    
    # Output particle identification figure for first frame 
    f = tp.locate(frames[start], diameter = particle_diameter, separation = particle_diameter*2, maxsize = gyration_brightness, minmass = min_brightness)
    tp.annotate(f, frames[start])
    
    # Continue if user is happy with identification
    #print("Continue with data collection? (y/n)")
    #continue_bool = input()
    #if(continue_bool == "y"):
            
    # Trajectory segmentation
    f = tp.batch(frames[start:stop], diameter = particle_diameter, separation = particle_diameter*1.5, maxsize = gyration_brightness, minmass = min_brightness, engine = py_engine, processes = windows_os)
    t = tp.link_df(f, search_range = 10, memory = 5)
    t1 = tp.filter_stubs(t, threshold = survival_threshold)
    # d = tp.compute_drift(tf, 30)
    # t1 = tp.subtract_drift(tf, d)
            
    # Compare the number of particles in the unfiltered and filtered data.
    print('Unfiltered counts:', t['particle'].nunique())
    print('Final filtered counts:', t1['particle'].nunique())
    
    # Plot filtered particle trajectories
    tp.plot_traj(t1, mpp = 4.47047253)
    
    # Track pixel displacement and calculate velocity
    data = []
    for item in set(t1.particle):
        sub = t1[t1.particle==item]
        dvx = np.diff(sub.x)
        dvy = np.diff(sub.y)
        for x, y, dx, dy, frame, mass, size, ecc, signal, raw_mass, ep in \
        zip(sub.x[:-1], sub.y[:-1], dvx, dvy, sub.frame[:-1], sub.mass[:-1], sub['size'][:-1], sub.ecc[:-1], sub.signal[:-1], sub.raw_mass[:-1], sub.ep[:-1]):
            data.append({'particle': item * 1000 + filename, # Creates unique particle id
                        'dx': dx,
                        'dy': dy,
                        'x': x,
                        'y': y,
                        'frame': frame,
                        'size': size,
                        'ecc': ecc,
                        'signal': signal,
                        'mass': mass,
                        'raw_mass': raw_mass,
                        'ep': ep,
                        'video': filename
                        })
    df = pd.DataFrame(data)

    # Export to csv if user is happy with trajectories
    print("PRINTING CSV DATA")
    #continue_bool = input()
    #if(continue_bool == "y"):
    df.to_csv('../csvs/'+ csv_name +'.csv', mode='a', header = False) 

# Windows multiprocessing check
if __name__ == "__main__":
    
    csv_name = 'All_Particles_2'
    
    # Define the column names
    columns = [ 'row_num',
                'particle',
                'dx',
                'dy',
                'x',
                'y',
                'frame',
                'size',
                'ecc',
                'signal',
                'mass',
                'raw_mass',
                'ep',
                'video']

    # Create an empty DataFrame with the specified column names
    df = pd.DataFrame(columns=columns)

    # Write the DataFrame to a CSV file
    df.to_csv('../csvs/' + csv_name + '.csv', index=False)
    
    # Specify the folder containing the AVI video files
    folder_path = '../test_video/'

    # Loop through each AVI video file in the folder
    for filename in os.listdir(folder_path):
        if filename.endswith('.avi'):
            # Create a pims Video object for the current video file
            video = pims.Video(os.path.join(folder_path, filename))
            
            # Remove file extension from filename
            tempTuple = os.path.splitext(filename)
            filename = tempTuple[0]
            
            # Track particles
            print("PROCESSING VIDEO: " + filename)
            compute_traj(video, filename, csv_name)
