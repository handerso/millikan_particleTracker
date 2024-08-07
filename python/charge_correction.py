#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 04/25/2023 for The Ohio State University, Department of Physics
Please refer to Millikan_Experimental_Report.pdf
Credit to ricktjwong
@author: Hyun Wallace Anderson
"""

import math
import csv
import pandas as pd
import matplotlib.pyplot as plt
import trackpy as tp

# Datasheet to be transformed:
filename = "Particle_Data"
        
# Constants
PI = math.pi
G = 9.8067 # (m/s^2)
PLATE_DISTANCE = 3.8965E-3 # (m)
DENSITY_OIL = 860 # (kg/m^3)
B = 6.17E-6 # correction factor
N = 1.825E-5 # Viscosity of air at 20 deg_C (N*s/m^2)
CHARGE_REAL = 1.60217663E-19 # coulombs

# Independant variables
voltage = 450.0 # Voltage (V)
p = 73.818545 # Barometric Pressure (cm/Hg)

# Calibration constants
frame_rate = 30.893 # (fps)
distance_scale = 0.4259259259E-03 # Microscope increments in m
pixel_scale = 4.47047253E-06 # (m/px)
    
# df = pd.read_csv('../csvs/' + filename + '.csv')

# open the CSV file for reading and writing
with open('../' + filename + '.xlsm', 'r+') as csvfile:
    reader = csv.reader(csvfile)
    writer = csv.writer(csvfile)

    # iterate over each row in the CSV file
    for row in reader:
        # update the value in the second column of each row
        ve_fin = row[5]
        vg_fin = row[6]
        row[7] = math.sqrt((9 * N * vg_fin) / (2 * G * DENSITY_OIL))
        radius = row[7]
        row[8] = ((6 * PI * PLATE_DISTANCE) / voltage)*math.sqrt((9 * pow(N, 3)) / (2 * DENSITY_OIL * G))*pow(1 + (B / (p * radius)), -3/2)*((ve_fin + vg_fin) * math.sqrt(vg_fin))
        row[9] = row[8] / CHARGE_REAL

        # write the updated row back to the CSV file
        writer.writerow(row)

# df.to_csv('../' + filename + '_corrected.csv')