#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 04/25/2023 for The Ohio State University, Department of Physics
Calculates radius and charges from tracked velocity data in csv.
Please refer to Millikan_Experimental_Report.pdf
Credit to ricktjwong
@author: Hyun Wallace Anderson
"""

import math
import pandas as pd
import os
# print(os.getcwd())
# Datasheet to be transformed:
filename = "All_Particles"

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
turning_treshold = 10 # +/- frames around turning point
frame_rate = 30.893 # (fps)
distance_scale = 0.4259259259E-03 # Microscope increments in m
pixel_scale = 4.47047253E-06 # (m/px) using glass reference
pixel_scale_microscope = 4.61508209E-6 # (m/px) using microscope units

# Errors
frame_rate_err = 0.4465 # (fps)
pixel_scale_err = pixel_scale # 2.70502868E-8 (m/px) taking into account increment etching width
plate_err = 4.41E-7 # meters
voltage_err = 0.05 # Volts
    
df = pd.read_csv('Documents/Course Work/Physics 5700/Oil Drop Lab/millikan-master/csvs/' + filename + '.csv')

def get_ep(i):
    grouped = df.groupby('particle')
    return float(grouped['ep'].mean().iloc[0])

# Only consider particles which change direction within X frames from the midpoint of the video
def returnParticles(middle):
    filtered = []
    turning_indices = []
    for i in pivoted_y:
        turning_index = pivoted_y[i].idxmin()
        if abs(turning_index - middle) < turning_treshold:
            filtered.append(i)
            turning_indices.append(turning_index)
    return filtered, turning_indices

def filterParticles():
    turning_indices = []
    for i in pivoted_y:
        turning_index = pivoted_y[i].idxmin()
        turning_indices.append(turning_index)
    middle = int(sum(turning_indices)/len(turning_indices))
    return middle

def returnFinal():
    final = []
    for i, j in zip(filtered_particles, turning_indices):
        particle = pivoted_dy[i][pd.Series.notnull(pivoted_dy[i])]
        particle_ep = pivoted_ep[i][pd.Series.notnull(pivoted_ep[i])]
        bef = particle[:j-particle.index[0]]
        aft = particle[j-particle.index[0]:]
        # ep = abs(particle_ep.mean())
        ep = abs(particle_ep.std())
        if len(bef) > 0 and len(aft) > 0 and ep < 0.25:
            v_g = bef.sum()/len(bef)
            v_e = aft.sum()/len(aft) # TO-DO: if dx is less than dy
            if v_g < 0 and v_e > 0 and abs(v_e) > abs(v_g): # ensure v_g is negative, v_e is positive, and v_g is smaller in magnitude compared to v_e):
                ve_fin = abs(v_e*frame_rate*pixel_scale)
                vg_fin = abs(v_g*frame_rate*pixel_scale)
                radius = math.sqrt((9 * N * abs(vg_fin)) / (2 * G * DENSITY_OIL))
                charge = ((6 * PI * PLATE_DISTANCE) / voltage)*math.sqrt((9 * pow(N, 3)) / (2 * DENSITY_OIL * G))*pow(1 + (B / (p * radius)), -3/2)*((ve_fin + vg_fin) * math.sqrt(abs(vg_fin)))
                quantization = abs(charge / CHARGE_REAL)
                
                ve_err = math.sqrt(pow(v_e*pixel_scale*frame_rate_err, 2) + pow(v_e*pixel_scale_err*ep*frame_rate, 2))
                vg_err = math.sqrt(pow(v_g*pixel_scale*frame_rate_err, 2) + pow(v_g*pixel_scale_err*ep*frame_rate, 2))
                plate_d_prop = ((6 * PI * plate_err) / voltage)*math.sqrt((9 * pow(N, 3)) / (2 * DENSITY_OIL * G)) * pow(1 + (B / (p * radius)), -3/2)*((ve_fin + vg_fin) * math.sqrt(abs(vg_fin)))
                voltage_prop = voltage_err*((6 * PI * PLATE_DISTANCE) / (pow(voltage, 2)))*math.sqrt((9 * pow(N, 3)) / (2 * DENSITY_OIL * G)) * pow(1 + (B / (p * radius)), -3/2)*((ve_fin + vg_fin) * math.sqrt(abs(vg_fin)))
                ve_prop = (27 * math.sqrt(3/2) * PI * PLATE_DISTANCE * pow(abs(vg_fin),3/4) * (3 * ve_fin + vg_fin) * math.sqrt(pow(N,3) / (G * DENSITY_OIL))) / (math.sqrt(ve_fin) * voltage * math.sqrt(3 * math.sqrt(abs(vg_fin)) + (math.sqrt(2) * B * math.sqrt(G * DENSITY_OIL / N) / p) * ((3 * p * math.sqrt(abs(vg_fin)) + (math.sqrt(2) * B * math.sqrt(G * DENSITY_OIL / N)))))) * ve_err
                vg_prop = (27 * math.sqrt(3) * PI * PLATE_DISTANCE * p * math.sqrt(ve_fin) * math.sqrt(pow(N,3)/(G * DENSITY_OIL))) * ((3 * B * (vg_fin + ve_fin) * math.sqrt(G * DENSITY_OIL / N)) + (4 * B * abs(vg_fin) * math.sqrt(G * DENSITY_OIL / N)) + (6 * math.sqrt(2) * p * pow(abs(vg_fin), 3/2))) / ((2 * voltage * pow(abs(vg_fin), 1/4) * math.sqrt((3 * math.sqrt(abs(vg_fin))) + (math.sqrt(2) * B * math.sqrt(G * DENSITY_OIL / N) / p))) * pow(((3 * p * math.sqrt(abs(vg_fin))) + (math.sqrt(2) * B * math.sqrt(G * DENSITY_OIL / N))) ,2)) * vg_err
                charge_err = math.sqrt(pow(plate_d_prop, 2) + pow(voltage_prop, 2) + pow(vg_prop,2))

                final.append({'particle': i, 
                            'turning_point': j, 
                            'v_e (px/frame)': v_e, 
                            'v_g (px/frame)': v_g,
                            'v_e (m/s)': ve_fin, 
                            'v_g (m/s)': vg_fin,
                            'radius': radius,  
                            'charge': charge,
                            'plate err': plate_d_prop,
                            'voltage err': voltage_prop,
                            'v_e prop': ve_prop,
                            'v_g prop': vg_prop,
                            'v_e err': ve_err,
                            'v_g err': vg_err,
                            'ep': ep,
                            'charge err': charge_err,
                            'quantization': quantization,
                            })        
    return final

pivoted = df.pivot(index = 'frame', columns = 'particle')
pivoted_y = pivoted['y']
pivoted_dy = pivoted['dy']
pivoted_ep = pivoted['ep']

filtered_particles, turning_indices = returnParticles(250)
print (len(filtered_particles))
final = returnFinal()
df = pd.DataFrame(final)
df.to_csv('Documents/Course Work/Physics 5700/Oil Drop Lab/millikan-master/python/'+ filename +'_transformed.csv')