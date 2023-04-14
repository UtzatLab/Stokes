# -*- coding: utf-8 -*-
"""
Created on Tue Apr 11 17:33:11 2023

@author: alexm
"""

import Simple_Stokes as SS
import matplotlib.pyplot as plt
import numpy as np

path = r"C:\Users\alexm\Desktop\Power Meter Logs"

filenames = []
for i in range(1, 9):
    filenames.append(f"DATA0{i}.CSV")
    
SS.PM100D_stokes(path, filenames, 22.5, 180, flip = -1)

#%%

spacing = 22.5

power_arr = np.array([1.103, 1.263, 1.155, 1.223, 1.161, .557, .0827, .438, 1.089, 1.267, 1.166, 1.245, 1.165, .5661, .1098, .4673, 1.104])
power_arr = power_arr / np.max(power_arr)

ang_arr = np.asarray([i*spacing for i in range(len(power_arr))])

SS.simple_stokes(ang_arr, power_arr, spacing, 360)

