# -*- coding: utf-8 -*-
"""
Created on Thu Mar 23 09:12:12 2023

@author: alexm
"""

import numpy as np
import matplotlib.pyplot as plt
import random
from matplotlib.ticker import LogLocator

import Mueller_Matrices as MM


I = 1
Q = random.uniform(-1, 1)
U = random.uniform(-1, 1)
V = random.uniform(-1, 1)

p = np.sqrt(Q**2 + U **2 + V**2)
P = random.random()

Q/=p/P
U/=p/P
V/=p/P


S = MM.Stokes(I, Q, U, V).S

ang_arr = np.linspace(0,360,100)
rad_arr = np.deg2rad(ang_arr)

power_arr_0 = ang_arr * 0
power_arr_1 = ang_arr * 0 


M_HLP = MM.Mueller()
M_VLP = MM.Mueller()

M_HLP.linpol(0)
M_VLP.linpol(90)

MHLP = M_HLP.M
MVLP = M_VLP.M

M_qwp_ang = MM.Mueller()


for _a in range(len(ang_arr)): 
    a = ang_arr[_a]
    M_qwp_ang.qwp(a, 90)
    MQWP = M_qwp_ang.M
    power = np.matmul(MHLP, np.matmul(MQWP, S))
    power_arr_0[_a] = power[0]
    power = np.matmul(MVLP, np.matmul(MQWP, S))
    power_arr_1[_a] = power[0]


fig, ax = plt.subplots(subplot_kw={'projection': 'polar'}, figsize=(5,5), dpi = 300)

ax.plot(rad_arr, power_arr_0, color = "#DD7E68", label = r"$LP(0^o), QWP(\theta)$")
ax.plot(rad_arr, power_arr_1, color = "#68CFE0", label = r"$LP(90^o), QWP(\theta)$")

Title = np.asarray([[round(I, 3)],[round(Q, 3)],[round(U, 3)],[round(V, 3)]])

ax.set_title(f"{Title}", va='bottom', fontsize = 14)
plt.legend()
angle = np.deg2rad(67.5)
ax.legend(loc="lower left", bbox_to_anchor=(.5 + np.cos(angle)/2, .5 + np.sin(angle)/2), fontsize = 14)

plt.show()
