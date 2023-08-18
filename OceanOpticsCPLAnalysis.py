# -*- coding: utf-8 -*-
"""
Created on Fri May 19 16:57:26 2023

@author: alexm

Working code, lacks use notes. Takes spectrometer output as a function of QWP 
angle to generate Stokes parameters as a function of wavelength. Includes other 
functions I used on occasion but are generally not used.
Consult Alex with any clarification needed.

"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import Simple_Stokes as SS
from pandas import read_csv
import os
import imageio

degcheck = np.load(r"C:\Users\alexm\Desktop\Galvos\Montana_CPL\Optimization_params\20230816\White_Excitation_Path_16-17-37_QWP_pos.npy")
check = np.load(r"C:\Users\alexm\Desktop\Galvos\Montana_CPL\Optimization_params\20230816\White_Excitation_Path_16-17-37_full_spectra.npy")

check = check[:5*16*3648*5]
check = np.reshape(check, (5,-1, 3648))


int_arr_arr = check
deg_a = np.linspace(0, 360, 16, endpoint = False)
deg_arr = np.tile(deg_a, 5)

mpl.rcParams['figure.dpi'] = 300


#deg_arr = np.load(r"C:\Users\alexm\Desktop\Galvos\CPL Data\2023\08\10\CdSe-CdS_glassvial_ci_Right excitation_10rot_int11000ms_17-00-06_d.npy")
#int_arr_arr = np.load(r"C:\Users\alexm\Desktop\Galvos\CPL Data\2023\08\10\CdSe-CdS_glassvial_ci_Right excitation_10rot_int11000ms_17-00-06_i.npy")


n = len(int_arr_arr[0])
# Format: int_arr_arr[scan number][wavelength number]

wavelengths = np.load(r"C:\Users\alexm\Desktop\Galvos\CPL Data\wavelengths.npy")
def wvl_ind(goal): 
    for i in range(3648): 
        if wavelengths[i] <= goal and wavelengths[(i+1)%3648] > goal: 
            return i
        
pol = 4

spacing = 1
cmap = mpl.colormaps['viridis'] # Continuous 
#cmap = mpl.colormaps['hsv'] # Cyclic

fig, ax = plt.subplots()
for i in range(n):
    if i%spacing == 0: 
        ax.plot(wavelengths, int_arr_arr[pol][i], color = cmap(i/n), alpha = 0.5) 
        plt.hlines(0, 0, 700000, color = "k")
        # plt.vlines(498.5, 0, np.max(int_arr_arr[i]))
        # color = cmap(deg_mod_norm[i]) to plot by angle
        # color = cmap(i/n) to plot by time 
        
plt.xlim(400, 700)
plt.show()
        
#%%
deg_arr = deg_arr[:n]
int_arr_arr = int_arr_arr[:n, :]
tint_arr_arr = int_arr_arr[:n, wvl_ind(505):wvl_ind(515)]
integrated_arr = np.zeros(n)

for i in range(n): 
    integrated_arr[i] = np.sum(tint_arr_arr[i])
    
imax = np.max(integrated_arr)
imin = np.min(integrated_arr)
integrated_arr = integrated_arr #/ imax
    
SS.uneven_stokes(deg_arr, integrated_arr, alpha = 1, fit_top=False)
plt.ylim(0, 1*imax)
plt.ylim(imin, imax)


#%%

nw = len(wavelengths)

I_arr = np.zeros(nw)
Q_arr = np.zeros(nw)
U_arr = np.zeros(nw)
V_arr = np.zeros(nw)
R2_arr = np.zeros(nw)
dop_arr = np.zeros(nw)

int_arr_arr_t = np.transpose(int_arr_arr[pol])
# Format: int_arr_arr[wavelength number][scan number]

xlim = [400, 700]

for w in range(wvl_ind(xlim[0]), wvl_ind(xlim[1])): 
    fit = SS.uneven_stokes(deg_arr, int_arr_arr_t[w], dontplot=True, ret = np.pi/2)
    I_arr[w] = fit["I"]
    Q_arr[w] = fit["Q"]
    U_arr[w] = fit["U"]
    V_arr[w] = fit["V"]
    R2_arr[w] = fit["R2"]
    dop_arr[w] = fit["dop"]

#%%
Imax = I_arr #np.max(I_arr) #

title = r"White light, initial LCP"

fig, ax = plt.subplots()
ax.set_xlim(xlim[0], xlim[1])
ax.set_xlim(420, 700)
ax.set_xlabel("Wavelength (nm)")
#ax.set_ylim(-0.5, 1)
#ax.set_ylim(-0.05, 0.05)
#ax.set_ylim(-0.2,0.2)
ax.set_ylim(-1.1, 1.1)

ylim = ax.get_ylim()
ax.plot(wavelengths, I_arr/np.max(I_arr) , c="gray", label = "Int (a.u.)", ls = (0, (3, 1, 1, 1)))
ax.plot(wavelengths, dop_arr, c="y", label = "dop")
#ax.plot(wavelengths, I_arr/Imax, c="gray", label = "I")
ax.plot(wavelengths, Q_arr/Imax, c="r", label = "Q")
ax.plot(wavelengths, U_arr/Imax, c="g", label = "U")
ax.plot(wavelengths, V_arr/Imax, c="b", label = "V")
ax.plot(wavelengths, R2_arr, c="k", label = "R2")
ax.plot(wavelengths, wavelengths*0, c="k", ls = "dashed")

plt.title(title)

h, l = plt.gca().get_legend_handles_labels()
order = [0,2,3,4,1,5]

plt.legend([h[i] for i in order], [l[i] for i in order], loc= "lower right")

plt.show()

#%%

Imax = I_arr #np.max(I_arr) #

title = r"CdSe-CdS in glass pipette, HLP excitation"

fig, ax = plt.subplots()
ax.set_xlim(xlim[0], xlim[1])
ax.set_xlabel("Wavelength (nm)")
#ax.set_ylim(-0.5, 1)
ax.set_ylim(-0.5, 0.5)
#ax.set_ylim(-0.2,0.2)
ax.set_ylim(-1, 1)

size = 1

ylim = ax.get_ylim()
ax.plot(wavelengths, I_arr/np.max(I_arr) * ylim[1], c="gray", label = "Int (a.u.)", ls = (0, (3, 1, 1, 1)))
ax.scatter(wavelengths, dop_arr, c="y", label = "dop", alpha = R2_arr, edgecolor = "none", s=size)
#ax.plot(wavelengths, I_arr/Imax, c="gray", label = "I")
ax.scatter(wavelengths, Q_arr/Imax, c="r", label = "Q", alpha = R2_arr, edgecolor = "none", s=size)
ax.scatter(wavelengths, U_arr/Imax, c="g", label = "U", alpha = R2_arr, edgecolor = "none", s=size)
ax.scatter(wavelengths, V_arr/Imax, c="b", label = "V", alpha = R2_arr, edgecolor = "none", s=size)
#ax.plot(wavelengths, R2_arr, c="k", label = "R2")
ax.plot(wavelengths, wavelengths*0, c="k", ls = "dashed")

plt.title(title)

h, l = plt.gca().get_legend_handles_labels()
order = [0,2,3,4,1,5]

plt.legend([h[i] for i in order], [l[i] for i in order], loc= "upper right")

plt.show()

#%%

def extrema(array): 
    min_ = np.min(array[wvl_ind(xlim[0]):wvl_ind(xlim[1])])
    max_ = np.max(array[wvl_ind(xlim[0]):wvl_ind(xlim[1])])
    if max_ > -min_: 
        return max_
    else:
        return min_

Im = extrema(I_arr)
Qm = extrema(Q_arr)
Um = extrema(U_arr)
Vm = extrema(V_arr)

print(f"{Im}, {Qm}, {Um}, {Vm}")

#%%

average = np.average(V_arr[wvl_ind(505): wvl_ind(525)]/I_arr[wvl_ind(505): wvl_ind(525)])
print(average)
print(np.arcsin(average) * 28.648 * 1000)

#%%

print(np.min(Q_arr))

#%%

nw = len(wavelengths)

directory = r"C:\Users\alexm\Desktop\Galvos\CPL Data\2023\07\19"
filenames = []

for naa in range(58):
    
    int_arr_arr2 = np.reshape(int_arr_arr[16*(naa):16*(naa+1),:], (16, -1, 3648))
    int_arr_arr2 = np.sum(int_arr_arr2, axis = 1)
    
    I_arr = np.zeros(nw)
    Q_arr = np.zeros(nw)
    U_arr = np.zeros(nw)
    V_arr = np.zeros(nw)
    R2_arr = np.zeros(nw)
    dop_arr = np.zeros(nw)
    
    int_arr_arr_t = np.transpose(int_arr_arr2)
    # Format: int_arr_arr[wavelength number][scan number]
    
    for w in range(wvl_ind(475), wvl_ind(580)): 
        fit = SS.uneven_stokes(deg_arr, int_arr_arr_t[w], dontplot=True, ret = np.pi/2)
        I_arr[w] = fit["I"]
        Q_arr[w] = fit["Q"]
        U_arr[w] = fit["U"]
        V_arr[w] = fit["V"]
        R2_arr[w] = fit["R2"]
        dop_arr[w] = fit["dop"]
    
    
    Imax = I_arr #np.max(I_arr) #
    
    fig, ax = plt.subplots()
    ax.plot(wavelengths, I_arr/Imax, c="gray", label = "I")
    ax.plot(wavelengths, Q_arr/Imax, c="r", label = "Q")
    ax.plot(wavelengths, U_arr/Imax, c="g", label = "U")
    ax.plot(wavelengths, V_arr/Imax, c="b", label = "V")
    ax.plot(wavelengths, R2_arr, c="k", label = "R2")
    ax.plot(wavelengths, dop_arr, c="y", label = "dop", alpha = 0.3)
    ax.plot(wavelengths, wavelengths*0, c="k", ls = "dashed")
    
    #ax.plot(wavelengths, spectrum/100000, c = "m", label = "laser")
    #plt.vlines([508.5, 510.5], -1, 1)
    #plt.hlines(0, 0, 10000)
    
    plt.xlim(480, 580)
    plt.xlabel("Wavelength (nm)")
    
    plt.ylim(-0.1, 1)
    plt.ylim(-0.02, 0.02)
    plt.title(f"{naa}")
    plt.legend(loc= "upper right")
    #plt.show()

    plt.savefig(os.path.join(directory, f"{naa}.png"))
    filenames.append(os.path.join(directory, f"{naa}.png"))
    
images = []
for file_name in filenames:
    images.append(imageio.imread(file_name))

imageio.mimsave(os.path.join(directory, "RR_by_rotation.gif"), images)

for file_name in filenames: 
    os.remove(file_name)
    
#%%

fig, ax = plt.subplots()
ax.plot(wavelengths, np.abs(np.rad2deg(0.5 * np.arctan2(U_arr, Q_arr))), c="k", label = "angle")
plt.xlim(600, 700)
#plt.xlabel("Wavelength (nm)")
plt.ylabel(r"$\Psi$ (deg)")
plt.ylim(0, 180)
plt.show()

#%%

fig, ax = plt.subplots(subplot_kw={'projection': 'polar'}, dpi = 300)
ax.set_rorigin(xlim[0] - 10)
ax.set_xlim(0, np.pi)
ax.set_ylim(xlim[0], xlim[1])
color_arr = np.linspace(0, 1, len(wavelengths))
for c in range(len(color_arr)): 
    if wavelengths[c] > ax.get_ylim()[1] or wavelengths[c] < ax.get_ylim()[0]: 
        color_arr[c] = np.nan
m = np.nanmin(color_arr)
color_arr = color_arr - m
color_arr = np.nan_to_num(color_arr)
ax.scatter(((0.5 * np.arctan2(U_arr, Q_arr)))%np.pi, wavelengths, alpha = 1, s=1.25, cmap = "jet", c=color_arr, edgecolors='none', label = "angle")
#ax.plot(np.linspace(0, np.pi * 2, 1000), np.zeros(1000) + 510, label = "Emission Max", c = "g", alpha = 0.5, ls = "dashed")
#ax.plot(np.linspace(0, np.pi * 2, 1000), np.zeros(1000) + 484, label = "Laser", c = "c")


ax.set_title(title + "\n" + r"$\Psi$ (deg)", pad = -30)
ax.set_rticks(np.arange(xlim[0], xlim[1]+1, 100))
ax.tick_params(axis="both", pad = 2)
ax.set_xlabel("Wavelength (nm)", labelpad = -30, loc = "right")

plt.show()

aaa = np.abs(np.rad2deg(0.5 * np.arctan2(U_arr, Q_arr)))
#61810

    
#%%

M = np.array([[ 0.37069409, -0.11630367,  0.0618421,  -0.03039849],
              [-0.15832925,  0.17280737, -0.09488918,  0.04164628],
              [ 0.00446851,  0.03991156,  0.05940509,  0.01040715],
              [-0.00953322, -0.0058859,  -0.00271213,  0.05697619]])

'''M = np.array([[ 0.59209853, -0.12056972, 0.05405228, -0.02867919-8.68707335e-11j], 
 [-0.15595252, 0.40365554, -0.12993679, 0.06065545-5.87834469e-11j],
 [ 0.01714056, 0.06362684, 0.25854014, 0.0141345],
 [-0.013382, -0.01098193, -0.00685854, 0.23949]])'''

Minv = np.linalg.inv(M)

S_arr = np.array([[I_arr[i], Q_arr[i], U_arr[i], V_arr[i]]  for i in range(len(wavelengths))])

Strue = S_arr * 0

for i in range(len(wavelengths)): 
    output = np.matmul(Minv, S_arr[i])
    for j in range(4): 
        Strue[i][j] = output[j]

S_arr_t = np.transpose(Strue)

Imax = S_arr_t[0] #np.max(I_arr) #

title = r"RR-Perovksites, RCP Excitation"

fig, ax = plt.subplots()
ax.set_xlim(xlim[0], xlim[1])
ax.set_xlabel("Wavelength (nm)")

#ax.set_ylim(-0.002, 0.002)

ylim = ax.get_ylim()
ax.plot(wavelengths, S_arr_t[0]/Imax, c="gray", label = "I")
ax.plot(wavelengths, S_arr_t[1]/Imax, c="r", label = "Q")
ax.plot(wavelengths, S_arr_t[2]/Imax, c="g", label = "U")
ax.plot(wavelengths, S_arr_t[3]/Imax, c="b", label = "V")
ax.plot(wavelengths, wavelengths*0, c="k", ls = "dashed")
ax.plot(wavelengths, (S_arr_t[1]/Imax)**2 + (S_arr_t[2]/Imax)**2 + (S_arr_t[3]/Imax)**2, c = "y")

#ax.set_ylim(-0.005,0.005)
ax.set_ylim(-1, 2)
ax.set_xlim(480, 600)

plt.legend()

plt.show()


M = np.array([[  1,   -0.31,  0.17, -0.08],
              [-0.43,  0.47, -0.26,  0.11],
              [ 0.01,  0.11,  0.16,  0.03],
              [-0.03, -0.02, -0.01,  0.15]])

