# -*- coding: utf-8 -*-
"""
Created on Thu Aug 17 09:32:34 2023

@author: alexm

Work in progess/example code. Takes measured intensities as a function of 
original polarization, QWP angle, and wavelength to generate the Mueller matrix 
of a setup for white light. Currently intended for use with the ocean optics flame. 
Consult Alex with any clarifiaction needed. 
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import Simple_Stokes as SS
from pandas import read_csv
import os
import imageio

degcheck = np.load(r"C:\Users\alexm\Desktop\Galvos\Montana_CPL\Optimization_params\20230816\White_Excitation_Path_16-17-37_QWP_pos.npy")
int_arr_arr = np.load(r"C:\Users\alexm\Desktop\Galvos\Montana_CPL\Optimization_params\20230816\White_Excitation_Path_16-17-37_full_spectra.npy")

int_arr_arr = np.reshape(int_arr_arr, (5,-1, 3648))

deg_a = np.linspace(0, 360, 16, endpoint = False)
deg_arr = np.tile(deg_a, 5)

mpl.rcParams['figure.dpi'] = 300

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

pol_list = ["Horizontal", 
            "Vertical", 
            "Diagonal", 
            "RCP", 
            "LCP"]

pol_dict = {"Horizontal": {"I": [], 
                           "Q": [], 
                           "U": [], 
                           "V": [], 
                           "dop": [], 
                           "R2": []}, 
            "Vertical": {"I": [], 
                         "Q": [], 
                         "U": [], 
                         "V": [], 
                         "dop": [], 
                         "R2": []}, 
            "Diagonal": {"I": [], 
                         "Q": [], 
                         "U": [], 
                         "V": [], 
                         "dop": [], 
                         "R2": []}, 
            "RCP": {"I": [], 
                    "Q": [], 
                    "U": [], 
                    "V": [], 
                    "dop": [], 
                    "R2": []}, 
            "LCP": {"I": [], 
                    "Q": [], 
                    "U": [], 
                    "V": [], 
                    "dop": [], 
                    "R2": []}}

nw = len(wavelengths)
xlim = [400, 700]

for i in range(len(pol_list)):   
    I_arr = np.zeros(nw)
    Q_arr = np.zeros(nw)
    U_arr = np.zeros(nw)
    V_arr = np.zeros(nw)
    R2_arr = np.zeros(nw)
    dop_arr = np.zeros(nw)
    
    int_arr_arr_t = np.transpose(int_arr_arr[i])

    for w in range(wvl_ind(xlim[0]), wvl_ind(xlim[1])): 
        fit = SS.uneven_stokes(deg_arr, int_arr_arr_t[w], dontplot=True, ret = np.pi/2)
        I_arr[w] = fit["I"]
        Q_arr[w] = fit["Q"]
        U_arr[w] = fit["U"]
        V_arr[w] = fit["V"]
        R2_arr[w] = fit["R2"]
        dop_arr[w] = fit["dop"]
        
    pol_dict[pol_list[i]]["I"] = I_arr
    pol_dict[pol_list[i]]["Q"] = Q_arr
    pol_dict[pol_list[i]]["U"] = U_arr
    pol_dict[pol_list[i]]["V"] = V_arr
    pol_dict[pol_list[i]]["dop"] = dop_arr
    pol_dict[pol_list[i]]["R2"] = R2_arr
    
    Imax = I_arr #np.max(I_arr) #
    
    title = f"White light, initial {pol_list[i]}"
    
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

Hintt = 300/250
Vintt = 400/250

Cintt = 0.95


MM = np.zeros((4,4,3648))

MM[0][0] = 0.5 * (pol_dict["Horizontal"]["I"]/Hintt + pol_dict["Vertical"]["I"]/Vintt)
MM[0][1] = 0.5 * (pol_dict["Horizontal"]["I"]/Hintt - pol_dict["Vertical"]["I"]/Vintt)
MM[0][2] = pol_dict["Diagonal"]["I"] - MM[0][0]
MM[0][3] = pol_dict["RCP"]["I"]/Cintt - MM[0][0]

MM[1][0] = 0.5 * (pol_dict["Horizontal"]["Q"]/Hintt + pol_dict["Vertical"]["Q"]/Vintt)
MM[1][1] = 0.5 * (pol_dict["Horizontal"]["Q"]/Hintt - pol_dict["Vertical"]["Q"]/Vintt)
MM[1][2] = pol_dict["Diagonal"]["Q"] - MM[1][0]
MM[1][3] = pol_dict["RCP"]["Q"]/Cintt  - MM[1][0]

MM[2][0] = 0.5 * (pol_dict["Horizontal"]["U"]/Hintt + pol_dict["Vertical"]["U"]/Vintt)
MM[2][1] = 0.5 * (pol_dict["Horizontal"]["U"]/Hintt - pol_dict["Vertical"]["U"]/Vintt)
MM[2][2] = pol_dict["Diagonal"]["U"] - MM[2][0]
MM[2][3] = pol_dict["RCP"]["U"]/Cintt  - MM[2][0]

MM[3][0] = 0.5 * (pol_dict["Horizontal"]["V"]/Hintt + pol_dict["Vertical"]["V"]/Vintt)
MM[3][1] = 0.5 * (pol_dict["Horizontal"]["V"]/Hintt - pol_dict["Vertical"]["V"]/Vintt)
MM[3][2] = pol_dict["Diagonal"]["V"] - MM[3][0]
MM[3][3] = pol_dict["RCP"]["V"]/Cintt  - MM[3][0]



color = "b"
lw = 1

norm = MM[0][0]

fig, ((m00, m01, m02, m03), (m10, m11, m12, m13), (m20, m21, m22, m23), (m30, m31, m32, m33),) = plt.subplots(4, 4, dpi = 300, sharey = True, sharex = True, figsize=(6,4))

m00.plot(wavelengths, MM[0][0]/norm, c = color, lw = lw)
m01.plot(wavelengths, MM[0][1]/norm, c = color, lw = lw)
m02.plot(wavelengths, MM[0][2]/norm, c = color, lw = lw)
m03.plot(wavelengths, MM[0][3]/norm, c = color, lw = lw)

m10.plot(wavelengths, MM[1][0]/norm, c = color, lw = lw)
m11.plot(wavelengths, MM[1][1]/norm, c = color, lw = lw)
m12.plot(wavelengths, MM[1][2]/norm, c = color, lw = lw)
m13.plot(wavelengths, MM[1][3]/norm, c = color, lw = lw)

m20.plot(wavelengths, MM[2][0]/norm, c = color, lw = lw)
m21.plot(wavelengths, MM[2][1]/norm, c = color, lw = lw)
m22.plot(wavelengths, MM[2][2]/norm, c = color, lw = lw)
m23.plot(wavelengths, MM[2][3]/norm, c = color, lw = lw)

m30.plot(wavelengths, MM[3][0]/norm, c = color, lw = lw)
m31.plot(wavelengths, MM[3][1]/norm, c = color, lw = lw)
m32.plot(wavelengths, MM[3][2]/norm, c = color, lw = lw)
m33.plot(wavelengths, MM[3][3]/norm, c = color, lw = lw)

fig.tight_layout()
m00.set_ylim(-1,1)

plt.xlim(410, 700)
plt.xlabel("Wavelength")
fig.suptitle("Excitation Pathway", fontsize = 18, y = 1.05)
plt.show()



title = f"White light, LCP"

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
ax.plot(wavelengths, Q_arr/Imax, c="r", label = "Q")
ax.plot(wavelengths, U_arr/Imax, c="g", label = "U")
ax.plot(wavelengths, V_arr/Imax, c="b", label = "V")
ax.plot(wavelengths, wavelengths*0, c="k", ls = "dashed")

I_arrRe = MM[0][0]/norm - MM[0][3]/norm
Q_arrRe = MM[1][0]/norm - MM[1][3]/norm
U_arrRe = MM[2][0]/norm - MM[2][3]/norm
V_arrRe = MM[3][0]/norm - MM[3][3]/norm

ImaxRe = I_arrRe
alpha = 0.5

ax.plot(wavelengths, I_arrRe/np.max(I_arrRe) , c="gray", label = "Int (a.u.)", ls = (0, (3, 1, 1, 1)), alpha = alpha)
ax.plot(wavelengths, Q_arrRe/ImaxRe, c="r", label = "Q", alpha = alpha)
ax.plot(wavelengths, U_arrRe/ImaxRe, c="g", label = "U", alpha = alpha)
ax.plot(wavelengths, V_arrRe/ImaxRe, c="b", label = "V", alpha = alpha)
ax.plot(wavelengths, wavelengths*0, c="k", ls = "dashed")


plt.title(title)

h, l = plt.gca().get_legend_handles_labels()
order = [0,2,3,4,1,5]


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

