# -*- coding: utf-8 -*-
"""
Created on Fri Jul 14 08:16:12 2023

@author: alexm
"""

"""
Work in progress. This document calculates the correlation 
in error between the different stokes parameters. Consult Alex with
any clarification needed. 
"""

import numpy as np
from scipy.stats import poisson
from random import random
from random import seed 
import matplotlib.pyplot as plt
import Mueller_Matrices as MM

rot = 1
pp360 = 16
N = rot*pp360
sim_len = 100
offset = 345

Int_len = 50
Int_arr = np.zeros(Int_len)


I, Q, U, V = 1, 0, 0, -0.001
dop = Q**2+U**2+V**2
if dop >= 1: 
    Q /= dop
    U /= dop
    V /= dop

II_errcorr = np.zeros(Int_len)
QQ_errcorr = np.zeros(Int_len)
UU_errcorr = np.zeros(Int_len)
VV_errcorr = np.zeros(Int_len)

IQ_errcorr = np.zeros(Int_len)
IU_errcorr = np.zeros(Int_len)
IV_errcorr = np.zeros(Int_len)
QU_errcorr = np.zeros(Int_len)
QV_errcorr = np.zeros(Int_len)
UV_errcorr = np.zeros(Int_len)



for I in range(Int_len):
    Int = 10 ** (I/10) + 100    
    Int_arr[I] = Int
    
    rad_arr = np.linspace(0, rot*2*np.pi, N)
    power_arr = np.zeros(N)
    power_arr_n = np.zeros(N)
    
    Is, Q, U, V = 1, 0, 0, -0.001
    dop = Q**2+U**2+V**2
    if dop >= 1: 
        Q /= dop
        U /= dop
        V /= dop
    
    S = MM.Stokes(Is, Q, U, V).S * Int

    M_HLP = MM.Mueller()
    M_HLP.linpol(0)
    MHLP = M_HLP.M
    M_qwp_ang = MM.Mueller()

    for _a in range(N): 
        a = np.rad2deg(rad_arr[_a])
        M_qwp_ang.qwp(a, 90)
        MQWP = M_qwp_ang.M
        power = np.matmul(MHLP, np.matmul(MQWP, S))
        power_arr[_a] = power[0]
       
    Aerr = np.zeros(sim_len)
    Berr = np.zeros(sim_len)
    Cerr = np.zeros(sim_len)
    Derr = np.zeros(sim_len)
    
    for s in range(sim_len): 
        seed(s+offset)
        np.random.seed(seed=s+offset)
        
        for p in range(len(power_arr)): 
            r = poisson.rvs(power_arr[p], size=1)
            dns = 100
            det_noise = random() * 2*dns - dns
            det_noise = random() * dns
            power_arr_n[p] = r[0] + det_noise
        
        diff_arr = power_arr_n - power_arr
        
        M = np.zeros((4,4))
        Mlong = np.zeros((4, N))
        Slong = np.zeros(N)
        
        for _a in range(N):
            theta = rad_arr[_a]
            s2 = np.sin(2*theta)
            c4 = np.cos(4*theta)
            s4 = np.sin(4*theta)
        
            M[0][0] += 1
            M[0][1] += s2
            M[0][2] += c4
            M[0][3] += s4
            
            M[1][0] += s2
            M[1][1] += s2 * s2 
            M[1][2] += c4 * s2
            M[1][3] += s4 * s2
            
            M[2][0] += c4
            M[2][1] += s2 * c4
            M[2][2] += c4 * c4
            M[2][3] += s4 * c4
            
            M[3][0] += s4
            M[3][1] += s2 * s4
            M[3][2] += c4 * s4
            M[3][3] += s4 * s4
            
            Mlong[0][_a] = 1
            Mlong[1][_a] = s2
            Mlong[2][_a] = c4
            Mlong[3][_a] = s4
            
            Slong[_a] = power_arr[_a]
            
            
        Minv = np.linalg.inv(M)
        
        Mspec = np.matmul(Minv, Mlong)
        
        S1 = 2 * np.matmul(Mspec, power_arr)
        S2 = 2 * np.matmul(Mspec, diff_arr)
        
        Aerr[s] = S2[0]
        Berr[s] = S2[1]
        Cerr[s] = S2[2]
        Derr[s] = S2[3]
    
    Area = S1[0]
    Brea = S1[1]
    Crea = S1[2]
    Drea = S1[3]
    
    Irea = Area - Crea
    norm = Irea 
    
    Irea = Irea / norm
    Qrea = 2 * Crea / norm
    Urea = 2 * Drea / norm
    Vrea = Brea / norm
    
    Ierr = (Aerr - Cerr) / norm
    
    Qerr = 2 * Cerr / norm
    Uerr = 2 * Derr / norm
    Verr = Berr / norm

    II_errcorr[I] = np.average(Ierr*Ierr) - np.average(Ierr)*np.average(Ierr)
    QQ_errcorr[I] = np.average(Qerr*Qerr) - np.average(Qerr)*np.average(Qerr)
    UU_errcorr[I] = np.average(Uerr*Uerr) - np.average(Uerr)*np.average(Uerr)
    VV_errcorr[I] = np.average(Verr*Verr) - np.average(Verr)*np.average(Verr)
    
    IQ_errcorr[I] = np.average(Ierr*Qerr) - np.average(Ierr)*np.average(Qerr)
    IU_errcorr[I] = np.average(Ierr*Uerr) - np.average(Ierr)*np.average(Uerr)
    IV_errcorr[I] = np.average(Ierr*Verr) - np.average(Ierr)*np.average(Verr)
    QU_errcorr[I] = np.average(Qerr*Uerr) - np.average(Qerr)*np.average(Uerr)
    QV_errcorr[I] = np.average(Qerr*Verr) - np.average(Qerr)*np.average(Verr)
    UV_errcorr[I] = np.average(Uerr*Verr) - np.average(Uerr)*np.average(Verr)
#%%

Int_arrl = Int_arr # np.log10(Int_arr)

fig, ax = plt.subplots(dpi = 300) 
plt.xscale("log")
ax.plot(Int_arrl, II_errcorr, c = "k", label = "I:I")
ax.plot(Int_arrl, QQ_errcorr, c = "r", label = "Q:Q")
ax.plot(Int_arrl, UU_errcorr, c = "g", label = "U:U")
ax.plot(Int_arrl, VV_errcorr, c = "b", label = "V:V")
ax.plot(Int_arrl, np.linspace(0, 0.0001, len(Int_arrl)), alpha = 0, label = " ")
ax.plot(Int_arrl, IQ_errcorr, c = "#b3280a", label = "I:Q")
ax.plot(Int_arrl, IU_errcorr, c = "#005500", label = "I:U")
ax.plot(Int_arrl, IV_errcorr, c = "#000044", label = "I:V")
ax.plot(Int_arrl, QU_errcorr, c = "#efb102", label = "Q:U")
ax.plot(Int_arrl, QV_errcorr, c = "#9933ff", label = "Q:V")
ax.plot(Int_arrl, UV_errcorr, c = "#009999", label = "U:V")
ax.set_xlim(np.min(Int_arrl), np.max(Int_arrl))
ax.set_xlabel("Average Intensity")
ax.set_ylabel("Correlation of errors (<XY> -<X><Y>)")
plt.legend(bbox_to_anchor=(1.2, 0.9), loc="upper right")
ax.set_ylim(-2,2)
plt.title(f"Detector noise scale: {dns} counts. Initial seed = {offset}")
plt.show()


