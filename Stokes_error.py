# -*- coding: utf-8 -*-
"""
Created on Fri Jul 14 06:11:12 2023

@author: alexm
"""

"""
Work in progress. This document calculates the contributions to stokes 
parameters from error. Consult Alex with any clarification needed. 
"""

import numpy as np
from scipy.stats import poisson, norm
from random import random
from random import seed 
import matplotlib.pyplot as plt
import Mueller_Matrices as MM

rot = 1
pp360 = 16
N = rot*pp360
rad_arr = np.linspace(0, rot*2*np.pi, N)
power_arr = np.zeros(N)

Int_avg = 500_000 * 20

I, Q, U, V = 1, 0, 0, 0
dop = Q**2+U**2+V**2
if dop >= 1: 
    Q /= dop
    U /= dop
    V /= dop
    
S = MM.Stokes(I, Q, U, V).S * Int_avg

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

power_arr_n = np.zeros(N)

sim_len = 10000

Aerr = np.zeros(sim_len)
Berr = np.zeros(sim_len)
Cerr = np.zeros(sim_len)
Derr = np.zeros(sim_len)

power_arr = power_arr * (np.max(rad_arr) - rad_arr)

for s in range(sim_len): 
    """
    seed(s)
    np.random.seed(seed=s)
    """
    
    for p in range(len(power_arr)): 
        r = poisson.rvs(power_arr[p], size=1)
        dns = 100 * 20
        det_noise = random() * 2*dns - dns
        #det_noise = random() * dns
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

print(f"""Average Error, STD (True Stokes Value): \n
I: {np.average(Ierr)}, {np.std(Ierr)} ({Irea})\n
Q: {np.average(Qerr)}, {np.std(Qerr)} ({Qrea})\n
U: {np.average(Uerr)}, {np.std(Uerr)} ({Urea})\n
V: {np.average(Verr)}, {np.std(Verr)} ({Vrea})\n""")

Xax = np.linspace(0, sim_len, sim_len, endpoint = False)
#%%

fig, ax = plt.subplots(dpi = 300)
size = 1
alpha = 0.5
ax.scatter(Xax, Ierr, c = "k", label = fr"$I_{{error}}$", s = size, alpha = alpha)
ax.scatter(Xax, Qerr, c = "r", label = fr"$Q_{{error}}$", s = size, alpha = alpha)
ax.scatter(Xax, Uerr, c = "g", label = fr"$U_{{error}}$", s = size, alpha = alpha)
ax.scatter(Xax, Verr, c = "b", label = fr"$V_{{error}}$", s = size, alpha = alpha)

ax.set_xlim(0, sim_len)
edge = ax.get_ylim()
em = max(-edge[0], edge[1])
ax.set_ylim(-em, em)
#ax.set_ylim(-0.02, 0.02)

ax.set_xlabel("Simulation number")
ax.set_ylabel("Change attributable to noise")
plt.legend(loc = "upper right")

plt.title(rf"$S = {I, Q, U, V},$  $I_{{avg}} = {Int_avg}$,  Detector noise scale = {dns}")

plt.show()


#%%
fig, (axI, axQ, axU, axV) = plt.subplots(1, 4, dpi = 300, figsize = (8, 3), sharex = True, sharey = True)
fig.tight_layout()
bnnumb = 100
bn = np.linspace(-0.001, 0.001, bnnumb)

axI.hist(Ierr, bins = bn, color = "k")
axQ.hist(Qerr, bins = bn, color = "r")
axU.hist(Uerr, bins = bn, color = "g")
axV.hist(Verr, bins = bn, color = "b")
lim = axI.get_xlim()
em = max(-lim[0], lim[1])
axI.set_xlim(-em, em)


axI.set_title("I error \n" + rf"$\mu = {round(np.average(Ierr), 6)}$" + "\n" + rf"$\sigma = {round(np.std(Ierr), 6)}$")
axQ.set_title("Q error \n" + rf"$\mu = {round(np.average(Qerr), 6)}$" + "\n" + rf"$\sigma = {round(np.std(Qerr), 6)}$")
axU.set_title("U error \n" + rf"$\mu = {round(np.average(Uerr), 6)}$" + "\n" + rf"$\sigma = {round(np.std(Uerr), 6)}$")
axV.set_title("V error \n" + rf"$\mu = {round(np.average(Verr), 6)}$" + "\n" + rf"$\sigma = {round(np.std(Verr), 6)}$")

plt.show()

