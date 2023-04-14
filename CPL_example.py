# -*- coding: utf-8 -*-
"""
Created on Fri Apr 14 14:28:31 2023

@author: alexm
"""

import CPL

maxl = 783

lower = maxl - 50
upper = maxl + 50

File_dict = {"Achiral 1": ["E:\Alex\CPL Data\CPL_3_07.csv", 20], 
"Achiral 2": ["E:\Alex\CPL Data\CPL_3_17.csv", 20]}

for n,f in File_dict.items(): 
    a = CPL.CPL_Analysis(f[0], f[1], name = n)
    a.plot_all_spectrum(Inset = True, inset_xlim=(510, 513), inset_ylim = (33750, 34500), inset_yb = 0.003, inset_xb = 0.001)
    a.plot_spectrum(0, vlines=(a.wavelengths[lower],a.wavelengths[upper]))
    a.plotIntegratedIQUV(Bounds = (lower, upper), ylim = (0.96, 1))
    a.plotCPL(Intensity_weighted="Both", units="mdeg", ylim = (-100,100), xlim=(450, 600), legend=True)
    a.plotglum(ylim = (-1, 2.5), xlim=(450, 600), Intensity_weighted = "Both", color = "k", scaling = 2)
    a.plotSpectrumIQUV(xlim = (400, 800), Inset = True, inset_ylim=(-0.01, 0.01), inset_xlim=(450, 600), inset_y_size = 0.2, inset_x_size = -0.05)
