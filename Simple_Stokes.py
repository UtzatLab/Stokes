# -*- coding: utf-8 -*-
"""
Created on Tue Apr 11 15:48:19 2023

@author: alexm

Simple functions to quickly calculate Stokes parameters from various data types. 
For more complex opperations (Parameters by wavelength, CD/CPL, etc.) use the 
CPL library.

Recommended style: 
    import Simple_Stokes as SS
    
    SS.simple_stokes(~parameters~)

Currently supported data types: 
    simple_stokes: simple arrays of intensities with matching QWP angle arrays.
    PM100D_stokes: .CSV files generated from a Thorlabs PM100D power meter.
    
"""

import os
import numpy as np
from pandas import read_csv
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rcParams['figure.dpi'] = 300


def simple_stokes(degree_arr, power_arr, spacing, angle_extent, N_choice = 0): 
    """
    Basic stokes calculator which takes a degree and corresponding power array. 
    Plots data but does not execute plt.show(), so subsequent commands to modify 
    the graph can be written after calling the function.

    Parameters
    ----------
    degree_arr : Numpy array. 
        QWP angle in degrees. Use np.asarray(LIST) if data is type list. 
        DESCRIPTION.
    power_arr : Numpy array.
        Power for given QWP angles. 
    spacing : Integer.
        Spacing between angles of the QWP. Requires equally spaced measurements.
    angle_extent : Integer.
        Range of QWP angles used. 180 or 360. Note: Measurements for the final 
        angle do not need to be measured (i.e. measurements for 22.5 spacing 
        only need to go up to 157.5 for N = 8), but angle_extent should be 
        multiples of 180.
    N_choice : Integer, optional. 
        Option to manually set N instead of using the value picked by the code. 
        Not recommended to change and is primarily for debugging purposes. 
        The default is 0.

    Returns
    -------
    use_dict : Dictionary.
        Dictionary to return all relevant values of the calculation.
        Keys: "A", "B", "C", "D", "I", "Q", "U", "V", "dop" (degree of polarization), 
            "R2" (R-squared), "R2bar" (Corrected R-squared accounting for fitting multiple variables), 
            "Angles" (returns input degree_arr, used to match for graphing), 
            "Experiment" (returns input power_arr, used to match for graphing),
            "Reconstruction" (fitting of points at matching input angle array),
            "Fit_rad" (Continuous angle array for fitted intensities in radians), 
            "Fit_I" (Fitted intensities), 
            "Rounded Printout" (Summary of relevant data, rounded to fit as a graph title), 
            "Unrounded Printout" (Summary of relevant data, unrounded and printed to console by default)
    """
    
    rad_arr = np.deg2rad(degree_arr)
    
    A = 0
    B = 0
    C = 0 
    D = 0 
    
    _N = int(angle_extent/spacing)

    N = len(power_arr)
    
    loop = True
    while loop:
        if N%_N !=0: 
            N -= 1
        if N%_N == 0: 
            loop = False
    
    if N_choice: 
        N = N_choice
    
    for _a in range(N):
        theta = rad_arr[_a]
        s2 = np.sin(2*theta)
        c4 = np.cos(4*theta)
        s4 = np.sin(4*theta)
        
        p = power_arr[_a]
        
        A += p 
        B += p * s2
        C += p * c4
        D += p * s4


    A = A * 2 / N
    B = B * 4 / N 
    C = C * 4 / N 
    D = D * 4 / N

    I = A - C
    Q = 2 * C
    U = 2 * D
    V = B
    
    dop = np.sqrt(Q**2 + U**2 + V**2)/I
    
    re_arr = np.array([])
    
    for _a in range(len(power_arr)): 
        theta = rad_arr[_a]
        s2 = np.sin(2*theta)
        c4 = np.cos(4*theta)
        s4 = np.sin(4*theta)
        
        re = 1/2 * (A + B*s2 + C*c4 + D*s4)
        re_arr = np.append(re_arr, re)


    Sr = 0
    St = 0
    average = sum(power_arr)/len(power_arr)

    for _a in range(len(power_arr)): 
        Sr += (power_arr[_a] - re_arr[_a])**2
        St += (power_arr[_a] - average)**2
        
    R2 = (St - Sr)/St
    R2bar = 1 - (1 - R2) * (N-1)/(N-4)  # Modified R2 to account for fitting multiple variables
    
    
    '''Full simulation fit'''
    range_ = (N * spacing)*np.pi/180

    sim_ang_arr = np.linspace(0, range_, 360)
    sim_I = 1/2 * (I + Q*(np.cos(2*sim_ang_arr))**2 + U*np.cos(2*sim_ang_arr)*np.sin(2*sim_ang_arr) + V*np.sin(2*sim_ang_arr))

    data_print_round = rf"N = {N},  $R^2_{{adj}}$ = {round(R2bar, 6)}, dop = {round(dop, 3)}" + "\n" + rf"$S_{{norm}}$ = [{round(I/I, 3)}, {round(Q/I, 3)}, {round(U/I, 3)}, {round(V/I, 3)}]"
    data_print_unround = rf"N = {N},  $R^2_{{adj}}$ = R2bar, dop = {dop}" + "\n" + rf"$S_{{norm}}$ = [{I/I}, {Q/I}, {U/I}, {V/I}]"

    print(data_print_unround)
    
    fig, ax = plt.subplots(subplot_kw={'projection': 'polar'}, figsize=(5,5), dpi = 300)


    ax.plot(np.deg2rad(degree_arr), power_arr, color = "k", label = "Experiment")
    ax.plot(sim_ang_arr, sim_I, color = "r", label = "Fit")
    plt.ylim(0, 1.1)

    plt.title(data_print_round)
    
    use_dict = {"A": A, "B": B, "C": C, "D": D, 
                "I": I, "Q": Q, "U": U, "V": V, 
                "dop": dop, 
                "R2": R2, "R2bar": R2bar, 
                "Angles": degree_arr, "Experiment": power_arr,
                "Reconstruction": re_arr, 
                "Fit_rad": sim_ang_arr, "Fit_I": sim_I, 
                "Rounded Printout": data_print_round, "Unrounded Printout": data_print_unround}
    
    return use_dict
    

def PM100D_stokes(path, filenames, spacing, angle_extent, flip = 1, data_mod = False):
    """
    Stokes calculator for csv files generated with a Thorlabs PM100D power meter.
    Plots data but does not execute plt.show(), so subsequent commands to modify 
    the graph can be written after calling the function.

    Parameters
    ----------
    path : String.
        Path for folder where CSV files are stored.
    filenames : List of strings.
        File names of .CSV files. Include .CSV. Example generation: 
            filenames = []
            for i in range(1, 9):
                filenames.append(f"DATA0{i}.CSV")
    spacing : Integer.
        Spacing between angles of the QWP. Requires equally spaced measurements.
    angle_extent : Integer.
        Range of QWP angles used. 180 or 360. Note: Measurements for the final 
        angle do not need to be measured (i.e. measurements for 22.5 spacing 
        only need to go up to 157.5 for N = 8), but angle_extent should be 
        multiples of 180.
    flip : 1 or -1, optional
        Parameter that switches the fitting from a 0 degree linear polarizer to 
        a 90 degree linear polarizer. The default is 1.
    data_mod : Array-like, optional
        Used to manually access and modify data if there is a csv file issues. 
        Typically should not be used and is for debugging purposes only.
        The default is False.

    Returns
    -------
    use_dict : Dictionary.
        Dictionary to return all relevant values of the calculation.
        Keys: "A", "B", "C", "D", "I", "Q", "U", "V", "dop" (degree of polarization), 
            "R2" (R-squared), "R2bar" (Corrected R-squared accounting for fitting multiple variables), 
            "Angles" (returns array of angles used, used to match for graphing), 
            "Spectra Dict" (returns dictionary of form "QWP angle (str)": [QWP angle (int), np array of values]),
            "Reconstruction" (fitting of points at matching input angle array),
            "Fit_rad" (Continuous angle array for fitted intensities in radians), 
            "Fit_I" (Fitted intensities), 
            "Rounded Printout" (Summary of relevant data, rounded to fit as a graph title), 
            "Unrounded Printout" (Summary of relevant data, unrounded and printed to console by default)

    """
    spectra_dict = {}
    angles = np.array([])
    
    I_len = []
    i = 0
    
    for f in filenames: 
        angle = (i) * spacing
        angles = np.append(angles, angle)
        spectra = read_csv(os.path.join(path, f), header=3, delim_whitespace=True)
        spectra = spectra.to_numpy()
        spectra = np.transpose(spectra)
        Intensity = spectra[0]*1000
        
        spectra_dict[str(angle)] = [angle, Intensity]
        I_len.append(len(Intensity))
        i += 1
    
    if data_mod: 
        angle_mod = str(data_mod[0])
        angle_replace = str(data_mod[1])
        spectra_dict[angle_mod][1] = spectra_dict[angle_replace][1]
        if len(data_mod) == 4: 
            I_len[data_mod[2]] = data_mod[3]
    
                
    rrange = min(I_len)

    A = 0
    B = 0
    C = 0 
    D = 0
    N = 0
        
    for k, v in spectra_dict.items(): 
        theta = np.deg2rad(v[0])
        s2 = np.sin(2*theta)
        c4 = np.cos(4*theta)
        s4 = np.sin(4*theta)
        
        for pi in range(rrange): 
            p = v[1][pi]
            A += p 
            B += p * s2
            C += p * c4
            D += p * s4
            N += 1


    A = A * 2 / N
    B = B * 4 / N 
    C = C * 4 / N 
    D = D * 4 / N

    I = A - C * flip
    Q = 2 * C * flip
    U = 2 * D * flip
    V = B * flip

    dop = np.sqrt(Q**2 + U**2 + V**2)/I

    re_arr = np.array([])

    _sum = 0
    _len = 0

    for k, v in spectra_dict.items(): 
        theta = np.deg2rad(v[0])
        s2 = np.sin(2*theta)
        c4 = np.cos(4*theta)
        s4 = np.sin(4*theta)
        
        re = 1/2 * (A + B*s2 + C*c4 + D*s4)
        re_arr = np.append(re_arr, re)
        
        _sum += np.sum(v[1])
        _len += len(v[1])


    Sr = 0
    St = 0
    average = _sum/_len

    re_arr_index = 0
    for k, v in spectra_dict.items():
        for i in v[1]: 
            Sr += (i - re_arr[re_arr_index])**2 
            St += (i - average)**2 
        re_arr_index += 1
        

    R2 = (St - Sr)/St
    R2bar = 1 - (1 - R2) * (N-1)/(N-4)  # Modified R2 to account for fitting multiple variables


    '''Full simulation fit'''

    sim_ang_arr = np.linspace(0, np.pi, 360)

    sim_I = 1/2 * (A + B*np.sin(2*sim_ang_arr) + C * np.cos(4*sim_ang_arr) + D * np.sin(4*sim_ang_arr))
    #sim_I = 1/2 * (I - Q*(np.cos(2*sim_ang_arr))**2 - U*np.cos(2*sim_ang_arr)*np.sin(2*sim_ang_arr) - V*np.sin(2*sim_ang_arr))

    data_print_round = rf"N = {N},  $R^2_{{adj}}$ = {round(R2bar, 6)}, dop = {round(dop, 3)}" + "\n" + rf"$S_{{norm}}$ = [{round(I/I, 3)}, {round(Q/I, 3)}, {round(U/I, 3)}, {round(V/I, 3)}]"
    data_print_unround = rf"N = {N},  $R^2_{{adj}}$ = R2bar, dop = {dop}" + "\n" + rf"$S_{{norm}}$ = [{I/I}, {Q/I}, {U/I}, {V/I}]"

    use_dict = {"A": A, "B": B, "C": C, "D": D, 
                "I": I, "Q": Q, "U": U, "V": V, 
                "dop": dop, 
                "R2": R2, "R2bar": R2bar, 
                "Angles": angles, "Spectra Dict": spectra_dict,
                "Reconstruction": re_arr, 
                "Fit_rad": sim_ang_arr, "Fit_I": sim_I, 
                "Rounded Printout": data_print_round, "Unrounded Printout": data_print_unround}  

    fig, ax = plt.subplots(subplot_kw={'projection': 'polar'}, figsize=(5,5), dpi = 300)
    
    for k, v in spectra_dict.items(): 
        ax.scatter(v[1]*0 + np.deg2rad(v[0]), v[1], color = "k", label = "Experiment", lw = 0.01)
    
    ax.plot(sim_ang_arr, sim_I, color = "r", label = "Fit")
    
    plt.title(data_print_round)
    
    print(data_print_unround)
    
    return use_dict