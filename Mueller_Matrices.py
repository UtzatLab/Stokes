# -*- coding: utf-8 -*-
"""
Created on Fri Apr 14 14:35:24 2023

@author: alexm
"""

import numpy as np

"""
Python object for common Mueller matrices. Example use: 
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
"""

class Mueller: 
    def __init__(self): 
        self.M = np.array([[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]])
    
    
    def set_full_matrix_list(self, elements_list):
        self.M = np.reshape(elements_list, (4,4))  
        
        
    def set_full_matrix(self, elements):
        self.M = elements
    
    
    def linpol(self, fax_angle_deg): 
        angle = np.deg2rad(fax_angle_deg)
        c = np.cos(2*angle)
        s = np.sin(2*angle)
        elements = 1/2 * np.array([[1,c,s,0],[c,c**2, c*s,0],[s,c*s,s**2, 0],[0,0,0,0]])
        self.set_full_matrix(elements)
        
        
    def qwp(self, fax_angle_deg, ret_deg):
        angle = np.deg2rad(fax_angle_deg)
        ret_rad = np.deg2rad(ret_deg)
        c = np.cos(2*angle)
        s = np.sin(2*angle)
        rc = np.cos(ret_rad)
        rs = np.sin(ret_rad)
        elements = np.array([[1,0,0,0],
                    [0,c**2+s**2*rc, c*s*(1-rc),s*rs],
                    [0,c*s*(1-rc),c**2*rc+s**2, -c*rs],
                    [0,-s*rs,c*rs,rc]])
        self.set_full_matrix(elements)
        
    def lens(self, a=1, b=0.995, c=-0.01, d=0, e=1, f=0.99):
        elements = np.array(
                    [[a,0,0,0],
                    [0,b,0,c],
                    [0,0,e,d],
                    [0,-c,d,f]])
        self.set_full_matrix(elements)
        
    def mirror(self, Rs = 0.98, Rp = 0.955, dphase = 0): 
        elements = 1/2 * np.array(
                    [[(Rs+Rp),(Rs-Rp),0,0],
                    [(Rs-Rp),(Rs+Rp),0,0],
                    [0,0, -2*np.sqrt(Rs*Rp)*np.cos(dphase),-2*np.sqrt(Rs*Rp)*np.sin(dphase)],
                    [0,0,2*np.sqrt(Rs*Rp)*np.sin(dphase),-2*np.sqrt(Rs*Rp)*np.cos(dphase)]])
        self.set_full_matrix(elements)
        
    def fmirror(self, Rs = 0.98, Rp = 0.955, dphase = 0): 
        elements = 1/2 * np.array(
                    [[(Rs+Rp),(Rp-Rs),0,0],
                    [(Rp-Rs),(Rs+Rp),0,0],
                    [0,0, -2*np.sqrt(Rs*Rp)*np.cos(-dphase),-2*np.sqrt(Rs*Rp)*np.sin(-dphase)],
                    [0,0,2*np.sqrt(Rs*Rp)*np.sin(-dphase),-2*np.sqrt(Rs*Rp)*np.cos(-dphase)]])
        self.set_full_matrix(elements)
        
    

class Stokes:
    def __init__(self, I, Q, U, V): 
        self.I = I 
        self.Q = Q
        self.U = U
        self.V = V
        self.S = np.array([[self.I], [self.Q], [self.U], [self.V]])