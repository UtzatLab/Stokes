# -*- coding: utf-8 -*-
"""
Created on Fri Apr 14 14:35:24 2023

@author: alexm
"""

import numpy as np

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
    

class Stokes:
    def __init__(self, I, Q, U, V): 
        self.I = I 
        self.Q = Q
        self.U = U
        self.V = V
        self.S = np.array([[self.I], [self.Q], [self.U], [self.V]])