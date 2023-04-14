# -*- coding: utf-8 -*-
"""
Created on Thu Mar 23 18:25:22 2023

@author: alexm
"""

import numpy as np
from pandas import read_csv
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid.inset_locator import InsetPosition
import matplotlib as mpl

mpl.rcParams['figure.dpi'] = 300


class CPL_Analysis: 
    def __init__(self, filename, deg_spacing, angle_extent = 180, name = None, LPA = 1): 
        """
        Class for analyzing Stokes parameters and CPL from files generated 
        by Ocean Optics.

        Parameters
        ----------
        filename : String.
            Path and filename for .csv file containing combined spectra for 
            QWP angles. Example file displayed in GitHub folder.
        deg_spacing : Float.
            Spacing between angles of the QWP. Requires equally spaced measurements.
        angle_extent : 180 or 360, optional.
            Range of QWP angles used. Note: Measurements for the final 
            angle do not need to be measured (i.e. measurements for 22.5 spacing 
            only need to go up to 157.5 for N = 8), but angle_extent should be 
            multiples of 180.
            The default is 180.
        name : String, optional
            Identifier for data file. Can be called to appear on graphs or printouts. 
            Primarily used to keep track when analyzing multiple files at once. 
            The default is None.
        LPA : 1 or -1, optional
            Parameter that switches the fitting from a 0 degree linear polarizer to 
            a 90 degree linear polarizer. The default is 1.

        Returns
        -------
        None.


        TODO: Add options for headers on CSV file. Currently requires headerless data.
        """

        self.name = name
        self.spectra = read_csv(filename, header=None)

        self.spectra = self.spectra.to_numpy()
        self.spectra = np.transpose(self.spectra)
        self.wavelengths = self.spectra[0]

        '''Removing duplicate wavelengths from data'''
        n = []

        for i in range(int(len(self.spectra))): 
            if self.spectra[i][0] != self.wavelengths[0]: 
                n.append(self.spectra[i])

        self.spectra = np.asarray(n)
        self.spectraT = np.transpose(self.spectra)
        
        
        self.deg_spacing = deg_spacing
        self.angle_extent = angle_extent
        
        self.count_int = int(self.angle_extent/self.deg_spacing)
        
        
        self.N = len(self.spectra)
        '''N must be an integer multiple of self.count_int, so remove points if not.'''
        run = 1
        count = 0
        while run: 
            if self.N%self.count_int != 0: 
                self.N -= 1
                count += 1
            run = 0
        print(f"Removed {count} data points based on degree spacing.")

        self.flip = LPA
        
        self.has_run_IntegratedIQUV = 0
        self.has_run_SpectrumIQUV = 0


    def plot_spectrum(self, angle_index = 0, xlim = False, ylim = False, vlines = False, legend = True): 
        """
        Plot representative spectrum from csv file.

        Parameters
        ----------
        angle_index : Integer, optional
            Index of spectrum to plot. The default is 0 (1st spectrum).
        xlim : Tuple of floats, optional
            Manually set wavelength axis limits. The default is False.
        ylim : Tuple of floats, optional
            Manually set intensity axis limits. The default is False.
        vlines : Tuple of float(s), optional
            Plot verticle line at given value. The default is False.
        legend : True or False, optional
            Toggles legend. The default is True.

        Returns
        -------
        None.

        """

        fig, ax = plt.subplots(figsize=(6,5))
        plt.plot(self.wavelengths, self.spectra[angle_index], label = self.name)
        plt.ylabel("Intensity (counts)", size = 14)
        plt.xlabel("Wavelength (nm)", size = 14)
        
        plt.xlim(self.wavelengths[0], self.wavelengths[-1])
        if xlim: 
            plt.xlim(xlim[0], xlim[1])
        
        if ylim: 
            plt.ylim(ylim[0], ylim[1])
        
        if vlines: 
            for v in vlines: 
                plt.axvline(x=v, color = "silver", ls = "--")
                
        if legend: 
            plt.legend()
        
        plt.show()
        
    
    def plot_all_spectrum(self, xlim = False, ylim = False, vlines = False, Inset = False, inset_xlim = False, inset_ylim = False, inset_yb = 0.01, inset_xb = 0.005, colormap = 'twilight_shifted'): 
        """
        Plots all spectra from csv file.

        Parameters
        ----------
        xlim : Tuple of floats, optional
            Manually set wavelength axis limits. The default is False.
        ylim : Tuple of floats, optional
            Manually set intensity axis limits. The default is False.
        vlines : Tuple of float(s), optional
            Plot verticle line at given value. The default is False.
        Inset : True or False, optional
            Toggle inset. The default is False.
        inset_xlim : Tuple of floats, optional
            Manually set inset wavelength axis limits. The default is False.
        inset_ylim : Tuple of floats, optional
            Manually set inset intensity axis limits. The default is False.
        inset_yb : Float, optional
            Factor that "squeezes" the given inset_ylim to fit inside the figure. 
            The lower limit is divided by this, and the upper limit is multiplied.
            The default is 1.01.
        inset_xb : Float, optional
            Factor that "squeezes" the given inset_xlim to fit inside the figure. 
            The lower limit is divided by this, and the upper limit is multiplied.
            The default is 1.005.
        colormap : String from matplotlib colormaps, optional
            Manually set colormap. The default is 'twilight_shifted'.

        Returns
        -------
        None.

        """
        fig, ax = plt.subplots(figsize=(6,5))
        
        n = len(self.spectra)
        
        cmap = mpl.colormaps[colormap]
        colors = cmap(np.linspace(0,1,n))
        
        
        for _a in range(len(self.spectra)): 
            plt.plot(self.wavelengths, self.spectra[_a], color = colors[_a])
            
        plt.ylabel("Intensity (a.u.)", size = 14)
        plt.xlabel("Wavelength (nm)", size = 14)
        
        plt.xlim(self.wavelengths[0], self.wavelengths[-1])
        if xlim: 
            plt.xlim(xlim[0], xlim[1])
        
        if ylim: 
            plt.ylim(ylim[0], ylim[1])
        
        if vlines: 
            for v in vlines: 
                plt.axvline(x=v, color = "silver", ls = "--")
        
        if Inset: 
        
            ax_inset = fig.add_subplot(111)
            ip = InsetPosition(ax, [0.5,0.5,0.5,0.5])
            ax_inset.set_axes_locator(ip)
            
            for _a in range(len(self.spectra)): 
                ax_inset.plot(self.wavelengths, self.spectra[_a], color = colors[_a])
            
            if inset_ylim:
                ax_inset.set_ylim(inset_ylim[0]/inset_yb, inset_ylim[1]*inset_yb)
                spacing = (inset_ylim[1]-inset_ylim[0])/4
                ax_inset.set_yticks(np.arange(inset_ylim[0], inset_ylim[1] + spacing, spacing))
                
            if inset_xlim:
                ax_inset.set_xlim(inset_xlim[0]/inset_xb, inset_xlim[1]*inset_xb)
                spacing = (inset_xlim[1]-inset_xlim[0])/3
                ax_inset.set_xticks(np.arange(inset_xlim[0], inset_xlim[1] + spacing, spacing))
                
            if inset_ylim:
                f0 = 1
                f1 = 1
                if inset_ylim[0] < 0: 
                    f0 = -1
                if inset_ylim[1] < 0: 
                    f1 = -1
                ax_inset.set_ylim(inset_ylim[0] - f0 * (inset_ylim[0]*inset_yb), inset_ylim[1]+f1*(inset_ylim[1]*inset_yb))
                ax_inset.locator_params(axis='y', nbins=5)

            if inset_xlim:
                f0 = 1
                f1 = 1
                if inset_xlim[0] < 0: 
                    f0 = -1
                if inset_xlim[1] < 0: 
                    f1 = -1
                ax_inset.set_xlim(inset_xlim[0] - f0 * (inset_xlim[0]*inset_xb), inset_xlim[1]+f1*(inset_xlim[1]*inset_xb))
                ax_inset.locator_params(axis='x', nbins=5)
            
        plt.show()
        
    
    def IntegratedIQUV(self, Normalized=True, units = "mdeg", Bounds = False):
        """
        Stokes parameter calculations for an integrated segment of wavelengths.
        Called by self.plotIntegratedIQUV, so it is unnecessary to call this 
        function before the other.

        Parameters
        ----------
        Normalized : True or False or Integer, optional
            Toggles if the data is normalized or not. If an integer is supplied, 
            the data will be normalized to the intensity at that index. This is 
            used if the max intensity is not where you want to normalize the data 
            (i.e. if the excition signal is greater than emission signal).
            The default is True.
        units : "mdeg", "deg", "rad", or "mrad", optional
            Units for CPL to be displayed in. The default is "mdeg".
        Bounds : Tuple of Integers, optional
            Gives wavelength bounds to integrate between. Must be supplied as 
            the index of the corresponding wavelength, not as the wavelength 
            itself. This can be checked in self.wavelengths. The default is False.

        Returns
        -------
        None.

        """
        
        self.angle_arr = np.asarray([i*self.deg_spacing for i in range(len(self.spectra))])
        self.angle_arr = np.deg2rad(self.angle_arr)
        self.power_arr = np.array([])

        for _a in self.spectra: 
            if Bounds: 
                power_sum = 0
                for i in range(Bounds[0], Bounds[1]): 
                    power_sum += _a[i]
                self.power_arr = np.append(self.power_arr, power_sum)
            if not Bounds: 
                self.power_arr = np.append(self.power_arr, np.sum(_a))
            
        if Normalized: 
            self.power_arr /= np.max(self.power_arr)
            if Normalized is int: 
                self.power_arr /= self.power_arr[Normalized]

        self.A_int = 0
        self.B_int = 0
        self.C_int = 0
        self.D_int = 0
        
        for _a in range(self.N): 
            theta = self.angle_arr[_a]
            s2 = np.sin(2*theta)
            c4 = np.cos(4*theta)
            s4 = np.sin(4*theta)
            
            p = self.power_arr[_a]
            
            self.A_int += p
            self.B_int += p*s2
            self.C_int += p*c4
            self.D_int += p*s4


        self.A_int = self.A_int * 2 / self.N
        self.B_int = self.B_int * 4 / self.N
        self.C_int = self.C_int * 4 / self.N
        self.D_int = self.D_int * 4 / self.N

        self.I_int = self.A_int - (self.C_int * self.flip)
        self.Q_int = 2 * self.C_int * self.flip
        self.U_int = 2 * self.D_int * self.flip
        self.V_int = self.B_int * self.flip

        self.P_int = np.sqrt(self.Q_int**2 + self.U_int**2 + self.V_int**2)/self.I_int

        
        '''Reconstructing the Intensity to check the fit'''
        self.reInt_arr = np.array([])
        for _a in range(len(self.angle_arr)): 
            theta = self.angle_arr[_a]
            s2 = np.sin(2*theta)
            c4 = np.cos(4*theta)
            s4 = np.sin(4*theta)
            reInt = 1/2 * (self.A_int + self.B_int*s2 + self.C_int*c4 + self.D_int*s4)
            self.reInt_arr = np.append(self.reInt_arr, reInt)

        '''Full simulation fit'''
        self.sim_ang_arr = np.linspace(0, 2*np.pi, 360)
        self.sim_I = 1/2 * (self.I_int 
                       + self.Q_int*(np.cos(2*self.sim_ang_arr))**2 * self.flip
                       + self.U_int*np.cos(2*self.sim_ang_arr)*np.sin(2*self.sim_ang_arr) * self.flip
                       + self.V_int*np.sin(2*self.sim_ang_arr)) * self.flip
        
        '''Least squares fit'''
        Sr = 0
        St = 0
        self.average = sum(self.power_arr)/len(self.power_arr)
        
        for _a in range(len(self.angle_arr)): 
            Sr += (self.power_arr[_a] - self.reInt_arr[_a])**2
            St += (self.power_arr[_a] - self.average)**2
            
        self.R2 = (St - Sr)/St
            
        
        self.CPL_int = (-1/2 * np.arcsin(self.V_int/self.I_int))
        
        if units == "mdeg" or units == "deg":
            self.CPL_int = np.rad2deg(self.CPL_int)
        if units == "mdeg" or units == "mrad":
            self.CPL_int *= 1000
        
        supported_units = ["mdeg", "deg", "rad", "mrad"]
        if units not in supported_units:
            print("Currently supported units: mdeg, deg, rad, mrad. Enter as a string.")
            return
        

        R2_string = r"$R^2$"
        
        self.Title_unrounded = f"S = [{self.I_int}, {self.Q_int}, {self.U_int}, {self.V_int}]  P = {self.P_int}  CPL = {self.CPL_int} {units}, {R2_string} = {self.R2}"
        self.Title = f"S = [{round(self.I_int, 3)}, {round(self.Q_int, 3)}, {round(self.U_int, 3)}, {round(self.V_int, 3)}]  P = {round(self.P_int, 4)}  CPL = {round(self.CPL_int, 3)} {units}"
        
        self.Title_unrounded_norm = rf"S = [{self.I_int/self.I_int}, {self.Q_int/self.I_int}, {self.U_int/self.I_int}, {self.V_int/self.I_int}]  P = {round(self.P_int, 7)}, $R^{{2}}$ = {round(self.R2, 7)}"
        self.Title_norm = rf"S = [{round(self.I_int/self.I_int, 3)}, {round(self.Q_int/self.I_int, 3)}, {round(self.U_int/self.I_int, 3)}, {round(self.V_int/self.I_int, 3)}]  P = {round(self.P_int, 3)}, $R^{{2}}$ = {round(self.R2, 4)}"
        
        #self.Title_unrounded = rf"S = [{{{self.I_int}}}, {{{self.Q_int}}}, {{{self.U_int}}}, {{{self.V_int}}}]  P = {{{self.P_int}}}  CPL = {{{self.CPL_int}}} {{{units}}}, $R^2$ = {{{self.R2}}}"
        
        print(self.Title_unrounded_norm)
        
        self._has_run_IntegratedIQUV = 1
        
        
    def plotIntegratedIQUV(self, Normalized = True, Bounds = False, Experimental=True, Reconstructed=False, Simulated=True, ylim=False, legend=True, Stokes = True, title = True): 
        """
        Plots polar data from self.IntegratedIQUV. Will call self.IntegratedIQUV, 
        so it is unnecessary to call that function before this one. 

        Parameters
        ----------
        Normalized : True or False or Integer, optional
            Toggles if the data is normalized or not. If an integer is supplied, 
            the data will be normalized to the intensity at that index. This is 
            used if the max intensity is not where you want to normalize the data 
            (i.e. if the excition signal is greater than emission signal).
            The default is True.
        Bounds : Tuple of Integers, optional
            Gives wavelength bounds to integrate between. Must be supplied as 
            the index of the corresponding wavelength, not as the wavelength 
            itself. This can be checked in self.wavelengths. The default is False.
        Experimental : True or False, optional
            Toggles plotting the experimental data. The default is True.
        Reconstructed : True of False, optional
            Toggles plotting the reconstructed data points. The default is False.
        Simulated : True or False, optional
            Toggles plotting the full fit data. The default is True.
        ylim : Tuple of floats, optional
            Manually set y limit. The default is False.
        legend : True or False, optional
            Toggles legend. The default is True.
        Stokes : True of False, optional
            Toggles title containing Stokes parameters. The default is True.
        title : True or False, optional
            Toggles title containing sample name and Stokes parameters. 
            Overwrites Stokes if True. The default is False.

        Returns
        -------
        None.

        """
        self.IntegratedIQUV(Normalized=Normalized, Bounds = Bounds)
            
        fig, ax = plt.subplots(subplot_kw={'projection': 'polar'}, figsize=(10,10))
        if Experimental:
            ax.plot(self.angle_arr, self.power_arr, color = "c", label = "Experimental", lw = 3)
        if Reconstructed: 
            ax.plot(self.angle_arr, self.reInt_arr, color = "k", label = "Reconstruction", lw = 3)
        if Simulated:
            ax.plot(self.sim_ang_arr, self.sim_I, color = "k", label = "Simulated", ls = "dashed", lw = 3)
            
            
        ax.tick_params(labelsize = 24, pad = 18)
        ax.tick_params(axis="y", labelsize = 20, colors = "gray")
        

        if legend: 
            plt.legend()
            angle = np.deg2rad(67.5)
            ax.legend(loc="lower left", bbox_to_anchor=(.5 + np.cos(angle)/2, .5 + np.sin(angle)/2), fontsize = 24)
        
        if ylim: 
            ymin, ymax = ylim[0], ylim[1]
            plt.ylim(ymin, ymax)
            
        if Stokes:
            ax.set_title(f"{self.Title_norm}", y = -.13, x=.5, fontsize = 24)
            
        if title: 
            ax.set_title(f"Sample: {self.name}\n{self.Title_norm}", y = -.2, x=.5, fontsize = 24)
            
        plt.show()
    
    
    
    def singleWvlIQUV(self, index, Normalized = True, units = "mdeg"):
        """
        Stokes parameter calculations for a single emission wavelength.
        Called by self.plotsingleWvlIQUV, so it is unnecessary to call this 
        function before the other.
        

        Parameters
        ----------
        index : Integer. 
            Index of desired wavelength. This can be checked in self.wavelengths.
        Normalized : True or False, optional
            Toggle normalizing data for stokes parameters. The default is True.
        units : "mdeg", "deg", "rad", or "mrad", optional
            Units for CPL to be displayed in. The default is "mdeg".

        Returns
        -------
        None.

        """
        
        self.angle_arr = np.asarray([i*self.deg_spacing for i in range(len(self.spectra))])
        self.angle_arr = np.deg2rad(self.angle_arr)
        self.power_arr = self.spectraT[index]
        

        if Normalized: 
            print("Normalized")
            self.power_arr /= np.max(self.power_arr)
        
                
        self.A_wvl = 0
        self.B_wvl = 0
        self.C_wvl = 0
        self.D_wvl = 0

        for _a in range(self.N): 
            angle = _a * self.deg_spacing
            theta = np.deg2rad(angle)
            s2 = np.sin(2*theta)
            c4 = np.cos(4*theta)
            s4 = np.sin(4*theta)
            
            p = self.power_arr[_a]
            
            self.A_wvl += p
            self.B_wvl += p*s2
            self.C_wvl += p*c4
            self.D_wvl += p*s4
                    
        self.A_wvl = self.A_wvl * 2 / self.N
        self.B_wvl = self.B_wvl * 4 / self.N
        self.C_wvl = self.C_wvl * 4 / self.N
        self.D_wvl = self.D_wvl * 4 / self.N


        self.I_wvl = self.A_wvl - self.C_wvl
        self.Q_wvl = 2 * self.C_wvl
        self.U_wvl = 2 * self.D_wvl
        self.V_wvl = self.B_wvl
        
        self.P_wvl = np.sqrt(self.Q_wvl**2 + self.U_wvl**2 + self.V_wvl**2)/self.I_wvl

        
        '''Reconstructing the Intensity to check the fit'''
        self.reWvl_arr = np.array([])
        for _a in range(len(self.angle_arr)): 
            theta = self.angle_arr[_a]
            s2 = np.sin(2*theta)
            c4 = np.cos(4*theta)
            s4 = np.sin(4*theta)
            reWvl = 1/2 * (self.A_wvl + self.B_wvl*s2 + self.C_wvl*c4 + self.D_wvl*s4)
            self.reWvl_arr = np.append(self.reWvl_arr, reWvl)

        '''Full simulation fit'''
        self.sim_ang_arr = np.linspace(0, 2*np.pi, 360)
        self.sim_I = 1/2 * (self.I_wvl
                       + self.Q_wvl*(np.cos(2*self.sim_ang_arr))**2 
                       + self.U_wvl*np.cos(2*self.sim_ang_arr)*np.sin(2*self.sim_ang_arr) 
                       + self.V_wvl*np.sin(2*self.sim_ang_arr))
        
        
        self.CPL_wvl = -1/2 * np.arcsin(self.V_wvl/self.I_wvl)
                
        if units == "mdeg" or units == "deg":
            self.CPL_wvl = np.rad2deg(self.CPL_wvl)
        if units == "mdeg" or units == "mrad":
            self.CPL_wvl *= 1000
            
        supported_units = ["mdeg", "deg", "rad", "mrad"]
        if units not in supported_units:
            print("Currently supported units: mdeg, deg, rad, mrad. Enter as a string.")
            return
        
        self.Title_unrounded = f"S = [{self.I_wvl}, {self.Q_wvl}, {self.U_wvl}, {self.V_wvl}]  P = {self.P_wvl}  CPL = {self.CPL_wvl} {units}"
        self.Title = f"S = [{round(self.I_wvl, 3)}, {round(self.Q_wvl, 3)}, {round(self.U_wvl, 3)}, {round(self.V_wvl, 3)}]  P = {round(self.P_wvl, 4)}  CPL = {round(self.CPL_wvl, 3)} {units}"
        print(self.Title_unrounded)
        
        
    def plotsingleWvlIQUV(self, index, Normalized = True, units = "mdeg", Experimental=True, Reconstructed=False, Simulated=True, ylim=False, legend=True, Stokes = True, title = False):
        """
        Plots polar data from self.singleWvlIQUV. Will call self.singleWvlIQUV, 
        so it is unnecessary to call that function before this one. 

        Parameters
        ----------
        index : Integer. 
            Index of desired wavelength. This can be checked in self.wavelengths.
        Normalized : True or False, optional
            Toggle normalizing data for stokes parameters. The default is True.
        units : "mdeg", "deg", "rad", or "mrad", optional
            Units for CPL to be displayed in. The default is "mdeg".
        Experimental : True or False, optional
            Toggles plotting the experimental data. The default is True.
        Reconstructed : True of False, optional
            Toggles plotting the reconstructed data points. The default is False.
        Simulated : True or False, optional
            Toggles plotting the full fit data. The default is True.
        ylim : Tuple of floats, optional
            Manually set y limit. The default is False.
        legend : True or False, optional
            Toggles legend. The default is True.
        Stokes : True of False, optional
            Toggles title containing Stokes parameters. The default is True.
        title : True or False, optional
            Toggles title containing sample name and Stokes parameters. 
            Overwrites Stokes if True. The default is False.

        Returns
        -------
        None.
        
        TODO: Clean up animation function.

        """
        self.singleWvlIQUV(index, Normalized = Normalized, units = units)
        
        wvl = self.wavelengths[index]
        
        fig, ax = plt.subplots(subplot_kw={'projection': 'polar'}, figsize=(10,10))

        
        if Experimental:
            ax.plot(self.angle_arr, self.power_arr, color = "c", label = "Experimental")
        if Reconstructed: 
            ax.plot(self.angle_arr, self.reWvl_arr, color = "k", label = "Reconstruction")
        if Simulated:
            ax.plot(self.sim_ang_arr, self.sim_I, color = "gray", label = "Simulated", ls = "dashed")
            #plt.fill_between(self.sim_ang_arr, self.sim_I, color = "gray", alpha=0.4)
            
        ax.tick_params(labelsize = 14)
        ax.tick_params(axis="y", labelsize = 12, colors = "gray")
        

        if legend: 
            plt.legend()
            angle = np.deg2rad(67.5)
            ax.legend(loc="lower left", bbox_to_anchor=(.5 + np.cos(angle)/2, .5 + np.sin(angle)/2), fontsize = 14)
        
        if ylim: 
            ymin, ymax = ylim[0], ylim[1]
            plt.ylim(ymin, ymax)
            
        if Stokes:
            ax.set_title(f"{self.Title}", y = -.05, x=.5, fontsize = 16)
            
        if title: 
            lam = r"$\lambda$"
            ax.set_title(f"Sample: {self.name}\n{lam} = {round(wvl,3)} nm\n{self.Title}", y = -.15, x=.5, fontsize = 16)
            
        plt.show()
        
        # plt.savefig('frame_'+str(index)+'_.png',dpi=300)
    
    
    def SpectrumIQUV(self, boxcar=0, Normalized = True): 
        """
        Stokes parameter calculations by wavelength.
        Called by self.plotSpectrumIQUV, so it is unnecessary to call this 
        function before the other.

        Parameters
        ----------
        boxcar : Integer, optional
            Boxcar smoothing width. The default is 0.
        Normalized : True or False or Integer, optional
            Toggles if the Stokes parameters are normalized. If an integer is supplied, 
            the data will be normalized to the intensity at that index. This is 
            used if the max intensity is not where you want to normalize the data 
            (i.e. if the excition signal is greater than emission signal).
            The default is True.

        Returns
        -------
        None.

        """
        
        self.A_arr = np.array([])
        self.B_arr = np.array([])
        self.C_arr = np.array([])
        self.D_arr = np.array([])

        for wvl in range(len(self.wavelengths)): 
            #Iterate wavelengths: 
            A = 0
            B = 0 
            C = 0
            D = 0
            for _a in range(self.N): 
                angle = _a * (self.deg_spacing)
                theta = np.deg2rad(angle)
                s2 = np.sin(2*theta)
                c4 = np.cos(4*theta)
                s4 = np.sin(4*theta)
                
                p = self.spectra[_a][wvl]
                
                A += p
                B += p*s2
                C += p*c4
                D += p*s4
                        
            A = A * 2 / self.N
            B = B * 4 / self.N
            C = C * 4 / self.N
            D = D * 4 / self.N
            
            self.A_arr = np.append(self.A_arr, A)
            self.B_arr = np.append(self.B_arr, B)
            self.C_arr = np.append(self.C_arr, C)
            self.D_arr = np.append(self.D_arr, D)


        self.I = self.A_arr - self.C_arr
        self.Q = 2 * self.C_arr
        self.U = 2 * self.D_arr
        self.V = self.B_arr
        
        
        #Boxcar smoothing 
        if boxcar: 
            kernel_size = boxcar
            kernel = np.ones(kernel_size) / kernel_size
    
            self.I = np.convolve(self.I, kernel, mode='same')
            self.Q = np.convolve(self.Q, kernel, mode='same')
            self.U = np.convolve(self.U, kernel, mode='same')
            self.V = np.convolve(self.V, kernel, mode='same')


        self.P = np.sqrt(self.Q**2 + self.U**2 + self.V**2)/self.I

        self.I_norm = self.I/np.max(self.I)
        self.Q_norm = self.Q/np.max(self.I)
        self.U_norm = self.U/np.max(self.I)
        self.V_norm = self.V/np.max(self.I)
        
        
        self.R2_arr = np.array([])
        
        for wvl in range(len(self.wavelengths)): 
            #Iterate wavelengths: 
            self.reInt_arr_s = np.array([])
            for _a in range(self.N): 
                angle = _a * (self.deg_spacing)
                theta = np.deg2rad(angle)
                s2 = np.sin(2*theta)
                c4 = np.cos(4*theta)
                s4 = np.sin(4*theta)
                
                reInt = 1/2 * (A + B*s2 + C*c4 + D*s4)
                
                self.reInt_arr_s = np.append(self.reInt_arr_s, reInt)
                
            '''Least squares fit'''
            Sr = 0
            St = 0
            self.average = sum(self.spectraT[wvl])/len(self.spectraT[wvl])
            
            for _a in range(self.N): 
                Sr += (self.spectraT[wvl][_a] - self.reInt_arr_s[_a])**2
                St += (self.spectraT[wvl][_a] - self.average)**2
                
            R2 = (St - Sr)/St
            self.R2_arr = np.append(self.R2_arr, R2)
        
        
        if isinstance(Normalized, int) and Normalized is not True:
            print(Normalized)
            self.div = self.I[Normalized]
    
            print(self.div)
            self.I_norm = self.I/self.div
            self.Q_norm = self.Q/self.div
            self.U_norm = self.U/self.div
            self.V_norm = self.V/self.div


    def plotSpectrumIQUV(self, boxcar = 0, Normalized = True, ylim = False, xlim = False, legend = True, title = False, Inset = False, inset_ylim = False, inset_xlim = False, inset_yb = 0.1, inset_xb = 0.025, inset_y_size = 0.2, inset_x_size = 0): 
        """
        Plots data from self.SpectrumIQUV. Will call self.SpectrumIQUV, 
        so it is unnecessary to call that function before this one. 

        Parameters
        ----------
        boxcar : Integer, optional
            Boxcar smoothing width. The default is 0.
        Normalized : True or False or Integer, optional
            Toggles if the Stokes parameters are normalized. If an integer is supplied, 
            the data will be normalized to the intensity at that index. This is 
            used if the max intensity is not where you want to normalize the data 
            (i.e. if the excition signal is greater than emission signal).
            The default is True.
        ylim : Tuple of floats, optional
            Manually set y limit. The default is False.
        xlim : Tuple of floats, optional
            Manually set x limit. The default is False.
        legend : True or False, optional
            Toggles legend. The default is True.
        title : True or False, optional
            Toggles title containing sample name. 
            Overwrites Stokes if True. The default is False.
        Inset : True or False, optional
            Toggles inset. Must enable to use inset modifiers. The default is False.
        inset_ylim : Tuple of floats, optional
            Manually set inset y limit. The default is False.
        inset_xlim : Tuple of floats, optional
            Manually set inset x limit. The default is False.
        inset_yb : Float, optional
            Factor that "squeezes" the given inset_ylim to fit inside the figure. 
            The lower limit is divided by this, and the upper limit is multiplied.
            The default is 0.1.
        inset_xb : Float, optional
            Factor that "squeezes" the given inset_xlim to fit inside the figure. 
            The lower limit is divided by this, and the upper limit is multiplied.
            The default is 0.025.
        inset_y_size : Float between -0.5 and 0.5, optional
            Change the size of the inset size in the y-direction. The default is 0.2.
        inset_x_size : Float between -0.5 and 0.5, optional
            Change the size of the inset size in the x-direction. The default is 0.

        Returns
        -------
        None.

        """

        self.SpectrumIQUV(boxcar = boxcar, Normalized = Normalized)
        
        fig, ax = plt.subplots(figsize=(6,5))
        
        if Normalized: 
            plt.plot(self.wavelengths, self.I_norm, color = "k", label = "I", linewidth = 2)
            plt.plot(self.wavelengths, self.Q_norm, color = "r", label = "Q", linewidth = 2)
            plt.plot(self.wavelengths, self.U_norm, color = "g", label = "U", linewidth = 2)
            plt.plot(self.wavelengths, self.V_norm, color = "b", label = "V", linewidth = 2)
            plt.ylabel("Normalized Stokes Parameters", size = 14)
        if not Normalized: 
            plt.plot(self.wavelengths, self.I, color = "k", label = "I", linewidth = 2)
            plt.plot(self.wavelengths, self.Q, color = "r", label = "Q", linewidth = 2)
            plt.plot(self.wavelengths, self.U, color = "g", label = "U", linewidth = 2)
            plt.plot(self.wavelengths, self.V, color = "b", label = "V", linewidth = 2)
            plt.ylabel("Unnormalized Stokes Parameters", size = 14)
            
        plt.xlabel("Wavelength (nm)", size = 14)
        
        if ylim: 
            ymin, ymax = ylim[0], ylim[1]
            plt.ylim(ymin, ymax)
        
        if xlim: 
            xmin, xmax = xlim[0], xlim[1]
            plt.xlim(xmin, xmax)
            
        if title: 
            plt.title(label = self.name)
            
        if legend: 
            plt.legend(loc = "upper right")
            if Inset: 
                plt.legend(loc = "upper left")
            
        if Inset:        
            ax_inset = fig.add_subplot(111)
            ip = InsetPosition(ax, [0.5-inset_x_size,0.5-inset_y_size,0.5+inset_x_size,0.5+inset_y_size])
            ax_inset.set_axes_locator(ip)
            
            if Normalized:
                ax_inset.plot(self.wavelengths, self.I_norm, color = "k", label = "I", linewidth = 2)
                ax_inset.plot(self.wavelengths, self.Q_norm, color = "r", label = "Q", linewidth = 2)
                ax_inset.plot(self.wavelengths, self.U_norm, color = "g", label = "U", linewidth = 2)
                ax_inset.plot(self.wavelengths, self.V_norm, color = "b", label = "V", linewidth = 2)
                
            if not Normalized: 
                ax_inset.plot(self.wavelengths, self.I, color = "k", label = "I", linewidth = 2)
                ax_inset.plot(self.wavelengths, self.Q, color = "r", label = "Q", linewidth = 2)
                ax_inset.plot(self.wavelengths, self.U, color = "g", label = "U", linewidth = 2)
                ax_inset.plot(self.wavelengths, self.V, color = "b", label = "V", linewidth = 2)
                
            if inset_ylim:
                f0 = 1
                f1 = 1
                if inset_ylim[0] < 0: 
                    f0 = -1
                if inset_ylim[1] < 0: 
                    f1 = -1
                ax_inset.set_ylim(inset_ylim[0] - f0 * (inset_ylim[0]*inset_yb), inset_ylim[1]+f1*(inset_ylim[1]*inset_yb))
                ax_inset.locator_params(axis='y', nbins=5)

            if inset_xlim:
                f0 = 1
                f1 = 1
                if inset_xlim[0] < 0: 
                    f0 = -1
                if inset_xlim[1] < 0: 
                    f1 = -1
                ax_inset.set_xlim(inset_xlim[0] - f0 * (inset_xlim[0]*inset_xb), inset_xlim[1]+f1*(inset_xlim[1]*inset_xb))
                ax_inset.locator_params(axis='x', nbins=5)
                    
        plt.show()
        

    def plotSpectrumR2(self, ylim = False, xlim = False, legend = True, title = False, boxcar = False):
        """
        TODO: Currently basically inoperable. Fix this.

        Parameters
        ----------
        ylim : TYPE, optional
            DESCRIPTION. The default is False.
        xlim : TYPE, optional
            DESCRIPTION. The default is False.
        legend : TYPE, optional
            DESCRIPTION. The default is True.
        title : TYPE, optional
            DESCRIPTION. The default is False.
        boxcar : TYPE, optional
            DESCRIPTION. The default is False.

        Returns
        -------
        None.

        """
        self.SpectrumIQUV()
       
        fig, ax = plt.subplots(figsize=(6,5))
       
        plt.ylim(-1, 1)
       
        plt.plot(self.wavelengths, self.R2_arr)
        plt.show()
       

    def plotCPL(self, boxcar = 0, Intensity_weighted = True, units = "mdeg", axline = True, ylim = False, xlim = False, color = 'c', legend = False, title = False):       
        """
        Plot CPL spectrum as generated by self.SpectrumIQUV. Will call self.SpectrumIQUV, 
        so it is unnecessary to call that function before this one. 

        Parameters
        ----------
        boxcar : Integer, optional
            Boxcar smoothing width. The default is 0.
        Intensity_weighted : True, False, or "Both", optional
            Toggles if the CPL spectrum will be multiplied by self.I_norm. This 
            is to flatten the noise from small number division in the spectrum 
            where the sample does not emit. Entering "Both" will plot both on 
            the sample graph. The default is True.
        units : "mdeg", "deg", "rad", or "mrad", optional
            Units for CPL to be displayed in. The default is "mdeg".
        axline : True or False, optional
            Toggles horizontal line at CPL = 0. The default is True.
        xlim : Tuple of floats, optional
            Manually set wavelength axis limits. The default is False.
        ylim : Tuple of floats, optional
            Manually set CPL axis limits. The default is False.
        color : Matplotlib color, optional
            Change plot color. The default is 'c'.
        legend : True or false, optional
            Toggles legend. If Intensity_weighted is not both, this displays 
            the sample name. If Intensity_weighted is "Both", this labels the 
            weighted and unweighted lines. The default is False.
        title : True or False, optional
            Toggles title of sample name. The default is False.

        Returns
        -------
        None.

        """
        self.SpectrumIQUV(Normalized = True, boxcar = boxcar) 

        self.CPL = -1/2 * np.arcsin(self.V/self.I)
        
        if units == "mdeg" or units == "deg":
            self.CPL = np.rad2deg(self.CPL)
        if units == "mdeg" or units == "mrad":
            self.CPL *= 1000
        supported_units = ["mdeg", "deg", "rad", "mrad"]
        if units not in supported_units:
            print("Currently supported units: mdeg, deg, rad, mrad. Enter as a string.")
            return
        
        fig, ax = plt.subplots(figsize=(6,5))
        
        if Intensity_weighted and Intensity_weighted != "Both": 
            plt.ylabel(f"CPL ({units}) (Intensity weighted)", size = 14)
            self.CPL_Intensity = self.CPL * self.I_norm
            plt.plot(self.wavelengths, self.CPL_Intensity, color = color, label = self.name)
            
        if not Intensity_weighted: 
            plt.ylabel(f"CPL ({units})", size = 14)
            plt.plot(self.wavelengths, self.CPL, color = color, label = self.name)
            
        if Intensity_weighted == "Both": 
            plt.ylabel(f"CPL ({units})", size = 14)
            self.CPL_Intensity = self.CPL * self.I_norm
            plt.plot(self.wavelengths, self.CPL_Intensity, color = color, label = "Intensity weighted")
            plt.plot(self.wavelengths, self.CPL, color = color, alpha = 0.3, label = "Unweighted")

        if axline: 
            plt.plot(self.wavelengths, self.wavelengths * 0, color = "k", lw = 0.7)
        
        if ylim: 
            ymin, ymax = ylim[0], ylim[1]
            plt.ylim(ymin, ymax)
        
        if xlim: 
            xmin, xmax = xlim[0], xlim[1]
            plt.xlim(xmin, xmax)
            
        if legend: 
            plt.legend()
            
        if title: 
            plt.title(label = self.name)
        
        
        plt.xlabel("Wavelength (nm)", size = 14)
        
        plt.show()
        
    
    def plotglum(self, boxcar = False, Intensity_weighted = True, scaling = 3, axline = True, ylim = False, xlim = False, color = 'c', legend = False, title = False): 
        """
        Plot g_lum spectrum as generated by self.SpectrumIQUV. Will call self.SpectrumIQUV, 
        so it is unnecessary to call that function before this one. 

        Parameters
        ----------
        boxcar : Integer, optional
            Boxcar smoothing width. The default is 0.
        Intensity_weighted : True, False, or "Both", optional
            Toggles if the CPL spectrum will be multiplied by self.I_norm. This 
            is to flatten the noise from small number division in the spectrum 
            where the sample does not emit. Entering "Both" will plot both on 
            the sample graph. The default is True.
        scaling : Float, optional
            Multiplies g_lum by 10 to the power of this float. The default is 3.
        units : "mdeg", "deg", "rad", or "mrad", optional
            Units for CPL to be displayed in. The default is "mdeg".
        axline : True or False, optional
            Toggles horizontal line at CPL = 0. The default is True.
        xlim : Tuple of floats, optional
            Manually set wavelength axis limits. The default is False.
        ylim : Tuple of floats, optional
            Manually set CPL axis limits. The default is False.
        color : Matplotlib color, optional
            Change plot color. The default is 'c'.
        legend : True or false, optional
            Toggles legend. If Intensity_weighted is not both, this displays 
            the sample name. If Intensity_weighted is "Both", this labels the 
            weighted and unweighted lines. The default is False.
        title : True or False, optional
            Toggles title of sample name. The default is False.

        Returns
        -------
        None.

        """
        
        self.SpectrumIQUV(Normalized = True, boxcar = boxcar) 
        
        glum_IW = (- 2 * self.V_norm)*10**scaling
        glum = (- 2 * self.V/self.I)*10**scaling
        
        fig, ax = plt.subplots(figsize=(6,5))
        
        if Intensity_weighted and Intensity_weighted != "Both":
            plt.plot(self.wavelengths, glum_IW, color = color, label = self.name)
        
        if not Intensity_weighted:
            plt.plot(self.wavelengths, glum, color = color, label = self.name)
            
        if Intensity_weighted == "Both": 
            plt.plot(self.wavelengths, glum_IW, color = color, label = "Intensity Weighted")
            plt.plot(self.wavelengths, glum, color = color, alpha = 0.3, label = "Unweighted")

        if axline: 
            plt.plot(self.wavelengths, self.wavelengths * 0, color = "k", lw = 0.7)
        
        if ylim: 
            ymin, ymax = ylim[0], ylim[1]
            plt.ylim(ymin, ymax)
        
        if xlim: 
            xmin, xmax = xlim[0], xlim[1]
            plt.xlim(xmin, xmax)
            
        if legend: 
            plt.legend()
            
        if title: 
            plt.title(label = self.name)
        
        if scaling != 0: 
            plt.ylabel(rf"$g_{{lum}} \times 10^{{{scaling}}}$", size = 14)
        if scaling == 0: 
            plt.ylabel(r"$g_{lum}$", size = 14)
        plt.xlabel("Wavelength (nm)", size = 14)
        
        plt.show()
        
    
    def quick_plot_all(self): 
        """
        Call to quickly view key data. 

        Returns
        -------
        None.

        """
        self.plotIntegratedIQUV(title = True)
        self.plotSpectrumIQUV(title = True)
        self.plotCPL(title = True)
        