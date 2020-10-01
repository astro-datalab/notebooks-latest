import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from astropy.io import fits, ascii
from astropy.table import Table
import batman
from ipywidgets import interactive, fixed
from IPython.display import display
import pdb
import os

plt.rcParams.update({'font.size': 18,'figure.figsize': (10,8)})

t = np.linspace(-1,1,100) ## time

repoURL = 'https://raw.githubusercontent.com/noaodatalab/notebooks-latest/master/06_EPO/e-TeenAstronomyCafe/'

def initial_lightcurve(radius=0.1):
    params = batman.TransitParams()
    params.t0 = 0.                       #time of inferior conjunction
    params.per = 48.                     #orbital period
    params.rp = radius                   #planet radius (in units of stellar radii)
    params.a = 15.                       #semi-major axis (in units of stellar radii)
    params.inc = 87.                     #orbital inclination (in degrees)
    params.ecc = 0.                      #eccentricity
    params.w = 90.                       #longitude of periastron (in degrees)
    params.u = [0.2, 0.0]                #limb darkening coefficients [u1, u2]
    params.limb_dark = "quadratic"       #limb darkening model
    
    return params


def plot_lc(radius=0.1,color='blue'):
    params = initial_lightcurve(radius=radius)
    t = np.linspace(-1, 1, 100)
    m = batman.TransitModel(params, t)    #initializes model
    flux = m.light_curve(params)          #calculates light curve
    plt.plot(t, flux * 100., color=color)
    plt.xlabel("Time from central transit (hours)")
    plt.ylabel("Relative brightness (%)")
    plt.ylim(98.5,100.2)
    return radius

def plot_initial_lc():
    plot_lc(radius=0.1)


def lc_jup_radii(radius=1.0):
    plot_lc(radius/10.)

def plot_interactive_rad():
    lc_rad = interactive(lc_jup_radii, radius=(0.0,1.2,0.05),color=fixed('blue'))
    
    display(lc_rad)
    return lc_rad


slopeStart = 0.0
slopeEnd = 0.3

class spectral_lc:
    def __init__(self):
        """
        Class for spectroscopic lightcuves
        """
        self.wavelengths =  np.array([  0.64   ,  0.61    , 0.57  ,   0.53   , 0.47     , 0.41 ])
        self.waveRange = np.max(self.wavelengths) - np.min(self.wavelengths)
        self.calc_radii()
        #self.radius_array = np.array([  0.08  ,  0.085  , 0.090  , 0.095   ,0.100   ,  0.105])
        self.colors_array = np.array([  'red' ,'orange','yellow' ,'green',  'blue',  'violet'])
        #self.limb_dark_arr = [   0.1 , 0.15   ,   0.2,   0.25     , 0.3     , 0.35]
    
    def calc_radii(self,Thickness=0.3):
       
        self.radius_array = 0.08 - 0.1 * Thickness * (self.wavelengths - np.max(self.wavelengths)) / self.waveRange
    
    def plot_lc_multicolor_loop(self,Thickness=slopeStart):
        self.calc_radii(Thickness=Thickness)
        for one_radius,one_color in zip(self.radius_array,self.colors_array):
            plot_lc(one_radius,color=one_color)
    
    def plot_lc_multicolor(self,Thickness=-0.3):
        lc_multi = interactive(self.plot_lc_multicolor_loop, Thickness=(slopeStart,slopeEnd,0.02))
        display(lc_multi)
    
    def spectrum_plot(self,Thickness=slopeStart):
        self.calc_radii(Thickness=Thickness)
        plt.plot(self.wavelengths,self.radius_array * 10.,color='black')
        for oneWave, oneRad, oneCol in zip(self.wavelengths,self.radius_array * 10.,self.colors_array):
            plt.plot([oneWave],[oneRad],'s',color=oneCol,markersize=12)
        
        plt.xlabel('Wavelength (microns)')
        plt.ylabel('Size (Earth Radii)')
        plt.ylim(0.7,1.15)
        
    def spectrum_plot_i(self):
        spec_p = interactive(self.spectrum_plot,Thickness=(slopeStart,slopeEnd,0.02))
        display(spec_p)
        
    def visualize_colors(self,Thickness=slopeStart):
        self.calc_radii(Thickness=Thickness)
        
        fig, ax = plt.subplots(figsize=(8,8))
        ax.set_aspect('equal')
        ax.set_xlim(-1.3,1.3)
        ax.set_ylim(-1.3,1.3)
        for one_radius,one_color in zip(self.radius_array,self.colors_array):
            circlePatch = plt.Circle((0., 0.), one_radius * 10.,linewidth=2,
                                     edgecolor=one_color,facecolor='none')
            ax.add_artist(circlePatch)

        for ind,angle in zip([0,-1],[0,1,0.7]):
            ax.plot([0,self.radius_array[ind] * np.cos(angle) * 10.],
                    [0,self.radius_array[ind] * np.sin(angle) * 10.],color=self.colors_array[ind])
            ax.text(0.2,np.sin(angle)* 0.5,"{:.2f} Earth Radii".format(self.radius_array[ind] * 10.))
        ax.set_xlabel('X Size (Earth Radii)')
        ax.set_ylabel('Y Size (Earth Radii)')
    
    def visualize_colors_i(self):
        size_p = interactive(self.visualize_colors,Thickness=(slopeStart,slopeEnd,0.02))
        display(size_p)

convertDict = {'H2O':'H$_2$O','CH4':'CH$_4$','CO2':'CO$_2$','Cloudy':'Cloudy'}

def show_example_spectra(atmospheres=['H2O','CH4','CO2','No Atmosphere']):
    dat = Table.read(repoURL+'09_Exoplanet_Spectra/Data/opacity_breakdown_gto_f_hd189733b.fits')
    dat['No Atmosphere'] = 0.019
    
    if len(atmospheres) > 1:
        fig, ax2D = plt.subplots(2,2,figsize=(12,12))
        axArr = np.ravel(ax2D)
    else:
        fig, oneAx = plt.subplots()
        axArr = [oneAx]
    
    
    
    for ind,atmosphere in enumerate(atmospheres):
        ax = axArr[ind]
        if atmosphere in list(convertDict.keys()):
            thisLabel = convertDict[atmosphere]
        else:
            thisLabel = atmosphere
    
        rad = np.sqrt(dat[atmosphere]) * 10.
        rad = (rad - np.mean(rad)) * 60. + np.mean(rad) ## exaggerate to see better
        ax.plot(dat['Wave'],rad,label=thisLabel)
        ax.set_xlim(2.5,5)
        ax.set_ylim(1.0,2.0)
        ax.legend()
        ax.set_xlabel('Wavelength (microns)')
        ax.set_ylabel('Size (Earth Radii)')



lc_offset = 0.015

class atmospheric_lc():
    def __init__(self,mysteryNum=1):
        
        datName = repoURL+'09_Exoplanet_Spectra/Data/mystery_lc_{}.fits'.format(mysteryNum)
        HDUList = fits.open(datName)
        
        #Remove path validation as the file will be stored in github.
        #if os.path.exists(datName):
        #    HDUList = fits.open(datName)
        #else:
        #    raise Exception("Mystery {} not found".format(mysterNum))
        
        self.wavelengths = HDUList['WAVE'].data
        self.lcData = HDUList['FLUX'].data
        self.radii = np.ones_like(self.wavelengths) * 1.5
        self.colors_array = cm.rainbow(np.linspace(0, 1, len(self.wavelengths)))
        self.orig_time = HDUList['TIME'].data
        self.time = np.linspace(-1,1,100) ## time
        
        HDUList.close()
        
    def plot_lc_engine(self):
        
        for ind,oneRad in enumerate(self.radii):
            flux_data = self.lcData[ind,:]
            
            params = initial_lightcurve(radius=oneRad/10.)
            m = batman.TransitModel(params, self.time)    #initializes model
            model_flux = m.light_curve(params)          #calculates light curve
            #pdb.set_trace()
            plt.plot(self.time,(flux_data - lc_offset * ind) * 100.,'o',
                     color=self.colors_array[ind])
            plt.plot(self.time,(model_flux - lc_offset * ind) * 100.,
                     color=self.colors_array[ind])
        plt.xlabel('Time from central transit (hours)')
        plt.ylabel('Relative Brightness (%)')
    
    def plot_lc(self):
        self.plot_lc_engine()
    
    def plot_lc_rad_by_rad(self,rad0,rad1,rad2,rad3,rad4,rad5,rad6,rad7,rad8,rad9):
        self.radii = [rad0,rad1,rad2,rad3,rad4,rad5,rad6,rad7,rad8,rad9]
        self.plot_lc_engine()
    
    def plot_lc_i(self):
        rs, re, ri = 1.0,2.0,0.05
        lc_multi = interactive(self.plot_lc_rad_by_rad,
                               rad0=(rs,re,ri),
                               rad1=(rs,re,ri),
                               rad2=(rs,re,ri),
                               rad3=(rs,re,ri),
                               rad4=(rs,re,ri),
                               rad5=(rs,re,ri),
                               rad6=(rs,re,ri),
                               rad7=(rs,re,ri),
                               rad8=(rs,re,ri),
                               rad9=(rs,re,ri))
        display(lc_multi)
    
    def plot_spectrum(self):
        plt.plot(self.wavelengths,self.radii)
        plt.xlabel('Time from Central Transit (hours)')
        plt.ylabel('Size (Earth Radii)')
        plt.ylim(1.0, 2.0)
## Now put in a background star
# fig, ax = plt.subplots(figsize=(8,8))
# ax.set_aspect('equal')
# ax.set_xlim(-1.3,1.3)
# ax.set_ylim(-1.3,1.3)
#
# circlePatch = plt.Circle((0., 0.), 10.,
#                          edgecolor=one_color,facecolor='yellow')
# ax.add_artist(circlePatch)
# circlePatch = plt.Circle((0., 0.), np.min(radius_array) * 10.,
#                          edgecolor='none',facecolor='black')
# ax.add_artist(circlePatch)
#
# sortInd = np.argsort(radius_array)
# for ind in sortInd:
#     thisRadius = radius_array[ind]
#     for max_radius,one_color in zip(radius_array,colors_array):
#
#         if thisRadius >= max_radius:
#             circlePatch = plt.Circle((0., 0.), max_radius * 10.,linewidth=3,
#                          edgecolor=one_color,facecolor='none',alpha=0.2)
#             ax.add_artist(circlePatch)
#
# ax.set_xlabel('X Size in Earths')
# ax.set_ylabel('Y Size in Earths')
