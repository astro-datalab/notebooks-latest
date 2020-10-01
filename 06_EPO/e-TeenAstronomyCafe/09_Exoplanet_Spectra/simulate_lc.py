import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits, ascii
from astropy.table import Table
import batman
from ipywidgets import interactive, fixed
from IPython.display import display
import pdb
from scipy.interpolate import interp1d
import exo_engine

nTime = 100
nWave = 10

def make_lc(showPlot=False):
    np.random.seed(0)

    dat = Table.read('opacity_breakdown_gto_f_hd189733b.fits')
    dat['No Atmosphere'] = 0.019
    orig_wave = dat['Wave']
    
    wave_array = np.linspace(2.5,5,nWave)
    fluxArray = np.zeros([nWave,nTime])
    t = np.linspace(-1, 1, nTime)
    
    for modInd,atmosphere in enumerate(['CO2','No Atmosphere','CH4','H2O']):
        
        orig_Rad = np.sqrt(dat[atmosphere]) * 10
        orig_Rad = (orig_Rad - np.mean(orig_Rad)) * 60. + np.mean(orig_Rad) ## exaggerate to see better
        orig_Rad = orig_Rad / 10. ## divide by star radius
        f = interp1d(orig_wave,orig_Rad)
        rad_array = f(wave_array)
    
        for ind,oneRad in enumerate(rad_array):
            rad = np.sqrt(dat[atmosphere]) * 10.
        
            params = exo_engine.initial_lightcurve(radius=oneRad)
            m = batman.TransitModel(params, t)    #initializes model
            flux = m.light_curve(params)          #calculates light curve
            flux_data = flux + np.random.randn(nTime) * 0.0007
            
            fluxArray[ind,:] = flux_data
            
            if showPlot == True:
                plt.plot(t,flux_data - 0.01 * ind)
        if showPlot == True:
            plt.show()
        
        hduFlux = fits.PrimaryHDU(fluxArray)
        hduFlux.name = 'FLUX'
        
        hduTime = fits.ImageHDU(t)
        hduTime.name = 'TIME'
        
        hduWave = fits.ImageHDU(wave_array)
        hduWave.name = 'WAVE'
        
        hduList = fits.HDUList([hduFlux,hduTime,hduWave])
        hduList.writeto('mystery_lc_{}.fits'.format(modInd+1),overwrite=True)
        
    
if __name__ == '__main__':
    make_lc()
    