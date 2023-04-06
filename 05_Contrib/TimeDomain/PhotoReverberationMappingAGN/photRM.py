"""
Module name: photRM

Authors: Isidora Jankov & Dr. Andjelka Kovačević

Contact: ijankov@proton.me
------------------------------------------------

Description:
------------
A set of functions for:

- generating artificial AGN light curves using a method based on the Damped 
  Random Walk process where real physical AGN quantities are taken into account 
  (see Kovačević et al. 2021).
- generating artificial AGN photometric light curves for testing photometric 
  reverberation mapping methods (see Jankov et al. 2022).
- implementation of photometric reverberation mapping formalism described in 
  Chelouche & Daniel (2012)
- utility functions for preparing these light curves for pyZDCF 
  (Jankov et al., in prep), ZDCF and PLIKE programs (Alexander 1997).
  
Python functions for generating artificial light curves (lc_conti, lc_response,
lc_merged), as well as filters_viz() are adapted from original python code 
written by Dr. Andjelka Kovačević. There were several improvements introduced
in this module:
    
    - lc_conti() was rewritten to utilize numpy arrays instead of nested for 
      loops for improved efficiency;
    - added functionality to lc_conti(): option for user-specified time lag or 
      bolometric luminosity;
    - generation of the Y band (continuum + emission line) artificial light
      curve was modularized for ease of use (lc_response() and lc_merged());
    - lc_two_bands() added as a wrap-up function for generation of photoRM-ready
      pair of artificial light curves;
    - filters_viz(): added option to choose a photometric system (including
      ZTF which is not supported by the 'speclite' module by default);
    - added documentation, comments and references.
    
The rest of the functions are written to support the photoRM procedure.


Attribution
-----------
If you use these functions for scientific work leading to a publication,
please cite:

- Kovačević et al. (2021, MNRAS, 505, 5012-5028)
- Jankov et al. (2022, Astronomische Nachrichten, 343, e210090)

where the method and code for generating artificial light curves is described, 
tested, and applied.
"""

# Imports:

# Standard libs
import os

# 3rd party libs
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import interpolate
from scipy.signal import find_peaks
from speclite import filters
import astropy.units as u

# Local libs
from add_asym import add_asym


      #----------------------------------------------------------------#
      #    Functions for generating artificial AGN light curves        #
      # (adapted from original python code by Dr. Andjelka Kovačević,  #
      #              improvements by Isidora Jankov)                   #
      #----------------------------------------------------------------#

def lc_conti(T, osc=True, A=0.14, noise=0.00005, z=0, method='Kelly1', 
             lag='random', lum='random'):
    """ 
    Generate one artificial light curve using a stochastic model based on the 
    Damped random walk (DRW) proccess. Parameters describing the model are 
    characteristic amplitude ("logsig2" in code) and time scale of the 
    exponentially-decaying variability ("tau" in code), both infered from 
    physical quantities such as bolometric luminosity of the AGN. It is also 
    possible to incorporate information about the underlying periodicity by 
    setting oscillation parameter to True and providing the amplitude value. 
    For further details regarding the model see Kovačević et al. (2021) and 
    references therein.
 
    
    Parameters:
    -----------
    T: int
        Total time span of the light curve. It is recommended to generate light
        curves to be at least 10 times longer than their characteristic 
        timescale (Kozłowski 2017). 
    osc: bool, default=True
        If True, light curve simulation will take an oscillatory signal into 
        account.
    A: float, default=0.14
        Amplitude of the oscillatory signal in magnitudes (used only if 
        oscillations=True).
    noise: float, default=0.00005
        Amount of noise to include in the light curve simulation.
    z: float, default=0
        Redshift.
    method: {'Kelly1', 'Kelly2'}, default='Kelly1'
        Method for calculating DRW model parameters (Kelly et al. 2009).
    lag: int, float or 'random', default='random'
        Time-lag in days. If set to 'random', it is infered from randomly 
        generated bolometric luminosity (in range 42.2 - 45.5). You can also 
        set it to a value (int or float) to be treated as a user defined 
        time-lag in days. In this case it is used to infer bolometric luminosity.
    lum: int, float or 'random', default='random'
        The logarithm of bolometric luminosity. If set to 'random', it is drawn 
        from a uniform distribution of values in range (42.2, 45.5).You can 
        also set it to a value (int or float) to be treated as a user defined 
        bolometric luminosity. In this case, a time-lag is infered from 
        this luminosity. It is advised to use values in range (40, 47.5). More 
        extreme values may or may not work with longer light curves (>2000 days).
        
    Returns:
    --------
    tt: np.array
        Days when the light curve was sampled.
    yy: np.array
        Magnitudes of the simulated light curve.
    err: np.array
        Light curve error (Ivezić et al. 2019)
    rblr: float
        BLR radius of the simulated light curve in days (i.e., time-lag).
        
    References:
    -----------
    * Ivezić, Ž., et al. 2019, ApJ, 873, 111 
      (https://iopscience.iop.org/article/10.3847/1538-4357/ab042c)
    * Kelly, B.C., Bechtold, J., & Siemiginowska, A. 2009, ApJ, 698, 895 
      (https://iopscience.iop.org/article/10.1088/0004-637X/698/1/895)
    * Kovačević, A., et al. 2021, 505, 5012-5028
      (https://academic.oup.com/mnras/article-abstract/505/4/5012/6292266?redirectedFrom=PDF)
    * Kozłowski, S. 2017, A&A, 597, A128 
      (https://www.aanda.org/articles/aa/full_html/2017/01/aa29890-16/aa29890-16.html)
    * Shankar, F., Weinberg, D. H. & Miralda-Escude, J. 2009, ApJ, 690, 20
      (https://iopscience.iop.org/article/10.1088/0004-637X/690/1/20/pdf)
    * Woo, J. H., & Urry, C. M. 2002, ApJ, 579, 530
      (https://iopscience.iop.org/article/10.1086/342878/fulltext/)
    """
    
    # Constants
    const1 = 0.455*1.25*1e38
    const2 = np.sqrt(1e09)
    meanmag = 23.
    
    # Generating survey days 
    tt = np.arange(1, T+1)
    times = tt.shape[0]
    
    # Generating L_bol
    
    # Case 1: random values
    if (lag == 'random') and (lum == 'random'):
        loglumbol = np.random.uniform(42.2,45.5,1)
        lumbol = np.power(10,loglumbol)
    
    # Case 2: infer luminosity from user-defined lag
    elif (lag != 'random') and (lum == 'random'):
        # check user input
        try:
            lag = float(lag)
        except ValueError:
            print("Keyword argument 'lag' must be of type: int or float")
            
        lumbol = ((lag)/(10**(1.527)))**(1/0.533)*1e44
    
    # Case 3: user-defined bolometric luminosity
    elif (lag == 'random') and (lum != 'random'):
        # check user input
        try:
            lum = float(lum)
        except ValueError:
            print("Keyword argument 'lum' must be of type: int or float")
        loglumbol = lum
        lumbol = np.power(10,loglumbol)
    
    # Pathological case
    elif (lag != 'random') and (lum != 'random'):
        print('You can not specify time-lag and bolometric luminosity at the' \
              'same time! One must be set to "random"!')
    
    # Pathological case
    else:
        try:
            lag = float(lag)
        except ValueError:
            print("Keyword 'lag' has invalid value.")
        
        try:
            lum = float(lum)
        except ValueError:
            print("Keyword 'lum' has invalid value.")
            
        print("Something is wrong with keywords 'lag' or 'lum': invalid input.")
     
    print("Properties of the simulated AGN object:")
    print(39*'-')
    print("log(L) = {:.2f}".format(np.log10(float(lumbol))))

    # Calculate supermassive black hole mass using L, 
    # Eddington ratio (eq. 18, Shankar+2009) and 
    # Eddington luminosity (Woo & Urry 2002)
    msmbh=np.power((lumbol*const2/const1),2/3.)
    print("MBH = {:.2e} M_sun".format(float(msmbh)))
    
    # Calculate DRW model parameters (Kelly et al. 2009): 
    # damping time scale (tau) & amplitude of correlation decay (sig)
    if method == 'Kelly1':
        # tau = f(L,Mbh)
        tau = 80.4*np.power(lumbol/1e45,-0.42)*np.power(msmbh/1e08,1.03) # Eq 31
        tau = tau*(1+z)
        logsig2 = -3.83-0.09*np.log(lumbol/1e45)-0.25*np.log(msmbh/1e08) # Eq 30
        sig = np.sqrt(np.power(10,logsig2))/np.sqrt(1+z)
        
    elif method == 'Kelly2':
         # tau = f(L,z)
        logtau = -8.13+0.24*np.log10(lumbol)+0.34*np.log10(1+z) # Eq 22
        tau = np.power(10,logtau)*(1+z)                         # Eq 17
        logsig2 = 8-0.27*np.log10(lumbol)+0.47*np.log10(1+z)    # Eq 25
        sig = np.sqrt(np.power(10,logsig2))/np.sqrt(1+z)        # Eq 18
        
    print("tau_DRW = {:.2f} days".format(float(tau)))
    print("sigma_DRW = {:.2f} mag^2/day".format(float(sig)))
    
    # Calculate the broad line region radius (Bentz et al. 2013)
    logrblr = 1.527 + 0.533*np.log10(lumbol/1e44)
    rblr = np.power(10,logrblr)
    rblr=rblr.item()
    print("Time-lag = {:.2f} days".format(rblr))
    
    
    # Calculating light curve magnitudes using DRW model
    
    ss = np.zeros(times)
    ss[0] = meanmag # light curve is initialized
    SFCONST2 = sig*sig
    ratio = -1/tau

    for i in range(1, times):
        ss[i] = np.random.normal(ss[i-1]*np.exp(ratio) + meanmag*(1-np.exp(ratio)),
                                     np.sqrt(10*0.5*tau*SFCONST2*((1-np.exp(2*ratio)))),1)
        
    # Calculating light curve error (Ivezic et al. 2019)    
    err = lc_err(ss)
    
    # Final light curve with oscillations
    if osc == True:
        # Calculate underlying periodicity
        conver=173.145 # convert from LightDays to AU
        lightdays=10
        P = np.sqrt(((lightdays*conver)**3)/(msmbh))
        # Calculating and adding oscillatory signal
        sinus=A*np.sin(2*np.pi*tt/(P*365))
        ss = ss + sinus
        yy = np.zeros(times)
        for i in range(times):
            # Adding noise to each magnitude value
            yy[i] = ss[i] + np.random.normal(0,((noise*ss[i])),1)
    
    # Final light curve without oscillations
    elif osc == False:
        yy = np.zeros(times)
        for i in range(0,times):
            # Adding noise to each magnitude value
            yy[i] = ss[i] + np.random.normal(0,((noise*ss[i])),1)
    
    return tt, yy, err, rblr

def lc_response(tc, fxc, rblr, norm=True, plot_kernel=False):
    """
    Model the response light curve representing the isolated line emission 
    lagging behind the continuum. The response is obtained by performing 
    convolution between the continuum flux (X band*) and a Gaussian transfer
    function (Chelouche & Daniel 2012). It is assumed that the continuum 
    emission in the Y band** is the same as the X band continuum emission.
    
    *  X band - hypothetical band covering only the continuum emission.
    ** Y band - hypothetical band covering emission line superimposed on the 
                continuum.

    Parameters
    ----------
    tc : np.array
        Temporal dimension of the pure continuum ligh curve.
    fxc : np.array
        Flux dimension of the pure continuum light curve.
    rblr : int or float
        Radius of the broad line region (days).
    norm : bool, default=True
        If True, norm the input continuum flux. Otherwise, leave it as is.
    plot_kernel : bool, default=False
        If True, plot the Gaussian kernel used in convolution with the 
        continuum flux (X-band) to obtain the emission line response function.

    Returns
    -------
    response: np.array
        The emission line response flux values.
    """
    # normalize flux
    if norm:
        fxc_norm = fxc/fxc.max()
    else:
        fxc_norm = fxc
        
    # Calculate one part of the Y band: the emission line (response)
    # --> Convolution between continuum flux (X band) and a Gaussian transfer function
    
    mu = rblr
    std = np.sqrt(mu)/2
    
    transfer_l = lambda t: 23*np.exp(-np.power(t - mu, 2.) / (np.power(std, 2.)))
    
    if plot_kernel:
        # Plot the Gaussian kernel
        plt.plot(transfer_l(np.arange(150)), label=r'$\psi \ (Y_l)$')
        plt.legend(fontsize=15)
        plt.title("Gaussian kernel")
        plt.xlabel("Time (days)")
    
    y_norm = np.convolve(np.ones_like(fxc_norm), transfer_l(tc), mode='full')
    valid_indices = (y_norm > 0.)
    y_norm = y_norm[valid_indices]
    response = np.convolve(fxc_norm, transfer_l(tc), 'full')[valid_indices]/y_norm
    
    return response[:len(tc)]
    
def lc_merged(tc, fxc, response, wl=0.2, wc=0.8, norm=True):
    """
    Perform the weighted sum of the pure continuum flux and the emission line 
    response flux to obtain the Y band* light curve.
    
    * Y band - hypothetical band covering emission line superimposed on the 
               continuum.

    Parameters
    ----------
    tc : np.array
        Temporal dimension of the pure continuum ligh curve.
    fxc : np.array
        Flux dimension of the pure continuum light curve.
    response: np.array
        The emission line response flux values.
    wl : float, default=0.2
        Weight of the emission line response in the total flux.
    wc : float, default=0.8
        Weight of the continuum in the total flux.
    norm : bool, default=True
        If True, norm the input continuum flux. Otherwise, leave it as is.

    Returns
    -------
    merged : np.array
        Light curve accounting for the hypothetical emission line and its 
        surrounding continuum (Y band light curve).

    """
    
    if norm:
        fxc_norm = fxc/fxc.max() # normalize flux
    else:
        fxc_norm = fxc
    
    # 3. Final Y band light curve
    # --> Light curve covering emission line and surrounding continuum
    # --> Weighted sum between Y band continuum and Y band emission line
    
    merged = wc*fxc_norm[:len(tc)] + wl*response[:len(tc)]
    
    return merged

def lc_two_bands(T, osc=True, A=0.14, noise=0.00005, z=0, method='Kelly1',
                 lag='random', lum='random', wl=0.2, wc=0.8, plot_kernel=False):
    """
    Generate two artificial light curves: 
        
        - one in the hypothetical X band covering solely the continuum emission;
        - one in the hypothetical Y band covering some emission line and its
          surrounding continuum;
    
    The X band light curve is generated using a stochastic model based on the 
    Damped random walk (DRW) proccess (see lc_conti() docs). The Y band light
    curve is modeled using a method described in Jankov et al. (2022) which is 
    based on the photoRM formalism by Chelouche & Daniel (2012).

    Parameters
    ----------
    T: int
        Total time span of the light curve. It is recommended to generate light
        curves to be at least 10 times longer than their characteristic 
        timescale (Kozłowski 2017). 
    osc: bool, default=True
        If True, light curve simulation will take an oscillatory signal into 
        account.
    A: float, default=0.14
        Amplitude of the oscillatory signal in magnitudes (used only if 
        oscillations=True).
    noise: float, default=0.00005
        Amount of noise to include in the light curve simulation.
    z: float, default=0
        Redshift.
    method: {'Kelly1', 'Kelly2'}, default='Kelly1'
        Method for calculating DRW model parameters (Kelly et al. 2009).
    lag: int, float or 'random', default='random'
        Time-lag in days. If set to 'random', it is infered from randomly 
        generated bolometric luminosity (in range 42.2 - 45.5). You can also 
        set it to a value (int or float) to be treated as a user defined 
        time-lag in days. In this case it is used to infer bolometric luminosity.
    lum: int, float or 'random', default='random'
        The logarithm of bolometric luminosity. If set to 'random', it is drawn 
        from a uniform distribution of values in range (42.2, 45.5).You can 
        also set it to a value (int or float) to be treated as a user defined 
        bolometric luminosity. In this case, a time-lag is infered from 
        this luminosity. It is advised to use values in range (40, 47.5). More 
        extreme values may or may not work with longer light curves (>2000 days).
    wl : float, default=0.2
        Weight of the emission line response in the total flux.
    wc : float, default=0.8
        Weight of the continuum in the total flux.
    plot_kernel : bool, default=False
        If True, plot the Gaussian kernel used in convolution with the 
        continuum flux (X-band) to obtain the emission line response function.

    Returns
    -------
    x_band: pd.DataFrame
        Time, flux and error of the continuum light curve packed into a 
        DataFrame.
    y_band: pd.DataFrame
        Time, flux and error of the continuum + em. line light curve packed 
        into a DataFrame.
    line_response: np.array
        Time, flux and error of the emission line response function packed into
        a DataFrame.
    """
    
    # Generate the pure continuum light curve (X band)
    tc, fxc, err_c, rblr = lc_conti(T,osc=osc, A=A, noise=noise, z=z, 
                                    method=method, lag=lag, lum=lum)
    
    # Calculate the response function of the emission line
    response = lc_response(tc, fxc, rblr, plot_kernel=plot_kernel)
    
    # Generate the light curve covering both continuum and the lagged emission
    # line (Y band).
    merged = lc_merged(tc, fxc, response, wl=wl, wc=wc)
    
    fxc_norm = fxc/fxc.max() # norm the continuum flux
    
    # Calculate light curve errors
    err_merg = np.sqrt(np.max(lc_err(merged)/merged)*merged)
    err_c = np.sqrt(np.max(lc_err(wc*fxc_norm)/wc*fxc_norm)*wc*fxc_norm)
    err_em = np.sqrt(np.max(lc_err(wl*response)/wl*response)*wl*response)
           
    # Create DataFrames for continuum LC (X band), continuum + line LC (Y band)
    # and for emission line response
    
    data_cont = {'t':tc,
                 'flux':fxc_norm,
                 'err':err_c}
    
    data_response = {'t':tc,
                     'flux':wl*response,
                     'err':err_em}
    
    data_merged = {'t':tc,
                   'flux':merged,
                   'err':err_merg}
    
    x_band = pd.DataFrame(data_cont)
    y_band = pd.DataFrame(data_merged)
    line_response = pd.DataFrame(data_response)
        
    return x_band, y_band, line_response


def lc_err(lc):
    """
    Calculates photometric uncertainty by implementing photometric error model
    of LSST (Ivezic et al. 2019).
    """
    gamma=0.039
    m5=24.7
    x=np.zeros(lc.shape)
    x=np.power(10, 0.4*(lc-m5))

    err = (0.005*0.005) + (0.04-gamma)*x + gamma*x*x
    return err
    
def lc_plot(tt, yy, T):
    """
    Simple plotting function.
    
    Parameters:
    -----------
    tt: np.array
        Days when the light curve was sampled.
    yy: np.array
        Light curve magnitudes.
    T: int
        Total time span of the light curve. 
    """
    
    fig = plt.figure(figsize=(15,5))
    
    ax = fig.add_subplot(111)
    ax.plot(tt, yy, 'ko', markersize = 1, label='1 day cadence')

    custom_xlim = (0, T)
    custom_ylim = (yy.min()-0.1, yy.max()+0.1)
    ax.set_xlabel('t [days]', fontsize = 18, labelpad=10)
    ax.set_ylabel('magnitude', fontsize = 18, labelpad=10)
    ax.tick_params(direction='in', pad = 5, labelsize=13)
    plt.setp(ax, xlim=custom_xlim, ylim=custom_ylim)
    ax.legend(fontsize=15)
    ax.grid(True)


       #----------------------------------------------------------#
       # Functions for processing the results from pyZDCF or ZDCF #
       #   programs and time-lag estimation (by Isidora Jankov)   #
       #----------------------------------------------------------#

def add_inverted_acf(acf):
    new_acf = pd.DataFrame(columns=['tau', '-sig(tau)', '+sig(tau)', 'dcf', '-err(dcf)', '+err(dcf)', '#bin'])
    new_acf['dcf'] = np.append(np.flip(acf['dcf']),acf['dcf'])
    new_acf['tau'] = np.append(np.flip(acf['tau']*(-1)), acf['tau'])
    new_acf['+sig(tau)'] = np.append(np.flip(acf['+sig(tau)']*(-1)), acf['+sig(tau)'])
    new_acf['-sig(tau)'] = np.append(np.flip(acf['-sig(tau)']*(-1)), acf['-sig(tau)'])
    new_acf['+err(dcf)'] = np.append(np.flip(acf['+err(dcf)']),acf['+err(dcf)'])
    new_acf['-err(dcf)'] = np.append(np.flip(acf['-err(dcf)']),acf['-err(dcf)'])
    new_acf['#bin'] = np.append(np.flip(acf['#bin']),acf['#bin'])
    return new_acf

def interp(a,b):
    # a: df with common grid of tau values
    # b: df which we want to force upon that grid
    
    new_b = pd.DataFrame(columns=['tau', '-sig(tau)', '+sig(tau)', 'dcf', '-err(dcf)', '+err(dcf)', '#bin'])

    for y in ['dcf','-sig(tau)','+sig(tau)','-err(dcf)','+err(dcf)','#bin']:
        if y == 'dcf':
            kind = 'quadratic'
        else:
            kind = 'nearest'
        f = interpolate.interp1d(b['tau'], b[y], kind=kind, fill_value="extrapolate")
        new_b[y] = f(a['tau'])
    
    new_b['tau'] = a['tau']
        
    return new_b

def delta_ccf(acf, ccf):
    """
    Subtract ACF from CCF and use error propagation (Laursen et al. 2019) to 
    estimate asymmetric errors in the resulting function.
    """
    
    delta = pd.DataFrame(columns=['tau', '-sig(tau)', '+sig(tau)', 'dcf', '-err(dcf)', '+err(dcf)', '#bin'])
    
    delta['tau'] = acf['tau']
    delta['#bin'] = acf['#bin']

    for i in ccf.index:
        x0 = [ccf.loc[i,'dcf'],-acf.loc[i,'dcf']]
        s1 = [ccf.loc[i,'-err(dcf)'], acf.loc[i,'-err(dcf)']]
        s2 = [ccf.loc[i,'+err(dcf)'], acf.loc[i,'+err(dcf)']]
        delta_val, err1, err2 = add_asym(x0,s1,s2,order=1)
        delta.loc[i,'dcf'] = delta_val
        delta.loc[i,'-err(dcf)'] = err1
        delta.loc[i,'+err(dcf)'] = err2
        
    for i in ccf.index:
        x0 = [ccf.loc[i,'tau'],-acf.loc[i,'tau']]
        s1 = [ccf.loc[i,'-sig(tau)'], acf.loc[i,'-sig(tau)']]
        s2 = [ccf.loc[i,'+sig(tau)'], acf.loc[i,'+sig(tau)']]
        delta_val, err1, err2 = add_asym(x0,s1,s2,order=1)
        delta.loc[i,'-sig(tau)'] = err1
        delta.loc[i,'+sig(tau)'] = err2
    
    # Change data types (must be done for PLIKE to work)
    delta['-sig(tau)'] = delta['-sig(tau)'].astype('int64')
    delta['+sig(tau)'] = delta['+sig(tau)'].astype('int64')
    delta['dcf'] = delta['dcf'].astype('float64')
    delta['-err(dcf)'] = delta['-err(dcf)'].astype('float64')
    delta['+err(dcf)'] = delta['+err(dcf)'].astype('float64')
    delta['#bin'] = delta['#bin'].astype('int64')
        
    return delta


def plot_ccf_acf(delta, ccf, acf, locator=10, save=False, peak=False, tau=0, err=(0,0), lims_x=(-20,80), lims_y=(-0.25,1.25)):
    """
    Plot CCF, ACF and their difference. Optionally, you can add peak location 
    (tau) and associated errors.
    
    Parameters:
    -----------
    delta: pd.DataFrame
    	Table returned by delta_ccf() function. Contains information about 
        CCF-ACF.
    ccf: pd.DataFrame
    	Table containing information about CCF.
    acf: pd.DataFrame
    	Table containing information about ACF.
    locator: int, default=10
        Parameter indicating a step for plotting x-axis ticks.
    save: bool, default=False
        If True, the resulting plot will be saved as ccf-acf.pdf
    peak: bool, default=False
        If True, plot will show peak location (but you also need to specify 
        'tau' keyword argument)
    tau: float, default=0
    	Time-lag value.
    err: tuple, default=(0,0)
    	Time-lag (-) and (+) error. Enter (-) error with '-' sign.
    lims_x: tuple
    	Lower and upper limit on x-axis of both plots.
    lims_y: tuple
    	Lower and upper limit on y-axis of both plots.
    """
    
    x1, x2 = lims_x
    y1, y2 = lims_y
    err_low, err_high = err
    
    fig = plt.figure(figsize=(9,6))
    plt.subplots_adjust(left=None, bottom=None, right=None, top=0.95, wspace=0,hspace=0)
    
    ax1 = fig.add_subplot(211)
    ax1.plot(acf['tau'], acf['dcf'], 'o-r', markersize=3, linewidth=1, label='ACF')
    ax1.plot(ccf['tau'], ccf['dcf'], 'o--k', markersize=3, linewidth=1, label='CCF')
    ax1.set_xlabel("Time")
    ax1.grid(which='major', axis='x', linestyle='--')
    ax1.xaxis.set_major_locator(plt.MultipleLocator(locator))
    ax1.yaxis.set_major_locator(plt.MultipleLocator(0.2))
    plt.setp(ax1.get_xticklabels(), visible=False)
    ax1.legend(loc='lower left', fontsize=13)
    #ax1.set_xlim(lim1,lim2)
    ax1.set_xlim(x1,x2)
    ax1.set_ylim(y1,y2)
    
    ax2 = fig.add_subplot(212)
    
    lower_error_x =  delta['-sig(tau)'] 
    upper_error_x =  delta['+sig(tau)']
    asymmetric_error_x = np.array(list(zip(lower_error_x, upper_error_x))).T
    
    lower_error_y =  delta['-err(dcf)']
    upper_error_y =  delta['+err(dcf)']
    asymmetric_error_y = np.array(list(zip(lower_error_y, upper_error_y))).T
    
    ax2.errorbar(delta['tau'], delta['dcf'], xerr=asymmetric_error_x, fmt='o-k',
                 markersize=2, linewidth=0.5, label= 'CCF - ACF', ecolor='gray', capsize=2)
    
    if peak == True:
        nearest_idx = np.abs(delta['tau'] - tau).argmin()
        ax2.vlines(tau, ymin=-1, ymax= delta['dcf'][nearest_idx],linestyles='dashed', colors='royalblue')
        ax2.axvspan(tau+err_low, tau+err_high, alpha=0.2)
    
    ax2.set_xlabel("Time (days)", fontsize=17, labelpad=8)
    ax2.grid(which='major', axis='x', linestyle='--')
    ax2.xaxis.set_major_locator(plt.MultipleLocator(locator))
    ax2.yaxis.set_major_locator(plt.MultipleLocator(0.1))
    ax2.legend(fontsize=13, loc='lower left')
    ax2.set_ylim(-0.25,0.25)
    ax2.set_xlim(x1,x2)
    if peak==True:
        ax2.text(0.72,0.05, r'$\tau$ = {:.1f} ({:.1f},+{:.1f}) d'.format(tau,err_low,err_high), transform=ax2.transAxes, size=12.5)
    
    fig.text(0.03, 0.5, "Correlation (arbit. units)", va='center', rotation='vertical',fontsize=18)
    if save==True:
        plt.savefig('ccf-acf.pdf',dpi=800)
    plt.show()
    
    
def peak_finder(y, x, interval, first=False):
    """
    Find highest positive peak in a user defined interval. Useful for finding
    coresponding time-lag in CCF. By setting first=True, you can also specify 
    that you want to select first peak found in a given interval, instead of 
    the highest one.
    """
    lim1, lim2 = interval
    condition = (x>=lim1) & (x<lim2)
    x = x[condition]
    y = y[condition]
    
    indexes, _ = find_peaks(y)

    peaks_x = []
    peaks_y = []
    for idx in indexes:
        if x[idx] > 0:
            peaks_x.append(x[idx])
            peaks_y.append(y[idx])

    peaks_x = np.asarray(peaks_x)
    peaks_y = np.asarray(peaks_y)
    tau = peaks_x[peaks_y.argmax()]
    
    print('Peak candidates (x-axis vals): ', peaks_x)
    
    # get the first peak only in a given interval
    if first:
        tau = peaks_x[0]
        peaks_x = peaks_x[0]
        peaks_y = peaks_y[0]
        print('First peak: ', tau)
        
    # get max. peak from peak candidates in a given interval
    else:
        tau = peaks_x[peaks_y.argmax()]
        print('Max. peak: ', tau)
    
    return tau, peaks_x, peaks_y



# Additional functions

def filters_viz(z=0, phot_sys='LSST', save=False, comp_spec_path='./data/'):
    """
    The function returns a plot of LSST/SDSS/ZTF broadband filter response
    curves with composite quasar spectrum at a given redshift. It also works
    if you use filter identifiers from speclite 
    (https://speclite.readthedocs.io/en/latest/filters.html).
    """
    
    # Read composite quasar spectrum (Vanden Berk et al. 2001)
    data = pd.read_csv(comp_spec_path+'comp_spec.txt', skiprows=23, header=None, sep=" ", skipinitialspace=True)
    data.columns = ['Wave', 'FluxD', 'e_FluxD']
    
    # Load LSST filters
    if phot_sys=='LSST':
        filt = filters.load_filters('lsst2016-*')
        
    # Load SDSS filters
    elif phot_sys=='SDSS':
        filt = filters.load_filters('sdss2010-*')
    
    # Load ZTF filters
    elif phot_sys=='ZTF':
        directory_name = './ZTF_data/ZTF_filters'
        fg_name = os.path.join(directory_name, 'ztf-g.ecsv')
        fr_name = os.path.join(directory_name, 'ztf-r.ecsv')
        fi_name = os.path.join(directory_name, 'ztf-i.ecsv')
        filt = filters.load_filters(fg_name,fr_name,fi_name)
    else:
        filt = filters.load_filters(phot_sys)
    
    
    # Plotting
    t1 = data['Wave'].values*u.Angstrom
    t1obs=(z+1)*t1
    filters.plot_filters(filt, wavelength_limits=(700, 11000), legend_loc='upper right')
    plt.plot(t1obs,(0.2*data['FluxD']/data['FluxD'].max()),label=r'$z={}$'.format(z),c='navy',linewidth=1.5)
    plt.legend(loc='upper left', fontsize=14)
    fig = plt.gcf()
    fig.set_size_inches(11,7)
    plt.xlabel(r"$\mathrm{Wavelength \ (\AA)}$",size=17, labelpad=7)
    plt.ylabel("Filter response",size=18, labelpad=9)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    
    CIV_wave = (z+1)*1549
    Ha_wave = (z+1)*6563
    Hb_wave = (z+1)*4861
    MgII_wave = (z+1)*2798
    CIII_wave = (z+1)*1908
    #Lya_wave = (z+1)*1216
    
    if CIV_wave < 11000:
        plt.annotate('C IV', xy =(CIV_wave, 0.1),
                    xytext =(CIV_wave, 0.1), size=13)
    if Ha_wave < 11000:
        plt.annotate(r'H$\alpha$', xy =(Ha_wave, 0.03),
                 xytext =(Ha_wave, 0.04), size=13)
    if Hb_wave < 11000:
        plt.annotate(r'H$\beta$', xy =(Hb_wave, 0.03),
                     xytext =(Hb_wave, 0.03), size=13)
    if MgII_wave < 11000:
        plt.annotate('Mg II', xy =(MgII_wave, 0.04),
                    xytext =(MgII_wave, 0.04), size=13)
    if CIII_wave < 11000:
        plt.annotate('C III', xy =(CIII_wave, 0.06),
                    xytext =(CIII_wave, 0.06), size=13)
    plt.grid(False)
        
    #if Lya_wave < 11000:
    #    plt.annotate(r'Ly$\alpha$ ', xy =(Lya_wave, 0.25),
    #                xytext =(Lya_wave, 0.21), size=13)
    
    if save == True:
        plt.savefig('filters.png',dpi=300)
    
