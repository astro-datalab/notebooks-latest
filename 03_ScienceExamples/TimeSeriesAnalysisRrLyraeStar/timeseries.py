__author__ = 'Robert Nikutta <robert.nikutta@noirlab.edu>, Astro Data Lab Team <datalab@noirlab.edu>'
__version__ = '20240427' # yyymmdd

import numpy as np
from astropy import timeseries
import astropy.units as u
import pylab as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, mark_inset
from matplotlib.ticker import MaxNLocator


def plot_scatter(x,y,ax=None):
    
    if ax is None:
        fig, ax = plt.subplots(figsize=(5,5))
    
    # apply units
    x = (x*u.deg).to('arcmin')
    y = (y*u.deg).to('arcsec')
    
    ax.plot(x-np.median(x),y-np.median(y),'bo',ms=3,alpha=0.6)
    ax.set_xlabel('RA offset (%s)' % x.unit)
    ax.set_ylabel('DEC offset (%s)' % y.unit)
    ax.axvline(0,c='0.5',lw=1)
    ax.axhline(0,c='0.5',lw=1)
    

# define a function to select only measurements in one band
def get_data(df,band='g'):
    sel = (df['filter'] == band)    
    t = df['mjd'][sel].values
    y = df['cmag'][sel].values
    dy = df['cerr'][sel].values
    return t,y,dy # return time, magnitudes in one band, uncertainties

# a reusable function to plot the lightcurve
def plot_raw_lightcurve(t,y,dy,title='',ax=None,lperc=13,rperc=99):
    
    if ax is None:
        fig, ax = plt.subplots()
    
    jd0 = t.min() # modified Julian date offset
    t = t-jd0 # first date of observations = 0
    
    axins = inset_axes(ax, 3.5, 1.3, loc=9) # make inset axes object
    
    left = np.percentile(t,lperc)
    right = np.percentile(t,rperc)
    ax.errorbar(t,y,yerr=dy,marker='.',ms=8,ls='none',color='g',lw=1,alpha=0.5,label='')
    axins.errorbar(t,y,yerr=dy,marker='.',ms=8,ls='none',color='g',lw=1,alpha=0.5,label='')
    axins.set_xlim(0.9999*left,1.0001*right)
    axins.xaxis.set_major_locator(MaxNLocator(4))

    # Main panel chores
    ax.set_xlabel('modified Julian date - %g (days)' % jd0)
    ax.set_ylabel('magnitude')
    ax.invert_yaxis()
    ax.set_title(title)
    #ax.legend(loc='lower left',frameon=True,ncol=1,markerscale=1.5)

    # draw bbox around inset; connect with parent axes
    mark_inset(ax, axins, loc1=1, loc2=4, fc="none", ec="0.5",lw=1,alpha=0.7)


def get_ls_periodogram(t,y,min_freq=1./1.,max_freq=1./0.1):
    
    """Compute Lomb-Scargle periodogram.
    
    Parameters
    ----------
    t : array
        Observation time array (e.g. MJD), ordered in ascending order.
    
    y : array
        Magnitude measurements at times ``t``.
        
    min_freq, max_freq : float or None
        The period finder can be guided by providing the min and max frequency
        in the ``y`` signal, in units 1/t. 
          min_freq = 1/longest expected period (in days)
          max_freq = 1/shortest expected perdiod (in days)
        The defaults are typical for RR Lyrae variability (RR Lyrae usually
        have a period of a fraction of one day).
        
    Returns
    -------
    period : array
        Phased period of the time-variable signal (fraction of the phase).
        
    power : array
        The periodogramm power as function if ``period``.
        
    """
    
    # Use astropy's LombScargle class
    ls = timeseries.LombScargle(t, y)

    # Compute the periodogram
    #   We guide the algorithm a bit:
    #     min_freq = 1/longest expected period (in days)
    #     max_freq = 1/shortest expected perdiod (in days)
    #   RR Lyrae usually have a period of a fraction of one day
    frequency, power = ls.autopower(minimum_frequency=min_freq,maximum_frequency=max_freq)
    period = 1./frequency # period is the inverse of frequency
    
    return period, power


def get_best_period(period,power):
    
    """Return the period with highest power."""
    
    return period[np.argmax(power)]


def plot_periodogram(period,power,best_period=None,title='',ax=None):

    """Plot a periodogram.
    
    Parameters
    ----------
    
    period, power : array
        The period and power 1-d arrays as returned by :func:`get_ls_periodogram()`
    
    best_period : float or None
        If float, the value of this ``best_period`` will be marked in the plot.
    
    title : str
        Title of the figure. Default: ''.
    
    ax : instance or None
        If instance of axis class, will plot to that object. If None, will generate a new figure and axis object.
    """
    
    if ax is None:
        fig, ax = plt.subplots()
        
    ax.plot(period,power,lw=0.1)
    ax.set_xlabel('period (days)')
    ax.set_ylabel('relative power')
    ax.set_title(title)
    
    if best_period is not None:
        ax.axvline(best_period,color='r');
        ax.text(0.03,0.93,'period = %.3f days' % best_period,transform=ax.transAxes,color='r')


def get_folded_phase(t,best_period):
    
    """Fold the observation times with the best period of the variable signal."""
    
    # light curve over period, take the remainder (i.e. the "phase" of one period)
    phase = (t / best_period) % 1
    
    return phase


def plot_folded_lightcurve(t,y,best_period,dy=None,ax=None):
    
    """Plot folded lightcurve.
    
    Parameters
    ----------
    
    t, y : array
        Time and magnitude 1-d arrays
        
    best_period : float
        True period of the signal.
        
    dy : array or None
        If array, the values are the uncertainies on ``y``, and the plot will show errorbars.
        If None, the plot will have no errorbars.
        
    ax : instance or None
        If instance of axis class, will plot to that object. If None, will generate a new figure and axis object.
    """

    phase = get_folded_phase(t,best_period)
    
    if ax is None:
        fig, ax = plt.subplots()
        
    marker = '.'
    ms = 10
    lw = 1
    color = 'g'
    alpha = 0.6
    
    
    if dy is not None:
        ax.errorbar(phase,y,yerr=dy,marker=marker,ms=ms,ls='none',lw=lw,color=color,alpha=alpha)
    else:
        ax.plot(phase,y,marker=marker,ms=ms,ls='none',lw=lw,color=color,alpha=alpha)
            
    ax.invert_yaxis()
    ax.set_xlabel('phase (days)')
    ax.set_ylabel('magnitude')


def do_everything(ra,dec,t,y,dy):
    # make fig & axes
    fig, axes = plt.subplots(2,2,figsize=(14,10))
    
    # plot data in sky
    plot_scatter(ra,dec,ax=axes[0,0])
    
    # plot lightcurve data
    plot_raw_lightcurve(t,y,dy,ax=axes[0,1])
    
    # plot periodogram
    period, power = get_ls_periodogram(t,y)
    best_period = get_best_period(period,power)
    plot_periodogram(period,power,best_period,ax=axes[1,0])
    
    # plot folded lightcurve
    plot_folded_lightcurve(t,y,best_period,dy=dy,ax=axes[1,1])


