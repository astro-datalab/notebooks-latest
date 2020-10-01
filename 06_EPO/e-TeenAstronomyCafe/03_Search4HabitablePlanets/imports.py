#%matplotlib inline
import matplotlib.pyplot as plt
from IPython.display import display
from ipywidgets import interactive

# Import standard Python science packages
import astropy.units as u
import numpy as np

# Import the transit modeling package
import batman

repoURL = 'https://raw.githubusercontent.com/DavidVargasMora/TACTests/master/'

def example_transit_plot(Time,time_1,flux_1):
    # Makes a new figure
    fig,ax = plt.subplots(2,1,figsize=(14,5))

    r_star,r_planet = 1.0, 0.1
    x = 24*(time_1-time_1[0])
    
    # Compute the planet's current position from the Time
    tmid = x[int(len(x)/2)]
    tdur = np.ptp(x[flux_1<1.])
    position = (Time-tmid)/(tdur/2.)
    
    star = plt.Circle((0,0), r_star,color='gold')
    ax[0].add_artist(star)
    planet = plt.Circle((position, 0.35), r_planet, color='black')
    ax[0].add_artist(planet)
    ax[0].set_xlim(-4*r_star, 4*r_star)
    ax[0].set_ylim(-1.2*r_star, 1.2*r_star)
    ax[0].axis('off')
    ax[0].set_aspect('equal')
    ax[0].set_title("Move the slider to watch the planet transit the star")

    x = 24*(time_1-time_1[0])
    idx = np.argmin(np.abs(x-Time))
    #idx = -1
    ax[1].plot(x[:idx], flux_1[:idx],lw=5,c='gold')
    ax[1].set_title("Observe how the star's brightness changes during the transit:")
    ax[1].scatter(Time,flux_1[np.argmin(np.abs(x-Time))],marker='x',s=200,c='black',zorder=10)
    ax[1].axvline(Time,linestyle='dashed',c='black',lw=2)
    ax[1].set_xlabel("Time (hours)")
    ax[1].set_ylabel("Star brightness")
    ax[1].set_xlim([-0.1,2.8])
    ax[1].set_ylim([0.988,1.002])
    ax[1].set_yticklabels(['{:.1%}'.format(xx) for xx in ax[1].get_yticks()])
    
    return

def transit_model_1_plot(period):
    """
    Function for interactive transit modeling
    """
    # Parameters of planet model
    params = batman.TransitParams()
    params.t0 = 0.                                     # time of inferior conjunction
    params.per = period                                # orbital period (days)
    params.rp = ((1.0*u.Rjup)/(1.203*u.Rsun)).si.value # planet radius (in units of stellar radii)
    params.a = 15.                                     # semi-major axis (in units of stellar radii)
    params.inc = 87.                                   # orbital inclination (in degrees)
    params.ecc = 0.                                    # eccentricity
    params.w = 90.                                     # longitude of periastron (in degrees)
    params.u = [0.1, 0.1]                              # limb darkening coefficients [u1, u2]
    params.limb_dark = "quadratic"                     # limb darkening model

    # Defining the timeframe and calculating flux
    time = np.linspace(-1, 1, 10000)*6        # Creates an array of 10000 points 
                                              # over a roughly 12-day period
    model = batman.TransitModel(params, time)
    flux_model = model.light_curve(params)

    # Plotting results
    plt.figure(figsize=(14,5))
    plt.plot(time-time[0], flux_model, 'r-')
    plt.title("Figure 3: Multiple Transits of a Hot Jupiter with a {:0.1f}-day Period".format(params.per))
    plt.xlabel("Time (days)")
    plt.ylabel("Star brightness")
    
    return
    #plt.show()

def transit_model_2_plot(radius):
    # Parameters of planet model
    rp = ((radius*u.Rjup)/(1.203*u.Rsun)).si.value
    params = batman.TransitParams()
    params.t0 = 0.                                 # time of inferior conjunction
    params.per = 3.                                # orbital period (days)
    params.rp =  rp                                # planet radius (in units of stellar radii)
    params.a = 15.                                 # semi-major axis (in units of stellar radii)
    params.inc = 87.                               # orbital inclination (in degrees)
    params.ecc = 0.                                # eccentricity
    params.w = 90.                                 # longitude of periastron (in degrees)
    params.u = [0.1, 0.1]                          # limb darkening coefficients [u1, u2]
    params.limb_dark = "quadratic"                 # limb darkening model

    # Defining the timeframe and calculating flux
    time = np.linspace(-0.025, 0.025, 100)*params.per # Creates an array of 100 points 
                                                      # near the transit time
    model = batman.TransitModel(params, time)            # Initializes model
    flux_model = model.light_curve(params)            # Calculates light curve

    # Loading observed data
    time_obs, flux_obs = np.loadtxt(repoURL+'03_Search4HabitablePlanets/Data/HD209458b.dat', unpack=True)
    
    # Convert times to hours (easier to read)
    time = 24*(time-time[0])
    time_obs = 24*(time_obs-time_obs[0])
    
    # Plotting results
    plt.figure(figsize=(14,5))
    plt.errorbar(time_obs, flux_obs, 0.001, marker='o', ls='', color='b', label='Observed')
    plt.plot(time, flux_model, 'r-', label='Model')
    plt.title("Figure 5: HD 209458 b")
    plt.legend(loc=4)
    plt.xlabel("Time (hours)")
    plt.ylabel("Star brightness")
    plt.gca().set_yticklabels(['{:.1%}'.format(xx) for xx in plt.gca().get_yticks()])
    
    return

def transit_model_3_plot(Rp,Rs):
    # Parameters of planet model
    params = batman.TransitParams()
    params.t0 = 0.                                   # time of inferior conjunction
    params.per = 3.                                  # orbital period (days)
    params.rp = ((Rp*u.earthRad)/(Rs*u.Rsun)).si.value # planet radius (in units of stellar radii)
    params.a = 15.                              # semi-major axis (in units of stellar radii)
    params.inc = 87.                                 # orbital inclination (in degrees)
    params.ecc = 0.                                  # eccentricity
    params.w = 90.                                   # longitude of periastron (in degrees)
    params.u = [0.1, 0.1]                            # limb darkening coefficients [u1, u2]
    params.limb_dark = "quadratic"                   # limb darkening model

    # Defining the timeframe and calculating flux
    time = np.linspace(-0.025, 0.025, 100)*params.per # Creates an array of 100 points 
                                                      # near the transit time
    model = batman.TransitModel(params, time)         # Initializes model
    flux_model = model.light_curve(params)            # Calculates light curve
    noise = np.random.normal(0, 0.001, flux_model.shape)
    flux_obs = flux_model + noise

    # Convert time to hours
    time = 24*(time-time[0])
        
    # Plotting results
    plt.figure(figsize=(14,5))
    plt.errorbar(time, flux_obs, 0.001, marker='o', ls='', color='b', label='Observed')
    plt.plot(time, flux_model, 'r-', label='Model')
    plt.title("Figure 6: Transit of a {:0.1f}-Earth-radius planet in front of a {:0.1f}-Sun-radius star".format(Rp, Rs))
    plt.legend(loc=4)
    plt.xlabel("Time (days)")
    plt.ylabel("Star brightness")
    plt.gca().set_yticklabels(['{:.1%}'.format(xx) for xx in plt.gca().get_yticks()])
    return

def planet_size_plot(m_planet,m_star,r_planet,r_star):
    fs = 22
    fig, ax = plt.subplots(figsize=(8, 8))
    star = plt.Circle((0,0), r_star,color='gold')
    ax.add_artist(star)
    ax.plot([0, r_star], [0, 0], 'k-')
    ax.text(r_star/2., 0, '$R_{s}$', va='top', ha='center',fontsize=fs)
    planet = plt.Circle((-0.35, 0.35), r_planet, color='black')
    ax.add_artist(planet)
    ax.plot([-0.35, -0.35+r_planet], [0.35, 0.35], 'k-',color='white')
    ax.text(-0.35+r_planet/2., 0.35, '$R_{p}$', va='top', ha='center',color='white',fontsize=fs)
    
    r_max = r_star if r_star>r_planet else 1.5*r_planet
    ax.set_xlim(-r_max, r_max)
    ax.set_ylim(-r_max, r_max)
    plt.axis('off')
    return

def plot_data(time_obs,flux_obs):
    # Convert time from days to hours (easier to read)
    time_obs = 24*(time_obs-time_obs[0])

    # Plotting results
    plt.figure(figsize=(14,5))
    plt.errorbar(time_obs, flux_obs, 0.001, marker='o', ls='', color='b', label='Observed')
    plt.title("Figure 4: Observations of HD 209458 b")
    plt.legend(loc=4)
    plt.xlabel("Time (hours)")
    plt.ylabel("Star brightness")
    plt.gca().set_yticklabels(['{:.1%}'.format(xx) for xx in plt.gca().get_yticks()])

def plot_multiple_transits(time_2,flux_2):
    plt.figure(figsize=(14,5))
    plt.plot(time_2-time_2[0], flux_2)
    plt.title("Figure 2: Multiple Transits of a Hot Jupiter")
    plt.xlabel("Time (days)")
    plt.ylabel("Star brightness")
    plt.gca().set_yticklabels(['{:.1%}'.format(xx) for xx in plt.gca().get_yticks()])
    plt.show()
    