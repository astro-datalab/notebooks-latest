__version__ = '20221121'
__author__ = 'Giada Pastorelli <gpastorelli.astro@gmail.com>'

# imports
import pylab as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

#Reusable function to plot the Hess diagrams 
def plot_2D_hess(H, xed, yed, 
                 ax = None, cmap = None, norm = None):
    '''
    Plot 2D Hess diagrams
    
    Parameters
    ----------
    H: ndarray
        2D Hess diagram    
    xed: array
        Edges along x axis
    yed: array
        Edges along y axis
    ax: Axes object
        If None get current axes to make the plot
        otherwise make plot in the specified ax
    cmap: str or Colormap
        If not specified, default is viridis
    norm: str or Normalize
        Normalization method to scale data.
        If given can be an instance of Normalize 
        or scale name, i.e. "linear", "log", "symlog" 

    Returns
    -------
    cc: AxesImage
        Axes Image from imshow
    '''

    ax = ax or plt.gca()
    
    if cmap is None:
        cmap ='viridis'
        
    cc = ax.imshow( H.T, 
                    cmap=cmap,
                    norm=norm,
                    interpolation='None',
                    origin='lower', aspect='auto',  
                    extent=[xed[0], xed[-1], yed[0], yed[-1]]
                    )
    return cc    

# Function to make a customizable colorbar for the Hess diagrams 
def make_cbar(ax, cc, label, boundaries=None, extend ='both'):
    '''
    Use a divider to make it easy to plot a colorbar

    Parameters
    ----------
    ax: Axes object
        
    cc: AxesImage

    label: str
        colorbar label
    boundaries: None or a sequence
        Optional, only used in particular cases to 
        specify the boundaries of the colorbar
    extend: str
        'neither', 'both', 'min', 'max'
         Make pointed end(s) for out-of-range values (unless 'neither')

    Returns
    -------
    cb: Colorbar
        colorbar
    '''
    
    divider = make_axes_locatable(ax)
    cax = divider.append_axes('right', size='5%', pad=0.05)
            
    cb = plt.colorbar(cc, cax=cax, orientation='vertical', extend=extend, 
                            boundaries= boundaries, 
                            extendfrac = 0.1,
                            label =  label )
    return cb