#
# lensdemo_funcs.py
#
# Function module for strong lensing demos
#
# Intended for use with lensdemo_script.py
#
# Copyright 2009 by Adam S. Bolton
# Creative Commons Attribution-Noncommercial-ShareAlike 3.0 license applies:
# http://creativecommons.org/licenses/by-nc-sa/3.0/
# All redistributions, modified or otherwise, must include this
# original copyright notice, licensing statement, and disclaimer.
# DISCLAIMER: ABSOLUTELY NO WARRANTY EXPRESS OR IMPLIED.
# AUTHOR ASSUMES NO LIABILITY IN CONNECTION WITH THIS COMPUTER CODE.
#

import numpy as N

def xy_rotate(x, y, xcen, ycen, phi):
    """
    NAME: xy_rotate

    PURPOSE: Transform input (x, y) coordiantes into the frame of a new
             (x, y) coordinate system that has its origin at the point
             (xcen, ycen) in the old system, and whose x-axis is rotated
             c.c.w. by phi degrees with respect to the original x axis.

    USAGE: (xnew,ynew) = xy_rotate(x, y, xcen, ycen, phi)

    ARGUMENTS:
      x, y: numpy ndarrays with (hopefully) matching sizes
            giving coordinates in the old system
      xcen: old-system x coordinate of the new origin
      ycen: old-system y coordinate of the new origin
      phi: angle c.c.w. in degrees from old x to new x axis

    RETURNS: 2-item tuple containing new x and y coordinate arrays

    WRITTEN: Adam S. Bolton, U. of Utah, 2009
    """
    phirad = N.deg2rad(phi)
    xnew = (x - xcen) * N.cos(phirad) + (y - ycen) * N.sin(phirad)
    ynew = (y - ycen) * N.cos(phirad) - (x - xcen) * N.sin(phirad)
    return (xnew,ynew)

def gauss_2d(x, y, par):
    """
    NAME: gauss_2d

    PURPOSE: Implement 2D Gaussian function

    USAGE: z = gauss_2d(x, y, par)

    ARGUMENTS:
      x, y: vecors or images of coordinates;
            should be matching numpy ndarrays
      par: vector of parameters, defined as follows:
        par[0]: amplitude
        par[1]: intermediate-axis sigma
        par[2]: x-center
        par[3]: y-center
        par[4]: axis ratio
        par[5]: c.c.w. major-axis rotation w.r.t. x-axis
        
    RETURNS: 2D Gaussian evaluated at x-y coords

    NOTE: amplitude = 1 is not normalized, but rather has max = 1

    WRITTEN: Adam S. Bolton, U. of Utah, 2009
    """
    (xnew,ynew) = xy_rotate(x, y, par[2], par[3], par[5])
    r_ell_sq = ((xnew**2)*par[4] + (ynew**2)/par[4]) / N.abs(par[1])**2
    return par[0] * N.exp(-0.5*r_ell_sq)

def sie_grad(x, y, par):
    """
    NAME: sie_grad

    PURPOSE: compute the deflection of an SIE potential

    USAGE: (xg, yg) = sie_grad(x, y, par)

    ARGUMENTS:
      x, y: vectors or images of coordinates;
            should be matching numpy ndarrays
      par: vector of parameters with 1 to 5 elements, defined as follows:
        par[0]: lens strength, or 'Einstein radius'
        par[1]: (optional) x-center (default = 0.0)
        par[2]: (optional) y-center (default = 0.0)
        par[3]: (optional) axis ratio (default=1.0)
        par[4]: (optional) major axis Position Angle
                in degrees c.c.w. of x axis. (default = 0.0)

    RETURNS: tuple (xg, yg) of gradients at the positions (x, y)

    NOTES: This routine implements an 'intermediate-axis' convention.
      Analytic forms for the SIE potential can be found in:
        Kassiola & Kovner 1993, ApJ, 417, 450
        Kormann et al. 1994, A&A, 284, 285
        Keeton & Kochanek 1998, ApJ, 495, 157
      The parameter-order convention in this routine differs from that
      of a previous IDL routine of the same name by ASB.

    WRITTEN: Adam S. Bolton, U of Utah, 2009
    """
    # Set parameters:
    b = N.abs(par[0]) # can't be negative!!!
    xzero = 0. if (len(par) < 2) else par[1]
    yzero = 0. if (len(par) < 3) else par[2]
    q = 1. if (len(par) < 4) else N.abs(par[3])
    phiq = 0. if (len(par) < 5) else par[4]
    eps = 0.001 # for sqrt(1/q - q) < eps, a limit expression is used.
    # Handle q > 1 gracefully:
    if (q > 1.):
        q = 1.0 / q
        phiq = phiq + 90.0
    # Go into shifted coordinats of the potential:
    phirad = N.deg2rad(phiq)
    xsie = (x-xzero) * N.cos(phirad) + (y-yzero) * N.sin(phirad)
    ysie = (y-yzero) * N.cos(phirad) - (x-xzero) * N.sin(phirad)
    # Compute potential gradient in the transformed system:
    r_ell = N.sqrt(q * xsie**2 + ysie**2 / q)
    qfact = N.sqrt(1./q - q)
    # (r_ell == 0) terms prevent divide-by-zero problems
    if (qfact >= eps):
        xtg = (b/qfact) * N.arctan(qfact * xsie / (r_ell + (r_ell == 0)))
        ytg = (b/qfact) * N.arctanh(qfact * ysie / (r_ell + (r_ell == 0)))
    else:
        xtg = b * xsie / (r_ell + (r_ell == 0))
        ytg = b * ysie / (r_ell + (r_ell == 0))
    # Transform back to un-rotated system:
    xg = xtg * N.cos(phirad) - ytg * N.sin(phirad)
    yg = ytg * N.cos(phirad) + xtg * N.sin(phirad)
    # Return value:
    return (xg, yg)
