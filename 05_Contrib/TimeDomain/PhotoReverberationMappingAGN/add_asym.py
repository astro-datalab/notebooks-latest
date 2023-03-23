import numpy as np
from numpy import pi,sqrt,exp

sqrt2pi = sqrt(2.*pi)

def add_asym(x0s,siglos,sighis,order=2,nmax=100,tol=1e-6,ohwell=False):
    """
    Purpose
    -------
    Add set of two or more numbers with asymmetric uncertainties.
    
    Attribution
    -----------
    If you use this function for scientific work leading to a publications,
    please cite
      Laursen et al. (2019, A&A, 627, 84),
    where the code is described, tested, and applied.

    Description
    -----------
    Each number is characterized by three numbers, viz. its central value x0,
    and its standard deviations siglo and sighi toward negative and positive
    values, respectively.
    
    That is, the numbers would be written as
      X = x0_{-siglo_x}^{+sighi_x}
      Y = y0_{-siglo_y}^{+sighi_y}
      Z = z0_{-siglo_z}^{+sighi_z}
      etc.

    The standard approach for this problem is to add the central values
    normally, and to add the lower and upper uncertainties separately in
    quadrature, i.e.  
      X + Y = {x0+y0}_{-siglo_tot}^{+sighi_tot},
    where
      siglo_tot^2 = siglo_x^2 + siglo_y^2;
      sighi_tot^2 = sighi_x^2 + sighi_y^2.
    This method has no statistical foundation, and is wrong.

    There is no one right method to solve this problem (as infinitely many PDFs
    could be described by the same three numbers), but one method that seems to
    give acceptable results for many different distributions is described in
    Barlow, R. (2003), http://arxiv.org/abs/physics/0306138v1.
    In short, the distributions are transformed --- linearly or quadratically
    --- to "proper" Gaussians, for which the normal, linear addition is valid.
    The total mean, variance, and skewness can thus be found, and these values
    are then used to transform back to get the three numbers for the final
    distribution.
    Whereas the "forward" transformation can be done analytically, the
    "backward" transformation must be done numerically by iteration.

    Parameters
    ----------
    x0s:    array of central values.
    siglos: array of lower uncertainties.
    sighis: array of upper uncertainties.
    order:  order of the transformation (1 for linear, 2 for quadratic, and
            0 for the standard, but wrong, method of adding upper and
            lower uncertainties separately in quadrature.
    nmax:   Maximum number of iterations (usually <10 are needed).
    tol:    Tolerance for the accepted result.
    ohwell: If True, suppress warning when using order=0.

    Returns
    -------
    numpy array containing the three numbers x0_tot, siglo_tot, and sighi_tot
    that best describe the PDF of the sum.

    Example
    -------
    Do the following addition, using a linear transformation:
    5_{-2}^{+1}  +  3_{-3}^{+1}  +  4_{-3}^{+2}.
    >>> x0 = [5,3,4]
    >>> s1 = [2,3,3]
    >>> s2 = [1,1,2]
    >>> add_asym(x0,s1,s2,order=1)
    array([ 10.92030132,   4.23748786,   2.94389109])
    
    That is, the sum would be written as 10.9_{-4.2}^{+2.9}.

    Do the same, using a quadratic transformation:
    >>> add_asym(x0,s1,s2,order=2)
    array([ 10.65168907,   4.47929965,   3.1759215 ])

    Notice that the result hinges somewhat on the chosen order, but that both
    find a lower central value than the "standard", but wrong, result of
    12_{-4.69}^{+2.45}. The reason is that in this case, all three addends have
    PDFs that are skewed toward lower values. Note also that both methods give
    more symmetric errors, in accord with the Central Limit Theorem. In
    contrast, the asymmetry of addends whose skewness are of the same sign will
    never decrease.
    """

    assert order in [0,1,2], "\nParameter `order` must be 1 or 2 (or 0, but... well, stick to 1 or 2)."
    x0s    = np.asarray(x0s)
    siglos = np.asarray(siglos)
    sighis = np.asarray(sighis)

    # ----- Standard, but wrong, method of adding lower and upper errors ------
    if order==0:                                  #   separately in quadrature |
        if not ohwell:                            #                            |
            print('This is WRONG! Mark my words, it is WROOOOOONG!!!')#        |
        x0    = np.sum(x0s)                       #                            |
        siglo = sqrt(np.sum(siglos**2))           #                            |
        sighi = sqrt(np.sum(sighis**2))           #____________________________|

    else:
      # ------------------   First, calculate total moments   -----------------
      mu,V,gamma = 0.,0.,0.                               #                    |
      if order==1:                                        #Cumul. mean,var,skew|
          sig    = (sighis+siglos) / 2.                   #"The mean" }_eq. 1  |
          alpha  = (sighis-siglos) / 2.                   #"The diff."}        |
          mu     = np.sum(x0s + (sighis-siglos) / sqrt2pi)#Biased mean         |
          V      = np.sum(sig**2 + alpha**2 * (1-2/pi))   #Variance            |
          t1     = sighis**3 - siglos**3                  #\                   |
          t2     = (sighis-siglos) * (sighis**2+siglos**2)# > 3 terms for gamma|
          t3     = (sighis-siglos)**3                     #/                   |
          gamma  = np.sum(2.*t1-1.5*t2+t3/pi) / sqrt2pi   #Skewness            |
      else:                                               #Same for order = 2  |
          sig    = (sighis+siglos) / 2.                   #"The mean" }_eq. 1  |
          alpha  = (sighis-siglos) / 2.                   #"The diff."}        |
          mu     = np.sum(x0s + alpha)                    #Biased mean         |
          V      = np.sum(sig**2 + 2*alpha**2)            #Variance            |
          gamma  = np.sum(6*sig**2 * alpha + 8*alpha**3)  #Skewness____________|

      # -----------------   Check if iteration is necessary  ------------------
      check = abs(gamma/mu**3) if mu != 0. else gamma #Avoid NaN for mu = 0    |
      if check < 1e-12:                           #If errors are too close to  |
          x0    = mu                              # being symmetric, the method|
          siglo = sqrt(V)                         # below gives an exception,  |
          sighi = siglo                           # so calculate normally______|

      else:
      # -----------  Use moments to find resulting distribution  --------------
          n = 0                                  #Counter for # of iterations  | 
          if order==1:                           #                             |
              D = 0.                             #Iterator; eq. to sighi-siglo |
              while True:                        #                             |
                  n   += 1                       #Update counter               |
                  S    = 2*V + D**2/pi           #S is equal to siglo^2+sighi^2|
                  Dold = D                       #                             |
                  D    = 2/(3.*S) * (sqrt2pi*gamma - D**3 * (1/pi-1)) #Iterate!|
                  if Dold != 0.:                 #                             |
                      if abs(D/Dold-1.) < tol: break #Accept if close enough   |
                  assert n<nmax, 'Too many iterations'#If it takes too long, fuck it
              S     = 2*V + D**2/pi              #Final update of S            |
              x0    = mu - D/sqrt2pi             #Biased mean                  |
              alpha = D / 2.                     #"The difference" }_ eq. 1    |
              sig   = sqrt(V - alpha**2*(1-2/pi))#"The mean"       }           |
              siglo = sig - alpha                #Lower error                  |
              sighi = sig + alpha                #Upper error                  |
                                                 #                             |
          elif order==2:                         #                             |
              alpha = 0.                         #Iterator; equal to (hi-lo)/2 |
              while True:                        #                             |
                  n    += 1                      #Update counter               |
                  aold  = alpha                  #                             |
                  alpha = gamma / (6*V - 4*alpha**2) #Iterate!                 |
                  if aold != 0.:                 #                             |
                      if abs(alpha/aold-1.) < tol: break#Accept if close enough|
                  assert n<nmax, 'Too many iterations'#If it takes too long, fuck it
              sig   = sqrt(V - 2*alpha**2)       #"The mean" (eq 1)            |
              x0    = mu  - alpha                #\                            |
              siglo = sig - alpha                # >Final values               |
              sighi = sig + alpha                #/____________________________|

    return np.array([x0, siglo, sighi])
#------------------------------------------------------------------------------
