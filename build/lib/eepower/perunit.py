###################################################################
#   PERUNIT.PY
#
#   A library of functions, constants and more
#   that are related to Per-Unit in Electrical Engineering.
#
#   February 13, 2019
#
#   Written by Joe Stanley
#
#   Special Thanks To:
#   Paul Ortmann - Idaho Power
#
#   Included Functions
#   - Per Unit Base Creator:        pu
#   - Per Unit Base Converter:      convert
#   DOCUMENTATION IS NOT UP TO DATE!!!!!!!!!!!!!
####################################################################

# Import libraries as needed:
import numpy as np

###################################################################
#   Define per unit base creator function
#
#   Calculate per unit base given voltage (V) and power (VA)
#   Inputs: v, s, phase, z
#   When calculating 3-phase per unit bases, use 3-phase-S-base
#   and line-to-line-voltage for s and v respectively.
#
#   Returns per unit base of z if z=True, pu base of i if z=False
#   Returns as value for 3-phase by default, can also provide
#   1-phase values.
###################################################################
def zpu(S,VLL=None,VLN=None,phase=3):
    if(VLL==None and VLN==None):
        raise ValueError("ERROR: One voltage must be provided.")
    if VLL!=None:
        return(VLL**2/S)
    else:
        return((np.sqrt(3)*VLN)**2/S)


def ipu(S,VLL=None,VLN=None,phase=3):
    if(VLL==None and VLN==None):
        raise ValueError("ERROR: One voltage must be provided.")
    if VLL!=None:
        return(S/(np.sqrt(3)*VLL))
    else:
        return(S/VLN)

###################################################################
#   Define per unit converter function
#
#   Converts a [quantity] to a new per-unit base (puB_new) given an
#   old per-unit base (puB_old).
#
#   Returns per-unit value in new base.
###################################################################
def convert(quantity, puB_old, puB_new):
    pu_new = quantity*puB_old/puB_new
    return(pu_new)
    

# Define Recomposition Function
def zrecompose(z_pu,S3phs,VLL=None,VLN=None):
    """
    zrecompose Function
    
    Function to reverse per-unit conversion and return the ohmic value
    of an impedance given its per-unit parameters of R and X (as Z).
    
    Parameters
    ----------
    z_pu:       complex
                The per-unit, complex value corresponding to the
                impedance
    --- MORE ---
    """
    # Evaluate the per-unit impedance
    zbase = zpu(S3phs,VLL,VLN)
    # Evaluate the impedance
    z = z_pu * zbase
    return(z)

# Define X/R Recomposition Function
def rxrecompose(x_pu,XR,S3phs,VLL=None,VLN=None):
    """
    rxrecompose Function
    
    Function to reverse per-unit conversion and return the ohmic value
    of an impedance given its per-unit parameters of X.
    
    Parameters
    ----------
    x_pu:       float
                The per-unit, complex value corresponding to the
                impedance
    --- MORE ---
    """
    # Ensure Absolute Value
    x_pu = abs(x_pu)
    # Find R from X/R
    r_pu = x_pu/XR
    # Compose into z
    z_pu = r_pu + 1j*x_pu
    # Recompose
    z = zrecompose(z_pu,S3phs,VLL,VLN)
    return(z)
    

# END OF FILE