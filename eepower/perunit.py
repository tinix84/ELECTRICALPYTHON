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
def pu(S3phs,VLL=None,VLN=None,phase=3,z=True):
    if(VLL==None and VLN==None):
        raise ValueError("ERROR: One voltage must be provided.")
	if z:
        if VLL!=None:
            return(VLL**2/S3phs)
        else:
            return((np.sqrt(3)*VLN)**2/S3phs)
	else:
        if VLL!=None:
            return(np.sqrt(3)*VLL/S3phs)
        else:
            return(VLN/S3phs)

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