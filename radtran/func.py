import numpy as np
from scipy import constants
import cgs

#=============================================================================
def planck(TTarg, wnarg):
    """
    Returns planck function [erg/sec/cm2]
    """
    exponent = cgs.HH*cgs.CC*wnarg/( cgs.KK*TTarg )
    planckfunction = 2.0*cgs.HH*cgs.CC**2*wnarg**3/(np.exp(exponent)-1)
    return planckfunction


