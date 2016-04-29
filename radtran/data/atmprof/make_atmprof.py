TEMP  = 300 # K
PSURF = 1.e2 # mbar
Z_NUM = 100
Z_MAX = 5000.0  # km
MU_ATM = 2.0 # g

import numpy as np
from scipy import constants
from scipy.optimize import fsolve
from input import *
import cgs


if __name__ == "__main__":

    layer_z = np.linspace(0, )
    SH = cgs.RR*TEMP/(MU_ATM*cgs.GG) # approximation

    layer_P = PSURF*np.exp(-layer_z/SH)
    layer_rho = 
    print "# isothermal atmospheric profile"
    print "# molecules: "
    print "# z(km)   p(mb)        T(K)	MU"
    for zi in xrange(Z_NUM):
        print layer_z[zi], layer_P[zi], TEMP, MU_ATM, 
    return TT, PP, rho
