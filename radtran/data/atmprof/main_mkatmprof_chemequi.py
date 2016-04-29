import numpy as np
from scipy import constants
from scipy.optimize import fsolve
import cgs


TEMP  = 500 # K
PSURF = 1.e4 # mbar
Z_NUM = 100
Z_MAX = 6000.0  # km
MU_ATM = 2.0 # g
molname = ["H2O"]
MIXRATIO = [98.0]
GRAV = cgs.GG*1.2


if __name__ == "__main__":

    layer_z = np.linspace(0, Z_MAX, Z_NUM)
    SH = cgs.RR*TEMP/(MU_ATM*GRAV) # approximation
#    print "SH", SH
    layer_P = PSURF*np.exp(-layer_z*1.0e5/SH)
    layer_n = layer_P*1.0e3/(cgs.RR*TEMP)*constants.N_A
    print "# isothermal atmospheric profile"
    print "# molecules: ", 
    for ii in xrange(len(molname)):
        print molname[ii],
    print ""
    print "# z(km)   p(mb)        T(K)	MU"
#    for zi in reversed(xrange(Z_NUM)):
    for zi in xrange(Z_NUM):
        print layer_z[zi], layer_P[zi], TEMP, MU_ATM, MIXRATIO
