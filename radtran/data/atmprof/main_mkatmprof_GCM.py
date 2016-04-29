import numpy as np
import pandas as pd
from scipy import constants
from scipy.optimize import fsolve
import cgs

MU_ATM = 28.8 # g
molname = ["CO2", "CH4"]
MIXRATIO = [3.300E+02, 1.700E+00]

#=============================================================================
def CCpressure(TT):
    """
    Approximate Clausius-Clapeyron equation. 
    Ref) A Short course in Cloud Physics eq(2.17)
    Latent heat is assumed to be constant.
    """
#    TTc = (TT-273.0)
#    return 611.2*np.exp(17.67*TTc/(TTc+243.5))
    return 2.53e11*np.exp(-5.42e3/TT)*10.0 # Pa => barye


def CCpressure_wagner(TT):

    AA = -7.764
    BB = 1.45838
    CC = -2.7758
    DD = -1.2330
    Pc = 2.212e7
    Tc = 647.3
    x = 1.0 - (TT/Tc)
    return Pc * np.exp(( AA * x + BB * x**1.5 + CC * x**3 + DD * x**6 )/(1-x)) * 10.0

#=============================================================================
if __name__ == "__main__":

    layer = pd.read_table("prof_GCM_S0.5_lon-92")

    layer_Pv = CCpressure(layer['temp'])
    layer_H2Oppm = layer_Pv/(layer['pres']*1e3)*(layer['rh']*1.0e-2)*1e6

    print "# molecules: H2O  CO2  CH4"
    print "#     z(km)      p(mb)        T(K)	MU    h2o    co2     ch4"
    for ii in xrange(len(layer['pres'])):
        print layer['alt'][ii], layer['pres'][ii], layer['temp'][ii], MU_ATM, layer_H2Oppm[ii], 
        for jj in xrange(len(molname)):
            print MIXRATIO[jj],
        print ""


