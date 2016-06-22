#=============================================================================
# Module
#=============================================================================

import numpy as np
from scipy import integrate
from scipy  import constants
import datetime
import util_interp
import errors
import cgs
import func

#=============================================================================
# functions
#=============================================================================

#=============================================================================
def run(layer_z, grid_wn, mesh_nXS):

    print mesh_nXS.shape
    layer_tau = np.zeros_like(layer_z)
    for wni in xrange(len(mesh_nXS)):
        print "-----------------------------------"
        print " wavelength : ", 1.0e4/grid_wn[wni]
        print "-----------------------------------"
        print " # altitude [km], optical depth from TOA"
        print "-----------------------------------"
        for zi in xrange(len(layer_z)):
            layer_z_upper  = layer_z[zi:]
            layer_nXS_upper = mesh_nXS[wni,zi:]
            layer_tau[zi] = integrate.trapz(layer_nXS_upper, x=layer_z_upper)
            print layer_z[zi]*1e-5, layer_tau[zi]

#    #------------------------------------------------
#    # save rawdata
#    #------------------------------------------------
#    data = np.dstack([grid_wn, grid_result])
#    np.savetxt(outdir+OUTFILE_TAG+".txt", data[0])
#    #------------------------------------------------
#




