#=============================================================================
# Module
#=============================================================================
import numpy  as np
from   scipy  import interpolate
import errors
import cgs
import util_interp


#=============================================================================
def calc_nXSofZ_molabs(layer_z, grid_wn, tuple_func_atmprof, molname, xsfile_tag, cntnm_on=True, xsfile_tag_cntnm=None):

    func_TofZ, func_PofZ, func_MUofZ, dict_func_NofZ = tuple_func_atmprof
    layer_P  = func_PofZ(layer_z)
    layer_T  = func_TofZ(layer_z)
    layer_n0 = layer_P/(cgs.RR*layer_T)*cgs.NA
    layer_n  = dict_func_NofZ[molname](layer_z)*1e-6*layer_n0

    #------------------------------------------------
    # set up lookup tables
    #------------------------------------------------
    xsfile = xsfile_tag + molname + ".txt"
    print "Reading ", xsfile, "..."
    WN_lookuptable, XS_lookuptable = read_lookuptable(molname, xsfile, (grid_wn[0], grid_wn[-1]))

    #------------------------------------------------
    # read lookup table
    #------------------------------------------------

    interpfunc = util_interp.interp_1d( WN_lookuptable[::-1], XS_lookuptable[::-1] )
    print "WN_lookuptable[::-1]", WN_lookuptable[::-1]
    print "grid_wn", grid_wn
    grid_xs    = interpfunc( grid_wn )

    mesh_n, mesh_xs = np.meshgrid( layer_n, grid_xs )
    mesh_nXS  = mesh_n * mesh_xs

    # returned numpy.array
    # WN_NUM x Z_NUM
    # [[ wn1 x z1, wn1 x z2, .., wn1 x zN ], 
    #  [ wn2 x z1, wn2 x z2, .., wn2 x zN ],... 
    return mesh_nXS 


#=============================================================================
def calc_nXSofZ_Rayleigh(layer_z, grid_wn, tuple_func_atmprof, polarization):
    """
    Rayleigh scattering
    """
    func_TofZ, func_PofZ, func_MUofZ, dict_func_NofZ = tuple_func_atmprof
    layer_P = func_PofZ(layer_z)
    layer_T = func_TofZ(layer_z)
    layer_n0 = layer_P/(cgs.RR*layer_T)*cgs.NA
#    print layer_n0

    grid_wl = 1.0/grid_wn
    grid_XS = 128.0*np.pi**5/(3.0*grid_wl**4)*polarization**2

    mesh_n0, mesh_XS = np.meshgrid(layer_n0, grid_XS)
#    print "mesh_n0", mesh_n0
#    print "mesh_XS", mesh_XS
    mesh_nXS = mesh_n0*mesh_XS
    return mesh_nXS # WN_NUM x Z_NUM



#=============================================================================
def read_lookuptable(molname, xsfile, wn_limit):
    """
    Extract Look-up Table
    """

    data = np.loadtxt(xsfile).T

    WN_grid = data[0]
    tau_grid = data[1]
    tau_grid[np.where(tau_grid <= 0.)] = 0.

    rho    = { 'CH4': 6.6050e-7, 'isoprene': 2.805e-6 }
    weight = { 'CH4': 15.0,      'isoprene': 68.0 }
    factor = 1./19.94/(rho[molname]*1e-3/weight[molname]*6.022e23)

    XS_grid = - 1.0 * np.log( 1. - tau_grid ) * factor
    XS_grid[np.where(XS_grid <= 1e-28)] = 1e-28

    return WN_grid, XS_grid




