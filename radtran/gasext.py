#=============================================================================
# Module
#=============================================================================
import sys, os
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/../common')

import netCDF4
import numpy  as np
from   scipy  import interpolate



import io_nc 
import errors
import cgs

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
    xsfile = xsfile_tag + molname + ".nc"
    WN_lookuptable, TT_lookuptable, PP_lookuptable, XS_lookuptable = read_lookuptable(xsfile, (grid_wn[0], grid_wn[-1]))

    #------------------------------------------------
    # check range of lookup table
    #------------------------------------------------
    if ( layer_P[-1] < PP_lookuptable[0] ):
        errors.exit_msg("Pressure grids in FILE_ATM smaller than the range of lookuptable.")
    if ( layer_P[0] > PP_lookuptable[-1] ):
        errors.exit_msg("Pressure grids in FILE_ATM is larger than the range of lookuptable.")

    #------------------------------------------------
    # CONTINUUM? (currently available only for CH4)
    #------------------------------------------------
    if (cntnm_on and (molname=="H2O")): 
        print "     including continuum"
        xsfile = xsfile_tag_cntnm + molname + ".npz"
        WN_lookuptable, TT_lookuptable, PP_lookuptable, XS_cntnm = read_lookuptable(xsfile, (grid_wn[0], grid_wn[-1]))
        XS_lookuptable = XS_lookuptable + XS_cntnm


    #------------------------------------------------
    # read lookup table
    #------------------------------------------------
    dict_griddata_logXSofWNTP = {}
    m_Tgrid, m_WNgrid, m_logPgrid = np.meshgrid(TT_lookuptable, WN_lookuptable, np.log(PP_lookuptable))
    dict_griddata_logXSofWNTP['coords'] = np.dstack([m_WNgrid.flatten(), m_Tgrid.flatten(), m_logPgrid.flatten()])[0]
    dict_griddata_logXSofWNTP[molname] = np.log(XS_lookuptable.flatten())

    mesh_n   = np.tile(layer_n, (len(grid_wn),1))
    mesh_logP, mesh_wn = np.meshgrid(np.log(layer_P), grid_wn)
    mesh_T, mesh_wn = np.meshgrid(layer_T, grid_wn)
    flat_points = np.dstack([mesh_wn.flatten(), mesh_T.flatten(), mesh_logP.flatten()])
    flat_logXS  = interpolate.griddata(dict_griddata_logXSofWNTP['coords'], dict_griddata_logXSofWNTP[molname], flat_points, method='nearest')[0]
    mesh_logXS = flat_logXS.reshape(len(grid_wn),len(layer_z))
    mesh_nXS = mesh_n*np.exp(mesh_logXS)

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

    grid_wl = 1.0/grid_wn
    grid_XS = 128.0*np.pi**5/(3.0*grid_wl**4)*polarization**2

    mesh_n0, mesh_XS = np.meshgrid(layer_n0, grid_XS)
    mesh_nXS = mesh_n0*mesh_XS
    return mesh_nXS # WN_NUM x Z_NUM



#=============================================================================
def read_lookuptable(xsfile, wn_limit):
    """
    Extract Look-up Table
    """
    WN_grid_org, PP_grid_org, TT_grid, XS_grid_org = io_nc.read_xstbl( xsfile )

#    if PP_grid_org.units=='mbar'
    PP_grid = PP_grid_org * 1.e3 # mbar => barye (x 1000)

    #------------------------------------------------
    # check wavenumber range
    #------------------------------------------------
    if ( (wn_limit[0] < min(WN_grid_org)) or (wn_limit[1] > max(WN_grid_org)) ) : 
        errors.exit_msg("Wavenumbers set are out of range of look-up tables.")
    #------------------------------------------------

    indx_min = nearestindex_WN(WN_grid_org, wn_limit[0])
    indx_max = nearestindex_WN(WN_grid_org, wn_limit[1])

    WN_grid = WN_grid_org[indx_min:indx_max+1]
    del WN_grid_org
    XS_grid = XS_grid_org[indx_min:indx_max+1]
    del XS_grid_org

    #### TEST (to be eventually removed)
    id1, id2, id3 = np.where(XS_grid <= 0)
    for i1, i2, i3 in zip(id1, id2, id3) :
        XS_grid[i1][i2][i3] = 1.e-100
    #### TEST

    return WN_grid, TT_grid, PP_grid, XS_grid



#=============================================================================
def nearestindex_WN(WN_lattice, wn) :
    """
    Returns the wavenumber index nearest to the given wavelength
    """
    idx = np.abs(WN_lattice-wn).argmin()
    return idx

