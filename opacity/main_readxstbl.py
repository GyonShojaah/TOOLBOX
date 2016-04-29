#=============================================================================
# Module
#=============================================================================

from setting_readxs import *
import numpy as np
import datetime
import util_interp
import errors
#import lowres

#=============================================================================
# functions
#=============================================================================

#=============================================================================
def read_lookuptable(xsfile, wn_min, wn_max):
    """
    Extract Look-up Table
    """
    data = np.load(xsfile)

    PP_grid  = data['P']*1.0e3 # mbar => barye (x 1000)
    TT_grid     = data['T']
    WN_grid_org = data['WN']
    XS_grid_org = data['XS']

    #------------------------------------------------
    # check wavenumber range
    #------------------------------------------------
    if ( (wn_min < min(WN_grid_org)) or (wn_max > max(WN_grid_org)) ) : 
        errors.exit_msg("Wavenumbers set are out of range of look-up tables.")
    #------------------------------------------------

    indx_min = nearestindex_WN(WN_grid_org, wn_min)
    indx_max = nearestindex_WN(WN_grid_org, wn_max)

    WN_grid = WN_grid_org[indx_min:indx_max+1]
    XS_grid = XS_grid_org[indx_min:indx_max+1]

    #### TEST (to be eventually removed)
    id1, id2, id3 = np.where(XS_grid <= 0)
    for i1, i2, i3 in zip(id1, id2, id3) :
        XS_grid[i1][i2][i3] = 1.e-48
    #### TEST

    return WN_grid, TT_grid, PP_grid, XS_grid



#=============================================================================
def nearestindex_WN(WN_lattice, wn) :
    """
    Returns the wavenumber index nearest to the given wavelength
    """
    idx = np.abs(WN_lattice-wn).argmin()
    return idx



#=============================================================================
# main
#=============================================================================

if __name__ == "__main__":

    now = datetime.datetime.now()
    print now.strftime("%Y-%m-%d %H:%M:%S")

    #------------------------------------------------
    # set up look-up tables
    #------------------------------------------------
    wn_min, wn_max, wn_num = 100.0, 33333.0, 10000
    if (USER_WN_ON):
        wn_min, wn_max, wn_num = WN_MIN, WN_MAX, WN_NUM
    else:
        tmp = np.load(XSFILE)
        wn_min_tmp, wn_max_tmp, wn_num_tmp = tmp['WN'][0], tmp['WN'][-1], len(tmp['WN'])
        if (wn_min_tmp > wn_min): wn_min = wn_min_tmp
        if (wn_max_tmp < wn_max): wn_max = wn_max_tmp
        if (wn_num_tmp < wn_num): wn_num = wn_num_tmp
        WN_lookuptable, TT_lookuptable, PP_lookuptable, XS_lookuptable = read_lookuptable(XSFILE, wn_min, wn_max)

    if (CNTNM_ON):
        WN_lookuptable, TT_lookuptable, PP_lookuptable, XS_CNTNM_lookuptable = read_lookuptable(XSFILE_CNTNM, wn_min, wn_max)
    #------------------------------------------------

    #------------------------------------------------
    # set up the outspectrum grid
    #------------------------------------------------
    WN_outspectrum = np.linspace(wn_min, wn_max, wn_num)
    #------------------------------------------------

    #------------------------------------------------
    # compute radiation at each wavenumber grid
    #------------------------------------------------
    for wni in range(wn_num):

        #--------------------------------------------
        # compute radiation
        wn = WN_outspectrum[wni]
        wni_near = nearestindex_WN(WN_lookuptable, wn)
        func_XSofTP = util_interp.interp_rect_spline(TT_lookuptable, PP_lookuptable, XS_lookuptable[wni_near])
        if (CNTNM_ON):
            func_XS_CNTNMofTP = util_interp.interp_rect_spline(TT_lookuptable, PP_lookuptable, XS_CNTNM_lookuptable[wni_near])

        #--------------------------------------------
        if (CNTNM_ON):
            result = func_XSofTP(TEMP, PRES*1e3) + func_XS_CNTNMofTP(TEMP, PRES*1e3)
        else:
            result = func_XSofTP(TEMP, PRES*1e3)
        # store the result
        if wni == 0 :
            outspectrum = np.array([result])
        else :
            outspectrum = np.r_[outspectrum, result]
        #--------------------------------------------

    #------------------------------------------------

    #------------------------------------------------
    # save rawdata
    #------------------------------------------------
    #np.savez(OUTFILE_TAG+".npz", WN=WN_outspectrum, SP=outspectrum)
    data = np.dstack([WN_outspectrum, outspectrum])
    np.savetxt(OUTFILE_TAG+".txt", data[0])
    #------------------------------------------------

    now = datetime.datetime.now()
    print now.strftime("%Y-%m-%d %H:%M:%S")

