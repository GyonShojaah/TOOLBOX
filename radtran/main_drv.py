#=============================================================================
# Module
#=============================================================================

from setting import *
from copy import deepcopy
import numpy as np
import datetime
import util_interp
import errors
import cgs
import cProfile
import gasext
import cld_simple

#------------------------------------------------
if (OBSMODE=="trans"):
    outdir = "out/transmission/"
    import obsmode_trans
#------------------------------------------------
elif (OBSMODE=="disk"):
    outdir = "out/emission/"
    import obsmode_disk
#------------------------------------------------
elif (OBSMODE=="nadir"):
    outdir = "out/nadir/"
    import obsmode_nadir
#------------------------------------------------
else:
    errors.exit_msg("INVALID OBSMODE")
#------------------------------------------------


#=============================================================================
# functions
#=============================================================================

def set_atmprof(infile):
    """
    Set up profile of amtmosphere
    """
    list_mol, dict_atmprof = read_atmprofile(infile)
    func_TofZ  = util_interp.interp_1d_boundary(dict_atmprof['z'], dict_atmprof['T'],  order=1, logx=False, logy=False)
    func_PofZ  = util_interp.interp_1d_boundary(dict_atmprof['z'], dict_atmprof['P'],  order=1, logx=False)
    func_MUofZ = util_interp.interp_1d_boundary(dict_atmprof['z'], dict_atmprof['MU'], order=1, logx=False, logy=False)
    dict_func_NofZ = {}
    for ii in xrange(len(list_mol)):
        dict_func_NofZ[list_mol[ii]] = util_interp.interp_1d_boundary(dict_atmprof['z'], dict_atmprof[list_mol[ii]], order=3, logx=False, logy=False)
    tuple_func_atmprof = (func_TofZ, func_PofZ, func_MUofZ, dict_func_NofZ)
    return dict_atmprof, list_mol, tuple_func_atmprof


#=============================================================================
def read_atmprofile(infile, key=0) :
    """
    Read profile of amtmosphere
    """
    with open(infile,'r') as f_linefile :
        for line in f_linefile :
            elements = line.rsplit()
            if elements[1] == "molecules:" :
                list_mol = elements[2:len(elements)]
            if elements[0] != '#' :
                break
    names = ['z', 'P', 'T', 'MU']
    formats = ['f8', 'f8', 'f8', 'f8']
    for ii in xrange(len(list_mol)) :
        names.append(list_mol[ii])
        formats.append('f8')

    dict_atmprof = np.loadtxt(infile, 
                      dtype={'names':tuple(names),
                             'formats':tuple(formats)}, 
                      comments='#')


    #------------------------------------------------
    # mbar => barye (x 1000)
    #------------------------------------------------
    dict_atmprof['P'] = dict_atmprof['P']*1e3

    #------------------------------------------------
    # km => cm
    #------------------------------------------------
    dict_atmprof['z'] = dict_atmprof['z']*1e5
    #------------------------------------------------

    #------------------------------------------------
    # check the order of altitude in FILE_ATM
    #------------------------------------------------
    reverse = 0
    z_old = dict_atmprof['z'][0]
    for zi in xrange(1,len(dict_atmprof['z'])):
        if (dict_atmprof['z'][zi] < z_old):
            reverse = reverse + 1
    if (reverse == len(dict_atmprof['z'])-1):
        dict_atmprof['z'] = dict_atmprof['z'][::-1]
        dict_atmprof['P'] = dict_atmprof['P'][::-1]
        dict_atmprof['T'] = dict_atmprof['T'][::-1]
        dict_atmprof['MU'] = dict_atmprof['MU'][::-1]
        for ii in xrange(len(list_mol)) :
            dict_atmprof[list_mol[ii]] = dict_atmprof[list_mol[ii]][::-1]
    elif (reverse != 0):
        errors.exit_msg("Check the order of FILE_ATM.")
    #------------------------------------------------

    if key :
        return list_mol, dict_atmprof[key]
    else :
        return list_mol, dict_atmprof



##=============================================================================
#def read_lookuptable(xsfile, wn_limit):
#    """
#    Extract Look-up Table
#    """
#    data = np.load(xsfile)
#
#    PP_grid     = data['P']*1.0e3 # mbar => barye (x 1000)
#    TT_grid     = data['T']
#    WN_grid_org = data['WN']
#    XS_grid_org = data['XS']
#
#    #------------------------------------------------
#    # check wavenumber range
#    #------------------------------------------------
#    if ( (wn_limit[0] < min(WN_grid_org)) or (wn_limit[1] > max(WN_grid_org)) ) : 
#        errors.exit_msg("Wavenumbers set are out of range of look-up tables.")
#    #------------------------------------------------
#
#    indx_min = nearestindex_WN(WN_grid_org, wn_limit[0])
#    indx_max = nearestindex_WN(WN_grid_org, wn_limit[1])
#
#    WN_grid = WN_grid_org[indx_min:indx_max+1]
#    del WN_grid_org
#    XS_grid = XS_grid_org[indx_min:indx_max+1]
#    del XS_grid_org
#
#    #### TEST (to be eventually removed)
#    id1, id2, id3 = np.where(XS_grid <= 0)
#    for i1, i2, i3 in zip(id1, id2, id3) :
#        XS_grid[i1][i2][i3] = 1.e-48
#    #### TEST
#
#    return WN_grid, TT_grid, PP_grid, XS_grid
#
#
#
##=============================================================================
#def nearestindex_WN(WN_lattice, wn) :
#    """
#    Returns the wavenumber index nearest to the given wavelength
#    """
#    idx = np.abs(WN_lattice-wn).argmin()
#    return idx



#=============================================================================
# main
#=============================================================================

def main():

    #------------------------------------------------
    # set up atmospheric profile
    #------------------------------------------------
    dict_atmprof, list_mol, tuple_func_atmprof = set_atmprof(FILE_ATM)

    #------------------------------------------------
    # set up wn- and z- grids
    #------------------------------------------------
    layer_z = np.linspace(dict_atmprof['z'][0], dict_atmprof['z'][-1], Z_NUM)
    if (OBSMODE=="nadir" and WN_NUM > 3):
        wn_num = 3
    else:
        wn_num = WN_NUM
    grid_wn = np.linspace(WN_MIN, WN_MAX, wn_num)

    #------------------------------------------------
    # set up extinction profile
    #------------------------------------------------
    mesh_nXS = np.zeros([wn_num, Z_NUM])
    if (MOLABS_ON):
        for molname in list_mol:
            print "adding extinction by ", molname, "..."
            mesh_nXS = mesh_nXS + gasext.calc_nXSofZ_molabs(layer_z, grid_wn, tuple_func_atmprof, molname, XSFILE_TAG, cntnm_on=CNTNM_ON, xsfile_tag_cntnm=XSFILE_TAG_CNTNM)
    if (RAYLEIGH_ON):
        print "adding extinction by Rayleigh scattering..."
        mesh_nXS     = mesh_nXS + gasext.calc_nXSofZ_Rayleigh(layer_z, grid_wn, tuple_func_atmprof, POLARIZATION)
    if (CLD_ON):
        print "adding extinction by clouds..."
        mesh_nXS     = mesh_nXS + cld_simple.calc_nXSofZ_cld(layer_z, grid_wn, tuple_func_atmprof, FILE_CLD)
    #------------------------------------------------

    mesh_nXS[np.where(mesh_nXS < 0)] = 0

    #------------------------------------------------
    # compute spectra
    #------------------------------------------------
    if (OBSMODE=="trans"):
        outfile = outdir+OUTFILE_TAG+".txt"
        obsmode_trans.run(layer_z, grid_wn, mesh_nXS, Rp, TAU_MIN, outfile)
    #------------------------------------------------
    elif (OBSMODE=="disk"):
        obsmode_trans.run(layer_z, grid_wn, mesh_nXS, Rp, TAU_MIN)
    #------------------------------------------------
    elif (OBSMODE=="nadir"):
        obsmode_nadir.run(layer_z, grid_wn, mesh_nXS)
    #------------------------------------------------


if __name__ == "__main__":

    now = datetime.datetime.now()
    print now.strftime("%Y-%m-%d %H:%M:%S")
    main()
    now = datetime.datetime.now()
    print now.strftime("%Y-%m-%d %H:%M:%S")
    print "\007" # beep 
#    cProfile.run('main()')
