#=============================================================================
# Module
#=============================================================================

from setting import *
import numpy as np
from scipy import integrate
from scipy  import constants
import util_interp
import errors
import cgs

RHOL = 1.0 # density of liquid/ice water
ALPHA_L = 7.0
ALPHA_S = 1.0
Qext = 2.0 # geometric optics

#=============================================================================
# functions
#=============================================================================

#=============================================================================
def calc_nXSofZ_cld(layer_z, grid_wn, tuple_func_atmprof, infile_cld):

    #------------------------------------------------
    # read water content profile
    #------------------------------------------------
    cld_profile = read_cldfile(infile_cld, tuple_func_atmprof)

    #------------------------------------------------
    # compute Qext, Qsca
    #------------------------------------------------
    # cldlayer_nXSext = np.zeros_like(cld_profile['z'])
    # cldlayer_nXSsca = np.zeros_like(cld_profile['z'])

    lattice_r = np.logspace(np.log10(D_MIN), np.log10(D_MAX), D_NUM)
    mesh_r,    mesh_zcld = np.meshgrid(lattice_r, cld_profile['z'])
    mesh_r,    mesh_lwc  = np.meshgrid(lattice_r, cld_profile['LWC'])
    mesh_r,    mesh_swc  = np.meshgrid(lattice_r, cld_profile['SWC'])

#    print "mesh_swc", mesh_swc

    mesh_dNdlogD_L = func_dNdlogD_gamma(mesh_r, D_EFF_L, mesh_lwc, ALPHA_L)
    mesh_dNdlogD_S = func_dNdlogD_gamma(mesh_r, D_EFF_S, mesh_swc, ALPHA_S)

    mesh_nXSext_L  = np.pi*mesh_r*mesh_dNdlogD_L # integration with r
    mesh_nXSext_S  = np.pi*mesh_r*mesh_dNdlogD_S # integration with r

#    print "mesh_nXSext_L", mesh_nXSext_L
#    print "mesh_nXSext_S", mesh_nXSext_S

    mesh_nXSext = mesh_nXSext_L + mesh_nXSext_S

    cldlayer_nXSext = Qext*integrate.trapz(mesh_nXSext, x=mesh_r, axis=1)
#    print "----------cloud data---------"
#    for ii in xrange(len(cld_profile['z'])):
#        print cld_profile['z'][ii], cldlayer_nXSext[ii]


#    print ""
    func_kext = util_interp.interp_1d(cld_profile['z'], cldlayer_nXSext, logx=False, bounds_error=True, kind="linear")
#    func_ksca = util_interp.interp_1d(cld_profile['z'], nXSsca_cldlayers, order=1, logx=False, logy=False, fill_value=0.0)

#    print "----------atmospheric data---------"
    layer_nXSext = np.zeros_like(layer_z)
    for zi in xrange(len(layer_z)):
        if  ( layer_z[zi] < cld_profile['z'][0] or cld_profile['z'][-1] < layer_z[zi] ):
            layer_nXSext[zi] = 0.0
        else:
            layer_nXSext[zi] = func_kext(layer_z[zi])

#        print layer_z[zi], layer_nXSext[zi]


#    nXSsca_layers = func_ksca(grid_z)
    mesh_nXSext = np.tile(layer_nXSext, (len(grid_wn),1))
    return mesh_nXSext


#=============================================================================
#def read_cldfile(infile):
#    names = ['z', 'LWC', 'R_eff']
#    formats = ['f8', 'f8', 'f8']
#    data = np.loadtxt(infile, 
#                      dtype={'names':tuple(names),
#                             'formats':tuple(formats)}, 
#                      comments='#')
#    data['z']   = data['z']*1.e5 # km => cm
#    data['LWC']   = data['LWC']*1.e-6 # g/m3 => g/cm3
#    data['R_eff'] = data['R_eff']*1.e-4 # um => cm
#    return data

#=============================================================================
def read_cldfile(infile, tuple_func_atmprof):

    names = ['z', 'LWC', 'LWC2', 'SWC', 'SWC2']
    formats = ['f8', 'f8', 'f8', 'f8', 'f8']
    cld_profile = np.loadtxt(infile, 
                      dtype={'names':tuple(names),
                             'formats':tuple(formats)}, 
                      comments='#')


    #------------------------------------------------
    # km => cm
    #------------------------------------------------
    cld_profile['z'] = cld_profile['z']*1e5
    #------------------------------------------------

    #------------------------------------------------
    # combine lwc and lwc2
    #------------------------------------------------
    cld_profile['LWC'] = ( cld_profile['LWC'] + cld_profile['LWC2'] ) * 0.5
    cld_profile['SWC'] = ( cld_profile['SWC'] + cld_profile['SWC2'] ) * 0.5
    #------------------------------------------------


    #------------------------------------------------
    # check the order of altitude
    #------------------------------------------------
    z_old = cld_profile['z'][0]
    reverse = 0
    for zi in xrange(1,len(cld_profile['z'])):
        if (cld_profile['z'][zi] < z_old):
            reverse = reverse + 1
    if (reverse == len(cld_profile['z'])-1):
        cld_profile['z'] = cld_profile['z'][::-1]
        cld_profile['LWC'] = cld_profile['LWC'][::-1]
        cld_profile['SWC'] = cld_profile['SWC'][::-1]
    elif (reverse != 0):
        errors.exit_msg("Check the order of FILE_CLD.")
    #------------------------------------------------

    # actually, 'P', 'cldh2o', 'R_eff'
    func_TofZ, func_PofZ, func_MUofZ, dict_func_NofZ = tuple_func_atmprof
    cldlayer_P   = func_PofZ(cld_profile['z'])
    cldlayer_T   = func_TofZ(cld_profile['z'])
    cldlayer_MU  = func_MUofZ(cld_profile['z'])

    cldlayer_rho = cldlayer_MU*cldlayer_P/(cgs.RR*cldlayer_T)
    cld_profile['LWC'] = cld_profile['LWC']*cldlayer_rho
    cld_profile['SWC'] = cld_profile['SWC']*cldlayer_rho

#    print "cld_profile['LWC']", cld_profile['LWC']
#    print "cld_profile['SWC']", cld_profile['SWC']


#    func_RHOofZ = util_interp.interp_1d_boundary(cld_profile['z'], cld_profile['LWC'], array_rho, logx=False, logy=False, order=1)
    return cld_profile



#=============================================================================
def func_dNdlogD_gamma(D, r_eff, wc, alpha):

#    if (wc==0.0):
#        dNdlogD = 0.0*D
#    else:
    b = (alpha + 3.0)/r_eff
    factor = (4.0*np.pi*RHOL/3.0)
    def integrand(Rarg):
        return Rarg**alpha*np.exp(-b*Rarg)*Rarg**2
    a = wc/(factor*integrate.quad(integrand, D_MIN, D_MAX)[0])
    dNdlogD = a*D**alpha*np.exp(-b*D)
    return dNdlogD



#=============================================================================
def func_NofD_exp(D, wc):

    # my thought
    # lambda0 = 
    # N0 = wc*lambda0**4.0/(np.pi*RHOL)
    # Kim & Del Genio (2013)
    N0      = 8.0 # 8.0 [1/cm4] = 8.e6 [1/m4]
    lambda0 = (np.pi*N0*RHOL/wc)**0.25
    NofD    = N0*np.exp(-lambda0*D)
    return NofD



