#=============================================================================
# Module
#=============================================================================

from setting_sp import *
import numpy as np
from scipy import integrate
from scipy  import constants
import util_interp
import errors
import cgs
#import func
#from copy import deepcopy
from bhmie_herbert_kaiser_july2012 import *

#=============================================================================
# functions
#=============================================================================

#=============================================================================
def profile_nXSofZ_cld(layer_z, grid_wn, infile_cld, infile_rfindx, alpha):

    #------------------------------------------------
    # read refraction index
    #------------------------------------------------
    func_rfindx_re, func_rfindx_im = read_rfindx(infile_rfindx)
    grid_rfindx    = func_rfindx_re(grid_wn) + func_rfindx_im(grid_wn)*1j

    #------------------------------------------------
    # read water content profile
    #------------------------------------------------
    if (CLD_GCM):
        print 'reading gcm output...'
        cld_profile = read_wcfile_gcm(FILE_HEIGHT, FILE_CLD, ilat=0)
    else:
        cld_profile = read_wcfile(infile_cld)

    #------------------------------------------------
    # prepare Qext, Qsca
    #------------------------------------------------
    lattice_r = np.logspace(np.log10(D_MIN), np.log10(D_MAX), D_NUM)
    lattice_Qext = np.zeros(D_NUM)
    lattice_Qsca = np.zeros(D_NUM)
    for ii in xrange(D_NUM):
        xx = 2.0*np.pi*wn*lattice_r[ii]
        s1,s2,qext,qsca,qback,gsca = bhmie(xx,rfindx,ANG_NUM)
        lattice_Qext[ii] = qext
        lattice_Qsca[ii] = qsca

    #------------------------------------------------
    # compute Qext, Qsca
    #------------------------------------------------
    cldlayer_nXSext = np.zeros(len(cld_profile['z']))
    cldlayer_nXSsca = np.zeros(len(cld_profile['z']))

    mesh_r,    mesh_zcld = np.meshgrid(lattice_r, cld_profile['z']))
    mesh_r,    mesh_lwc  = np.meshgrid(lattice_r, cld_profile['LWC'])
    mesh_Qext, mesh_lwc  = np.meshgrid(lattice_Qext, cld_profile['LWC'])
    mesh_dNdlogD = func_dNdlogD_gamma(mesh_r, D_EFF, mesh_lwc, alpha)
    mesh_nXSext  = np.pi*mesh_r**2*mesh_Qext*mesh_dNdlogD/mesh_r # integration with r
    cldlayer_nXSext = integrate.trapz(mesh_nXSext, x=mesh_r, axis=1)

    func_kext = util_interp.interp_1d(cld_profile['z'], cldlayer_nXSext, order=1, logx=False, logy=False, fill_value=0.0)
#    func_ksca = util_interp.interp_1d(cld_profile['z'], nXSsca_cldlayers, order=1, logx=False, logy=False, fill_value=0.0)

    layer_nXSext = func_kext(layer_z)
#    nXSsca_layers = func_ksca(grid_z)
 
    return layer_nXSext


#=============================================================================
#def read_wcfile(infile):
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
def read_wcfile_gcm(infile_height, infile_cld, atm_profile, ilat=0):

    # actually, 'P', 'cldh2o', 'R_eff'
    func_TofZ, func_PofZ, func_MUofZ, dict_func_NofZ = tuple_func_atmprof
    layer_P = func_PofZ(layer_z)
    layer_T = func_TofZ(layer_z)
    layer_MU = func_MUofZ(layer_z)
    layer_n0 = layer_MU*layer_P/(cgs.RR*layer_T)*constants.N_A

    data = {}
    cld_profile['z']     = np.loadtxt(infile_height).T[ilat]
    cld_profile['LWC']   = np.loadtxt(infile_cld).T[ilat]

    #------------------------------------------------
    # check the order of altitude
    #------------------------------------------------
    reverse = 0
    z_old = cld_profile['z'][0]
    for zi in xrange(1,len(cld_profile['z'])):
        if (cld_profile['z'][zi] < z_old):
            reverse = reverse + 1
    if (reverse == len(cld_profile['z'])-1):
        cld_profile['z'] = cld_profile['z'][::-1]
        cld_profile['LWC'] = cld_profile['LWC'][::-1]
    elif (reverse != 0):
        errors.exit_msg("Check the order of FILE_CLD.")

#    func_RHOofZ = util_interp.interp_1d_boundary(cld_profile['z'], cld_profile['LWC'], array_rho, logx=False, logy=False, order=1)
#    cld_profile['LWC'] = cld_profile['LWC']*func_RHOofZ(grid_z)

    return data



#=============================================================================
def read_rfindx(infile):

    names = ['wl', 're', 'im']
    formats = ['f8', 'f8', 'f8']
    rfindx = np.loadtxt(infile, 
                      dtype={'names':tuple(names),
                             'formats':tuple(formats)}, 
                      comments='#')
    wn_grid = 1.0/(rfindx['wl'][::-1]*1e-4)
    re_grid = rfindx['re'][::-1]
    im_grid = rfindx['im'][::-1]
    func_rfindx_re = util_interp.interp_1d_boundary(wn_grid, re_grid, order=3)
    func_rfindx_im = util_interp.interp_1d_boundary(wn_grid, im_grid, order=3)
    return func_rfindx_re, func_rfindx_im
    

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


#=============================================================================
def func_dNdlogD_gamma(D, r_eff, wc, alpha):

    if (wc==0.0):
        dNdlogD = 0.0*D
    else:
        b = (alpha + 3.0)/r_eff
        factor = (4.0*np.pi*RHOL/3.0)
        def integrand(Rarg):
            return Rarg**alpha*np.exp(-b*Rarg)*Rarg**2
        a = wc/(factor*integrate.quad(integrand, D_MIN, D_MAX)[0])
        dNdlogD = a*D**alpha*np.exp(-b*D)
    return dNdlogD



