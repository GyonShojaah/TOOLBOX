#=============================================================================
# Module
#=============================================================================

import numpy as np
from scipy import integrate
import datetime

import util_interp
import cgs
import func
import lowres
import errors

MU_ATM     = 28.966  # eventually becomes layer-dependent ?
#TAU_MIN = 1.0e-1
TAU_MIN    = 1.0e-2
TAU_MAX    = 1.0e+1
TAU_NUM    = 50
MU         = 0.5 # cos(60 degree)


#=============================================================================
def calc_sp(func_TofP, funclist_XSofTP, mol_name, atm_profile, wn):
    """
    Executes a the integral of the emission from each atmospheric layers. 
    """
    tau_top, tau_bottom = calc_TAUlimit(func_TofP, funclist_XSofTP, atm_profile, mol_name)
    intensity = 0.0

    #------------------------------------------------
    # check relevance of assumed atmosphere profile
    #------------------------------------------------
    if (tau_top > TAU_MAX):
        errors.exit_longmsg(["tau_top is greater than TAU_MAX. ", "=> add grids for upper atmosheric layers"])
    #------------------------------------------------

    #------------------------------------------------
    # contribution from surface
    #------------------------------------------------
    if (tau_bottom <= TAU_MAX): 
        intensity += func.planck(atm_profile['T'][-1], wn)*np.exp(-1.0*tau_bottom/MU)
    #------------------------------------------------

    #------------------------------------------------
    # contribution from upper atmosphere
    #------------------------------------------------
    if (tau_top > TAU_MIN):
        intensity += func.planck(atm_profile['T'][0], wn)*np.exp(-1.0*tau_top/MU)/MU*tau_top
    #------------------------------------------------

    #------------------------------------------------
    # contribution from atmospheric layers
    #------------------------------------------------
    if (tau_bottom > TAU_MIN): # 
        tau_limit = [np.maximum(tau_top, TAU_MIN), np.minimum(tau_bottom, TAU_MAX)]
        func_integrand = profile_integrand(func_TofP, funclist_XSofTP, mol_name, atm_profile, tau_limit, wn)
        intensity += integrate.quad(func_integrand, np.log(tau_limit[0]), np.log(tau_limit[1]), limit=50)[0]/MU
    #------------------------------------------------

    return intensity


#=============================================================================
def profile_integrand(func_TofP, funclist_XSofTP, mol_name, atm_profile, tau_limit, wn) :
    """
    Returns integral as a function of TAU
    """
    func_PoflogTAU = profile_PoflogTAU(func_TofP, funclist_XSofTP, atm_profile, mol_name, tau_limit)
    #------------------------------------------------
    def returnfunc_integrand(LOGTAUarg) :
        TAUarg = np.exp(LOGTAUarg)
        return func.planck(func_TofP(func_PoflogTAU(LOGTAUarg)), wn)*np.exp(-TAUarg/MU)*TAUarg
    #------------------------------------------------
    return returnfunc_integrand


#=============================================================================
def calc_TAUlimit(func_TofP, funclist_XSofTP, atm_profile, mol_name):
    """
    Returns the optical thickness (TAU) from TOA to the top layer (tau_top) and to the surface (tau_bottom)
    """
    func_nXSofP = profile_nXSofP(funclist_XSofTP, atm_profile, mol_name)
    factor    = cgs.RR/(MU_ATM*cgs.GG)
    #------------------------------------------------
    def func_dTAUdP(TAUarg, PParg):
        return  func_nXSofP(PParg)*func_TofP(PParg)*factor/PParg
    #------------------------------------------------
    tau_top    = func_nXSofP(atm_profile['P'][0])*atm_profile['T'][0]*factor
    tau_bottom = integrate.odeint(func_dTAUdP, tau_top, atm_profile['P'])[-1][0]
    return tau_top, tau_bottom



#=============================================================================
def profile_PoflogTAU(func_TofP, funclist_XSofTP, atm_profile, mol_name, tau_limit) :
    """
    Returns a FUNCTION that estimates the pressure [barye] with a given optical thickness measured from TOA
    """
    func_nXSofP = profile_nXSofP(funclist_XSofTP, atm_profile, mol_name)
    factor      = cgs.RR/(MU_ATM*cgs.GG)
    #------------------------------------------------
    def func_dPdlogTAU(PParg, logTAUarg):
        return 1.0/func_nXSofP(PParg)/func_TofP(PParg)/factor*PParg*np.exp(logTAUarg)
    #------------------------------------------------
    logTAU_lattice = np.linspace(np.log(tau_limit[0]), np.log(tau_limit[1]), TAU_NUM)
    P_MIN = calc_PofTAU(tau_limit[0], func_TofP, funclist_XSofTP, atm_profile, mol_name)
    PP_lattice  = integrate.odeint(func_dPdlogTAU, P_MIN, logTAU_lattice).T[0]
    tmp = util_interp.interp_1d_boundary(logTAU_lattice, PP_lattice, logx=False, order=3)
    return tmp

#=============================================================================
def calc_PofTAU(tau, func_TofP, funclist_XSofTP, atm_profile, mol_name):
    """
    Returns pressure corresponding to tau
    """
    func_nXSofP = profile_nXSofP(funclist_XSofTP, atm_profile, mol_name)
    #------------------------------------------------
    factor = cgs.RR/(MU_ATM*cgs.GG)
    def func_dPdlogTAU(PParg, logTAUarg):
        return 1.0/func_nXSofP(PParg)/func_TofP(PParg)/factor*PParg*np.exp(logTAUarg)
    #------------------------------------------------
    tau_top    = func_nXSofP(atm_profile['P'][0])*atm_profile['T'][0]*factor
    if (tau_top >= tau):
        PofTAU = atm_profile['P'][0]*(tau/tau_top)
    else:
        logTAU_lattice = np.linspace(np.log(tau_top), np.log(tau), TAU_NUM)
        PofTAU = integrate.odeint(func_dPdlogTAU, atm_profile['P'][0], logTAU_lattice)[-1][0]
    return PofTAU


#=============================================================================
def profile_nXSofP(funclist_XSofTP, atm_profile, mol_name):
    """
    Returns a FUNCTION that estimates cummurative ( (number density) x (cross section) ) from TOA with a given pressure [barye]
    """
    PP_layers = atm_profile['P']
    TT_layers = atm_profile['T']
    nXSofP_eachlayer = []
    for ii in xrange(len(funclist_XSofTP)) :
        nXSofP_eachlayer.append(funclist_XSofTP[ii](TT_layers, PP_layers)*atm_profile[mol_name[ii]])
    nXSofP_layers = sum(nXSofP_eachlayer)
    func_nXSofP = util_interp.interp_1d_boundary(PP_layers, nXSofP_layers, order=3)
    return func_nXSofP


