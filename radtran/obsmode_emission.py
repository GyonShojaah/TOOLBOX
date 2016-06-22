#=============================================================================
# Module
#=============================================================================
import numpy  as np
from   scipy  import integrate
import errors

#=============================================================================
def run(layer_z, grid_wn, mesh_nXS, tuple_func_atmprof, tau_min, outfile):
    """
    Obtain spectra
    """

    func_TofZ, func_PofZ, func_MUofZ, dict_func_NofZ = tuple_func_atmprof

    tau_top, tau_bottom = calc_TAUlimit(mesh_nXS, tuple_func_atmprof)
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
def calc_TAUlimit(layer_z, mesh_nXS, func_TofP, funclist_XSofTP, atm_profile, mol_name):
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






    #------------------------------------------------
    # compute tauchord until it becomes too opaque
    #------------------------------------------------
    list_grid_tauchord = []
    for ii in reversed(xrange(len(layer_z)-1)):
        tauchord_grid = calc_TAUchord(ii, layer_z, mesh_nXS, layer_z[-1], Rp)
        list_grid_tauchord.append(tauchord_grid)
    mesh_tauchord = np.array(list_grid_tauchord).T[:,::-1]

    #------------------------------------------------
    # check relevance of assumed atmosphere profile
    #------------------------------------------------
    if (np.any(mesh_tauchord[:,-1] > tau_min)):
        print mesh_tauchord[np.where(mesh_tauchord[:,-1] > tau_min)][-1]
        errors.exit_longmsg(["tauchord_min is greater than TAU_MIN. ", "=> add grids for upper atmosheric layers"])
    #------------------------------------------------

    #------------------------------------------------
    # compute transit radius
    #------------------------------------------------
    mesh_z        = np.tile(layer_z, (len(mesh_tauchord),1))
    mesh_fext     = (1.0 - np.exp(-mesh_tauchord))*(Rp+mesh_z[:,:-1])
    grid_Fext     = 2.0*integrate.trapz(mesh_fext, x=mesh_z[:,:-1], axis=1)
    grid_Rtransit = np.sqrt(Rp**2+grid_Fext)
    #------------------------------------------------

    grid_Rtransit = (grid_Rtransit-Rp)*1e-5 # cm => km
    grid_wl = 1.0e4/grid_wn
    #------------------------------------------------
    # save rawdata
    #------------------------------------------------
    data = np.dstack([grid_wl, grid_Rtransit])
    np.savetxt(outfile, data[0])
    #------------------------------------------------




#=============================================================================
def calc_TAUchord(zi, layer_z, mesh_nXS, z_max, Rp):

#    indx_zmin = np.where( layer_z > layer_z[zi] + DELTA )[0][0]
    Rarg  = layer_z[zi]
    indx_zmin = zi + 1
    z_min = layer_z[indx_zmin]
#    z_min = layer_z[indx_zmin]
#    z_grid = np.logspace(np.log10(z_min), np.log10(z_max), TAU_NUM)
#
    grid_offset = 2*(Rp+z_min)*np.sqrt(z_min-Rarg) - (4.0/3.0)*(z_min-Rarg)**1.5
    grid_offset = grid_offset * (mesh_nXS[:,zi] + mesh_nXS[:,indx_zmin])/(np.sqrt(2.0*Rp))

    layer_z_upper  = layer_z[indx_zmin:]
    mesh_z_upper   = np.tile(layer_z_upper, (len(mesh_nXS),1))
    mesh_nXS_upper = mesh_nXS[:,indx_zmin:]
    mesh_integrand = (Rp+mesh_z_upper)/np.sqrt(2*(mesh_z_upper-Rarg)*Rp)*mesh_nXS_upper
#    mesh_integrand = (np.sqrt(1-(Rp+Rarg)**2/(Rp+mesh_z_upper)**2) + (Rp+Rarg)**2/(Rp+mesh_z_upper)/np.sqrt((Rp+mesh_z_upper)**2-(Rp+Rarg)**2))*mesh_nXS_upper
#(Rp+mesh_z_upper)/np.sqrt(2*(mesh_z_upper-Rarg)*Rp)*mesh_nXS_upper
    grid_tauchord  = 2.0*integrate.trapz(mesh_integrand, x=mesh_z_upper, axis=1) + grid_offset
    return grid_tauchord


