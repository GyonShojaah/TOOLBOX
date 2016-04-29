#=============================================================================
# Module
#=============================================================================
import numpy  as np
from   scipy  import integrate
import errors

#=============================================================================
def run(layer_z, grid_wn, mesh_nXS, Rp, tau_min, outfile):
    """
    Obtain spectra
    """

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
#        errors.exit_longmsg(["tauchord_min is greater than TAU_MIN. ", "=> add grids for upper atmosheric layers"])
        print "CAUTION!: tauchord_min is greater than TAU_MIN."
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


