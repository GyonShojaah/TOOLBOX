#!/usr/bin/env python
import sys, os
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/../common')

import numpy as np
import sys
import ConfigParser

from scipy import constants
from scipy import special

from setup_makexstbl import *
import io_txt
import io_nc
from extract_lines import extract_hitran
from lineprofile import *
from molecules import molecules


# reference temperature
TRef, pRef = 296.0, 1013.25 # K, mbar

##############################################################################
def calc_xs ( molecule, linedata, WN_lattice, T_lattice, P_lattice  ) :
    """
    =======================================================================

    calc_xs (infile, wl_lattice, P_lattice, T_lattice) :

    -------- USAGE --------------------------------------------------------
    computation of line-by-line molecular absorption cross sections
    given 3D grids(wavelength, temperature and pressure)

    -------- INPUT --------------------------------------------------------
    line file  (str)                        input file extracted from HITRAN database by 'extract'
    WN_lattice (numpy.array, dim=(WN_NUM))  wavenumber grid  [cm^-1]
    T_lattice  (numpy.array, dim=(T_NUM))   temperature grid [K]
    P_lattice  (numpy.array, dim=(P_NUM))   pressure grid    [mbar]

    -------- OUTPUT --------------------------------------------------------
    xs_array   (numpy.array, dim=(WN_NUM, T_NUM, P_NUM)) cross section matrix  [cm^2 / cm]
    n[cm^-3] * sigma[] * l[cm] = [cm]

    =======================================================================
    """

    mass = molecules[molecule]['mass']
    xs = np.zeros([WN_NUM, T_NUM, P_NUM], dtype=np.float64)

    index_min = np.zeros(WN_NUM)
    index_max = np.zeros(WN_NUM)
    
    for count_T in xrange(T_NUM) :
        for count_P in xrange(P_NUM) :
                    
            strength_array  =  linestrength( molecule, T_lattice[count_T], TRef, linedata['position'], linedata['strength'], linedata['energy'] ) # unit: [cm^-1 / cm^-2]
            voigtfunc_array = voigtfunction( T_lattice[count_T], TRef, P_lattice[count_P]*1e3, pRef, linedata['position'], linedata['airWidth'], linedata['Tdep'], mass ) # unit: [cm]

            for count_WN in xrange(WN_NUM) :                

                if 'WN_CUTOFF' in globals():
                    wn_llim = WN_lattice[count_WN] - WN_CUTOFF
                    wn_ulim = WN_lattice[count_WN] + WN_CUTOFF
                    i_llim  = np.where(  wn_llim < linedata['position'] )[0][0]
                    i_ulim  = np.where(  linedata['position'] < wn_ulim )[0][-1]
                    if i_ulim == len( linedata['position'] ) - 1 :
                        intensity_array = strength_array[i_llim:]*voigtfunc_array(WN_lattice[count_WN])[i_llim:]
                    else :
                        intensity_array = strength_array[i_llim:i_ulim]*voigtfunc_array(WN_lattice[count_WN])[i_llim:i_ulim]
                    xs[count_WN][count_T][count_P] = np.sum(intensity_array)

                else :
                    intensity_array = strength_array*voigtfunc_array(WN_lattice[count_WN])
                    xs[count_WN][count_T][count_P] = np.sum(intensity_array)

    return xs





##############################################################################
if __name__ == "__main__":

    # read input HITRAN data
    infile = DATAFILE_DIR + DATAFILE
    if 'WN_CUTOFF' in globals():
        wn_limit = [ WN_MIN-WN_CUTOFF, WN_MAX+WN_CUTOFF ]
    else:
        wn_limit = [ 0., 10000. ]
    linedata = extract_hitran ( infile, wn_limit, MOLECULE, ISOTOPE )

    # filename
    outfile = OUTFILE_DIR + OUTFILE_TAG + '.nc'
    logfile = OUTFILE_DIR + OUTFILE_TAG + '.log'

    # make grid points
    WN_lattice = np.linspace(         WN_MIN,          WN_MAX,  num=WN_NUM) # cm-1
    P_lattice  = np.logspace(np.log10(P_MIN), np.log10(P_MAX),  num=P_NUM)  # mbar
    T_lattice  = np.linspace(          T_MIN,           T_MAX,  num=T_NUM)  # K

    # save information in log file
    io_txt.save_line( logfile, 'FILE: ' + OUTFILE_TAG + '.nc', addition=False )
    io_txt.save_line( logfile, '' )
    io_txt.save_line( logfile, 'wavenumber [cm^-1] : '+str(WN_MIN)+'-'+str(WN_MAX)+', '+str(WN_NUM)+' grid points' )
    io_txt.save_line( logfile, 'pressure [mbar]    : '+str( P_MIN)+'-'+str( P_MAX)+', '+str( P_NUM)+' grid points' )
    io_txt.save_line( logfile, 'temperature [K]    : '+str( T_MIN)+'-'+str( T_MAX)+', '+str( T_NUM)+' grid points' )
    io_txt.save_line( logfile, '' )
    io_txt.save_time( logfile, 'start : ' )

    # calculate cross section
    xs = calc_xs ( MOLECULE, linedata, WN_lattice, T_lattice, P_lattice  )

    # save output
    io_nc.save_xstbl( outfile, WN_lattice, P_lattice, T_lattice, xs )
    #np.savez(outfile, WN=WN_lattice, P=P_lattice, T=T_lattice, XS=xs)

    # save information in log file
    io_txt.save_time( logfile, 'end   : ' )



