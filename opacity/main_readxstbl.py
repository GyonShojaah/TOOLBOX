#!/usr/bin/env python
import sys, os
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/../common')

import errors
import numpy as np
import datetime

import util_interp

#=============================================================================
# functions
#=============================================================================

#=============================================================================
def read_lookuptable ( infile ):

    """
    Extract Look-up Table
    """

    ncfile_r = netCDF4.Dataset( infile, 'r', format='NETCDF3_64BIT')
    wn_r     = ncfile_r.variables['wn']
    pres_r   = ncfile_r.variables['pres']
    temp_r   = ncfile_r.variables['temp']
    xs_r     = ncfile_r.variables['xs']

    return wn_r, pres_r, temp_r, xs_r


#=============================================================================
# main
#=============================================================================

if __name__ == "__main__":


    filename = sys.argv[1]
    if '.nc' in filename :
        import io_nc
        WN_grid, P_grid, T_grid, XS_grid = io_nc.read_xstbl( filename )
    elif '.npz' in filename :
        data = np.load( filename )
        WN_grid = data['WN']
        P_grid = data['P']
        T_grid = data['T']
        XS_grid = data['XS']
    else :
        # error! 
        errors.exit_msg("Unknown file type of cross section tables.")


    print '--------------------------------------------------'

    print ' select pressure [mbar] : '
    for ii in xrange( len( P_grid ) ) :
        print '({0:d}) {1:e}'.format( ii, P_grid[ii] ), 
    print ''
    s_pres = raw_input()
    i_pres = int( s_pres )

    print '--------------------------------------------------'

    print ' select temperature [K] : '
    for ii in xrange( len( T_grid ) ) :
        print '({0:d}) {1:3f}'.format( ii, T_grid[ii] ), 
    print ''
    s_temp = raw_input()
    i_temp = int( s_temp )

    print '--------------------------------------------------'

    print ' type ouput file name'
    s_temp = raw_input()

    with open( s_temp , 'w') as f:

        for i_wn in xrange( len( WN_grid ) ):
            f.write( str(WN_grid[ i_wn ])+'\t'+str(XS_grid[ i_wn, i_temp, i_pres ])+'\n' )



