#=============================================================================
# Module
#=============================================================================

#from setting_readxstbl import *
import netCDF4
import numpy as np
import datetime
import util_interp
import sys
#import lowres

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
    WN_grid, P_grid, T_grid, XS_grid = read_lookuptable ( filename )

    print '--------------------------------------------------'

    print ' select pressure [mbar] : '
    for ii in xrange( len( P_grid ) ) :
        print '({0:d}) {1:e}'.format( ii, P_grid[ii] ), 
    print ''
    s_pres = raw_input()
    i_pres = int( s_pres )

    print '--------------------------------------------------'

    print ' select temperature [mbar] : '
    for ii in xrange( len( T_grid ) ) :
        print '({0:d}) {1:3f}'.format( ii, T_grid[ii] ), 
    print ''
    s_temp = raw_input()
    i_temp = int( s_temp )

    print '--------------------------------------------------'

    with open( filename+'.tmp' , 'w') as f:

        for i_wn in xrange( len( WN_grid ) ):
            f.write( str(WN_grid[ i_wn ])+'\t'+str(XS_grid[ i_wn, i_pres, i_temp ])+'\n' )



