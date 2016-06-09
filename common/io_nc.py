import netCDF4

#--------------------------------------------------------------------
def save_xstbl( outfile, wn_lattice, p_lattice, t_lattice, xs ):

    ncfile_w = netCDF4.Dataset( outfile,'w', format='NETCDF3_64BIT') 

    ncfile_w.createDimension( 'wn',    len(wn_lattice) )
    ncfile_w.createDimension( 'pres',  len(p_lattice)  )
    ncfile_w.createDimension( 'temp',  len(t_lattice)  )

    wn_w   = ncfile_w.createVariable( 'wn', 'float32', ('wn',))
    wn_w.title = 'Wavenumber'
    wn_w.units = 'cm^-1'
    wn_w[:] = wn_lattice

    pres_w = ncfile_w.createVariable( 'pres', 'float32', ('pres',))
    pres_w.title = 'Pressure'
    pres_w.units = 'mbar'
    pres_w[:] = p_lattice

    temp_w = ncfile_w.createVariable( 'temp', 'float32', ('temp',))
    temp_w.title = 'Temperature'
    temp_w.units = 'K'
    temp_w[:] = t_lattice

    xs_w   = ncfile_w.createVariable( 'xs', 'float32',  ('wn','pres','temp') )
    xs_w.title = 'Cross Section'
    xs_w.units = 'cm^2'
    xs_w[:] = xs

    ncfile_w.close()



#--------------------------------------------------------------------
def read_xstbl( infile ):

    ncfile_r = netCDF4.Dataset( infile,'r', format='NETCDF3_64BIT') 

    wn_grid = ncfile_r.variables['wn'][:]
    p_grid  = ncfile_r.variables['pres'][:]
    t_grid  = ncfile_r.variables['temp'][:]
    xs_grid = ncfile_r.variables['xs'][:,:,:]

    return wn_grid, p_grid, t_grid, xs_grid
