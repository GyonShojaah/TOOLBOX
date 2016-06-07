#------------------------------------------------
# input / output
#------------------------------------------------
MOLECULE='H2O'
ISOTOPE=1

DATAFILE_DIR  = '/Users/yuka/data/HITRAN2012/'
DATAFILE='01_hit12.par'

OUTFILE_DIR = 'xstbl/'
OUTFILE_TAG  = 'xstbl_HITRAN2012_00010-2000_m10000_H2O'
#OUTFILE_TAG  = 'tmp'

#------------------------------------------------
# wavelength range [cm^-1]
#------------------------------------------------
WN_MIN    = 10.0    # cm^-1
WN_MAX    = 10000.0 # cm^-1
WN_NUM    = 1000
WN_CUTOFF = 25.

#------------------------------------------------
# pressure range [mbar]
#------------------------------------------------         
P_MIN     = 1.0e-20 # mbar # might want to go to 1e-6
P_MAX     = 1.0e+10 # mbar # doesn't need to go above 1e4
#P_NUM     = 31
P_NUM     = 4

#------------------------------------------------
# temperature range [K]
#------------------------------------------------                  
#T_MIN     = 50   # K
T_MIN     = 150   # K
T_MAX     = 550  # K
T_NUM     = 8   # K
#T_NUM     = 31   # K
#T_NUM     = 4

