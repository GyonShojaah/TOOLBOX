#------------------------------------------------
# input / output
#------------------------------------------------
MOLECULE='H2O'
ISOTOPE=1

DATAFILE_DIR  = '../data/HITRAN2012/'
DATAFILE='01_hit12.par'

OUTFILE_DIR = '../xstbl/'
OUTFILE_TAG  = 'xstbl_HITRAN2012_02000-04000_m2001_c25_H2O'
#OUTFILE_TAG  = 'tmp'

#------------------------------------------------
# wavelength range [cm^-1]
#------------------------------------------------
WN_MIN    = 2000.0    # cm^-1
WN_MAX    = 4000.0 # cm^-1
#WN_NUM    = 2001
WN_NUM    = 2
WN_CUTOFF = 25.

#------------------------------------------------
# pressure range [mbar]
#------------------------------------------------         
P_MIN     = 1.0e-20 # mbar # might want to go to 1e-6
P_MAX     = 1.0e+5 # mbar # doesn't need to go above 1e4
#P_NUM     = 31
P_NUM     = 26

#------------------------------------------------
# temperature range [K]
#------------------------------------------------                  
#T_MIN     = 50   # K
T_MIN     = 150   # K
T_MAX     = 500  # K
T_NUM     = 8   # K
#T_NUM     = 31   # K
#T_NUM     = 4

