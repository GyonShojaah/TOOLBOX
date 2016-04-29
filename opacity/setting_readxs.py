#------------------------------------------------
# input parameters
#------------------------------------------------

TEMP = 200.0   # K
PRES = 1000.0 # mbar

#------------------------------------------------
# computation mode
#------------------------------------------------

# Use user-defined wavenumber grids ? 
# (If False, wavnumber grids of look-up tables are used)
USER_WN_ON = False

CNTNM_ON = True

#------------------------------------------------
# computational parameters
#------------------------------------------------

# CHECK CONVERGENCE BY INCREASING THIS NUMBER !!
# wavenumber grids
# (used if USER_WN_ON = True)
WN_MIN     = 7000
WN_MAX     = 8000
WN_NUM     = 1000

# spectral resolution
# (used if LOW_RES_ON = True)
RESOLUTION = 100


#------------------------------------------------
# input/output files
#------------------------------------------------

XSFILE      = "../xstbl/xstbl_HITRAN2012_00010-10000_m09991_c25_H2O.npz"
XSFILE_CNTNM      = "../xstbl/xstbl_cntnm_00010-10000_m09991_H2O.npz"
OUTFILE_TAG = "out/xstbl_00010-10000_m09991_H2O_T200K_P1000mbar_line+cntnm"
#OUTFILE_TAG = "profile_H2O_09000-10000_quad-50"


#------------------------------------------------
# planetary parameters (float)
#------------------------------------------------

Rp           = 6.371e8 # cm
#Rp           = 6.371e6 # m
POLARIZATION = 1.710e-24 # cm 
#POLARIZATION = 1.710e-30 # cm 
cosTH0       = 0.5 # cos(60 degree)


#------------------------------------------------
# internal settings (float)
#------------------------------------------------

TAU_MIN = 0.01
TAU_MAX = 10.0




