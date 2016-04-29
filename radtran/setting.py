#------------------------------------------------
# computation mode
#------------------------------------------------

OBSMODE     = "trans"
#OBSMODE     = "nadir"

# Include molecular absorption ?
MOLABS_ON   = True

# Include molecular absorption ?
CNTNM_ON   = False

# Include Rayleigh scattering ?
RAYLEIGH_ON = False

# Include cloud scattering ?
CLD_ON      = False

#------------------------------------------------
# computational parameters
#------------------------------------------------

# CHECK CONVERGENCE BY INCREASING THIS NUMBER !!
# grid point for integration
# (used if GRID_ON = True)
Z_NUM     = 1000
#Z_MAX     = 2.0e7 # cm

# CHECK CONVERGENCE BY DECREASING THIS NUMBER !!
# offset (float)
DELTA       = 1.0

# CHECK CONVERGENCE BY INCREASING THIS NUMBER !!
# wavenumber grids
# (used if USER_WN_ON = True)
WN_MIN      = 1000
WN_MAX      = 6000
WN_NUM = 10000
#WN_MIN = 1000
#WN_MAX = 10000
#WN_NUM = 9991


# CHECK CONVERGENCE BY INCREASING THIS NUMBER !!
# cloud particle grids
# (used if CLD_ON = True)
D_MIN   = 1e-8 # cm
D_MAX   = 1e-0 # cm
D_NUM   = 100
D_EFF_L   = 1.0e-3 # cm
D_EFF_S   = 3.0e-3 # cm
ANG_NUM = 10
RHOL    = 1.0 # g/cm


#------------------------------------------------
# input/output files
#------------------------------------------------

#FILE_ATM    = "data/atmprof/prof_segura2005_Temp_H2O_CH4"
FILE_ATM    = "data/atmprof/prof_segura2005_CH4"
#FILE_ATM    = "data/atmprof/prof_midlatsummer_ppmv"
#FILE_ATM    = "data/atmprof/prof_Gmin_T400K_H2O_CH4"
#FILE_ATM    = "data/atmprof/prof_ANN2261-2264aijlE_g10_R10_P10_tl_P4days"
#FILE_ATM    = "data/atmprof/profile_Neptune_T300K_P100bar"
#FILE_ATM    = "data/atmprof/profile_HATP11b_T500K"
#FILE_ATM    = "data/atmprof/profile_H2O_TEMP400_PSURF1e3_MIX1e-7"
#FILE_HEIGHT = "data/cld/"
#FILE_CLD    = "data/atmprof/prof_ANN2261-2264aijlE_g10_R10_P10_tl_P4days_cld"
#FILE_CLD     = "data/cld/WC2.DAT"
#FILE_HEIGHT = "data/cld/aqua_tl_P30days_92.5E_height"
#FILE_RFINDX = "data/cld/Hale.yml.txt"    # liquid water
#FILE_RFINDX = "data/cld//Warren.yml.txt" # ice water

#XSFILE_TAG  = "../xstbl/xstbl_HITRAN2012_00010-10000_m09991_c25_"
XSFILE_TAG  = "../xstbl/xstbl2_HITRAN2012_00010-10000_m09991_"
#XSFILE_TAG  = "lkuptbl/xstbl_HITRAN2012_20000-30000_m10001_c25_"
#XSFILE_TAG_CNTNM = "../xstbl/xstbl2_cntnm_00010-10000_m09991_"
XSFILE_TAG_CNTNM = "../xstbl/xstbl_cntnm_00010-10000_m09991_"
#OUTFILE_TAG = "cntnm_00010-10000_m10001_H2O_TEMP300_PSURF1e3_MIX1e-6"
#OUTFILE_TAG = "prof_Gmin_T400K_H2O_CH4_Rp3e9_01000-10000_m09991_trapz-1000_c25+cntnm"
#OUTFILE_TAG = "prof_GCM_S0.5_01000-10000_m09991_trapz-1000"
#OUTFILE_TAG = "Neptune_T300K_P100bar_00010-10000_m09991_quad-1000_line"
#OUTFILE_TAG = "prof_ANN1950-1952aijlE_g10_R10_P10_P15days_cld_trapz-1000_line"
#OUTFILE_TAG = "profile_HATP11b_T500K_00010-10000_m09991_quad-1000_line+cntnm"
#OUTFILE_TAG = "Earth_midlatsummer"
OUTFILE_TAG = "segura2005_CH4"

#------------------------------------------------
# planetary parameters (float)
#------------------------------------------------

Rp           = 6.371e8 # cm
#Rp            = 49.528e8 # cm
#Rp            = 30.0e8
#Rp            = 28.854e8
#Rp           = 6.371e6 # m
POLARIZATION = 1.710e-24 # cm^3 
#POLARIZATION = 1.710e-30 # m^3
cosTH0       = 0.5 # cos(60 degree)
ALPHA        = 1.0

#------------------------------------------------
# internal settings (float)
#------------------------------------------------

TAU_MIN = 0.02
#TAU_MAX = 10.0




