import numpy as np
import sys
import datetime
import ConfigParser
import MT_CKD
import util_interp

from scipy import constants
from scipy import special


OUTFILE_TAG  = 'xstbl2_cntnm_00010-10000_m09991_H2O'

WN_MIN   = 10.0   # cm^-1
WN_MAX   = 10000.0  # cm^-1
WN_NUM   = 9991
         
P_MIN    = 1.0e-5  # mbar # might want to go to 1e-6
P_MAX    = 1.0e+3  # mbar # doesn't need to go above 1e4
P_NUM    = 9

T_MIN    = 150  # K
T_MAX    = 350  # K
T_NUM    = 5    # K  per 25 K

CC = constants.c*1e2
HH = constants.h*1e7
KK = constants.k*1e7
GG = constants.g*1e2
AMU = 1.661e-24




##############################################################################
def calc_xs_cntnm(WN_lattice, T_lattice, P_lattice) :
    """
    =======================================================================

    calc_xs (infile, wl_lattice, P_lattice, T_lattice) :

    -------- USAGE --------------------------------------------------------
    computation of line-by-line molecular absorption cross sections
    give 3D (wavelength, temperature and pressure) grids

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
    xs = np.zeros([WN_NUM, T_NUM, P_NUM], dtype=np.float64)

    for count_T in xrange(T_NUM) :
        for count_P in xrange(P_NUM) :

            vi_out, csh2or_out, cfh2or_out = MT_CKD.get_cntnm(P_lattice[count_P], T_lattice[count_T])
            vi    = np.zeros(len(np.where(csh2or_out > 0)[0]))
            cntnm = np.zeros(len(np.where(csh2or_out > 0)[0]))
            jj = 0
            for ii in np.where(csh2or_out > 0)[0]:
                vi[jj]    = vi_out[ii]
                cntnm[jj] = csh2or_out[ii]
                jj        = jj + 1

            if any(cntnm[1:] <= 0 ):
                print "error"
                sys.exit()

            func_XSofV = util_interp.interp_1d_boundary(vi, cntnm, logx=False)
            grid_XS = func_XSofV(WN_lattice)

            for count_WN in xrange(WN_NUM):
                xs[count_WN][count_T][count_P] = grid_XS[count_WN]

    return xs


##############################################################################
def get_moldata(molecule, key):
    """
    =======================================================================

    get_moldata(molecule, key) :

    -------- USAGE --------------------------------------------------------
    get 'key' parameter of the molecule

    -------- INPUT --------------------------------------------------------
    molecule  (str)    molecule name           
    key       (str)    keyword ("ID, "mass", "beta", "freq", "dege")

    -------- OUTPUT --------------------------------------------------------
    valule    (float)  the value of the keyword

    =======================================================================
    """
    moldata = ConfigParser.ConfigParser()
    moldata.read('moldata.cfg')
    tmp = moldata.get(molecule,key)

    return map(float, tmp.split(","))




##############################################################################
if __name__ == "__main__":

    outfile = OUTFILE_TAG + '.npz'
    logfile = OUTFILE_TAG + '.log'

    now = datetime.datetime.now()
    with open(logfile, 'w') as f :
        f.write(now.strftime("%Y-%m-%d %H:%M:%S") + '\n')

    WN_lattice = np.linspace(         WN_MIN,          WN_MAX,  num=WN_NUM) # cm-
    P_lattice  = np.logspace(np.log10(P_MIN), np.log10(P_MAX),  num=P_NUM)  # mbar
    T_lattice  = np.linspace(          T_MIN,           T_MAX,  num=T_NUM)  # K

    xs = calc_xs_cntnm(WN_lattice, T_lattice, P_lattice)

    np.savez(outfile, WN=WN_lattice, P=P_lattice, T=T_lattice, XS=xs)

    now = datetime.datetime.now() 
    with open(logfile, 'a') as f :
        f.write('WLNUM: ' + str(WN_NUM) + '\n')
        f.write('PNUM : ' + str(P_NUM) + '\n')
        f.write('TNUM : ' + str(T_NUM) + '\n')
        f.write(now.strftime("%Y-%m-%d %H:%M:%S") + '\n')


