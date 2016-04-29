import numpy as np
import sys
import datetime
import util_interp
import MT_CKD

import ConfigParser
from cgs import *


from scipy import constants
from scipy import special


MOLNAME  = "H2O"
WN_CUTOFF = 25.0


DATAFILE_HEAD   = 'data/HITRAN2012/HITRAN2012_00000-40000_'

WN_MIN   = 10.0   # cm^-1
WN_MAX   = 10000.0  # cm^-1
WN_NUM   = 100
         


keywords = ['pressure:', 'temperature:', 'molecule:']

#CC = constants.c*1e2
#HH = constants.h*1e7
#KK = constants.k*1e7
#GG = constants.g*1e2
#AMU = 1.661e-24

##############################################################################
def calc_xs(infile, WN_lattice, T_lattice, P_lattice) :
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

    header, data = read_linefile(infile)
    mass = get_moldata(header['molecule:'], "mass")[0]
    xs = np.zeros([len(WN_lattice), len(T_lattice), len(P_lattice)], dtype=np.float64)

    index_min = np.zeros(WN_NUM)
    index_max = np.zeros(WN_NUM)

    for count_WN in xrange(WN_NUM) :
        
#        print "count_WN:", count_WN
#        now = datetime.datetime.now()
#        print now.strftime("%Y-%m-%d %H:%M:%S")

        # check the WN range
        if (data['position'][0] > (WN_lattice[count_WN]-WN_CUTOFF)):
            index_min = 0
        else:
            index_min = np.amin(np.where(data['position'] > (WN_lattice[count_WN]-WN_CUTOFF)))
        if (data['position'][-1] < (WN_lattice[count_WN]+WN_CUTOFF)):
            index_max = len(data['position'])-1
        else:
            index_max = np.amin(np.where(data['position'] > (WN_lattice[count_WN]+WN_CUTOFF) ))-1

#        print "wn", WN_lattice[count_WN]
#        print "index_min, index_max:", index_min, index_max

        if (index_max > index_min) :

            for count_T in xrange(len(T_lattice)) :
                for count_P in xrange(len(P_lattice)) :
                    
                    strength_array  =  linestrength( header['molecule:'], T_lattice[count_T], header['temperature:'], data['position'][index_min:index_max+1], data['energy'][index_min:index_max+1], data['strength'][index_min:index_max+1] ) # unit: [cm^-1 / cm^-2]
                    voigtfunc_array = voigtfunction( T_lattice[count_T], header['temperature:'], P_lattice[count_P]*1e3, header['pressure:'], data['position'][index_min:index_max+1], data['airWidth'][index_min:index_max+1], data['Tdep'][index_min:index_max+1], mass ) # unit: [cm]

                    intensity_array = strength_array*voigtfunc_array(WN_lattice[count_WN])
                    xs[count_WN][count_T][count_P] = np.sum(intensity_array)

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
    moldata.read('data/moldata.cfg')
    tmp = moldata.get(molecule,key)

    return map(float, tmp.split(","))


##############################################################################
def read_linefile(infile) :
    """
    =======================================================================

    read_linefile (infile) :

    -------- USAGE --------------------------------------------------------
    read header and data table of the line data file 
    which is extracted from HITRAN database by 'extract'

    -------- INPUT --------------------------------------------------------
    infile    (str)    line file extracted from HITRAN database by 'extract'

    -------- OUTPUT --------------------------------------------------------
    header    (dict)   molecule, P_0, T_0
    data      (dict)   wl, intensity
    mass      (float)  mass of the molecule

    =======================================================================
    """
    #read header
    header = {}
    with open(infile,'r') as f_linefile :
        for line in f_linefile :
            elements = line.rsplit()
            for j in xrange(len(keywords)) :
                if keywords[j] in elements :
                    if keywords[j] =='pressure:' :
                        if elements[elements.index('pressure:')+2] == 'mb':
                            header['pressure:'] = float(elements[elements.index(keywords[j])+1])*1.0e3
                    elif keywords[j] =='temperature:' :
                        header['temperature:'] = float(elements[elements.index(keywords[j])+1])
                    else :
                        header[keywords[j]] = elements[elements.index(keywords[j])+1]
            if elements[0] != '#' :
                break
    data = np.loadtxt(infile, 
                      dtype={'names':('position','strength','energy','airWidth','Tdep'),
                             'formats':('f8','f8','f8','f8','f8')}, 
                      comments='#')
    return header, data


##############################################################################
def voigtfunction(TT, T_0, PP, P_0, wn, airWidth, Tdep, mass) :
# def voigtfunction(TT, T_0, PP, P_0, wn, airWidth, selfWidth, Tdep, mass, mixratio) :
    """
    =======================================================================

    voigtfunction (TT, T_0, PP, P_0, wn)

    -------- USAGE --------------------------------------------------------
    compute broadening function (voigt profile)
    given temperature and pressure

    -------- INPUT --------------------------------------------------------
    TT        (float)  temperature                             [K]
    T_0       (float)  standard temperature                    [K]
    PP        (float)  total pressure                          [barye]
    P_0       (float)  standard pressure                       [barye]
    wn        (float)  wave number                             [cm-1]
    airWidth  (float)  air-broadening coefficients             [cm-1]
    Tdep      (float)  power of TT depedence of Lorentz width  [-]
    mass      (float)  molecular mass                          [AMU]

    -------- OUTPUT --------------------------------------------------------
    voigtfunc (func)   Voigt Function

    =======================================================================
    """
    #
    # QUESTION!!!
    # p_total = PP
    # p_partial = PP*mixratio
    # gamma_L = (airWidth*((p_total-p_partial)/P_0)+selfWidth*(p_partial/P_0))*(T_0/TT)**Tdep
    #
#    print mass, wn
    gamma_L = airWidth*(PP/P_0)*(T_0/TT)**Tdep
    gamma_D = wn*np.sqrt(2*np.log(2)*KK*TT/(mass*AMU*CC**2))
    #
    # Voigt function
    #
    # integral_-infty^+infty dx' g_L(x-x')*g_D(x'-x0)
    #   pressure broadening : g_L(x) = gamma_L/pi/((x-x0)^2+gamma_L^2)
    #   Doppler broadening  : g_D(x) = (1/gamma_D)*(ln(2)/pi)^(1/2)*exp(-ln(2)*((x-x0)/gamma_D))
    #
    # Voigt function is the real part of the Faddeeva function
    # w(x+iy) = V(x,y) + iL(x,y)
    #
    def return_voigtfunc(x) : 
        realvalue = np.sqrt(np.log(2))/gamma_D * (x-wn)
        imagvalue = np.sqrt(np.log(2))*gamma_L/gamma_D
        complexvalue = realvalue + imagvalue*1j
        faddeeva = np.sqrt(np.log(2))/gamma_D/np.sqrt(np.pi) * special.wofz(complexvalue)
        return faddeeva.real
    return return_voigtfunc




def linestrength(molecule, TT, T_0, wn, energy, strength) :
    """
    =======================================================================

    linestrength(TT, T_0, wn, energy) :

    -------- USAGE --------------------------------------------------------
    compute line strength
    given temperature and pressure

    -------- INPUT --------------------------------------------------------
    TT        (float)  temperature                             [K]
    T_0       (float)  standard temperature                    [K]
    wn        (float)  wavenumber of the transition            [cm-1/cm-2]
    energy    (float)  energy of the transition                [cm-1/cm-2]
    strength  (float)  line strength at standard condition     [cm-1/cm-2]

    -------- OUTPUT --------------------------------------------------------
    new_strength  (float)  new line strength                   [-]

    -------- REFERENCE -----------------------------------------------------
    Norton and Rinsland (1991)

    =======================================================================
    """

    beta = get_moldata(molecule, "beta")
    freq = get_moldata(molecule, "freq")
    dege = get_moldata(molecule, "dege")

    fac_exponent  = np.exp(-HH*CC*energy/(KK*TT))/np.exp(-HH*CC*energy/(KK*T_0))
    fac_exponent *= (1 - np.exp(-HH*CC*wn/(KK*TT)))/(1 - np.exp(-HH*CC*wn/(KK*T_0)))
    fac_Q_rot = (T_0/TT)**beta
    fac_Q_vib = 1.0
    for i in xrange(len(freq)) :
        fac_Q_vib *= (( 1 - np.exp(-(HH*CC*freq[i])/(KK*TT)) )/( 1 - np.exp(-(HH*CC*freq[i])/(KK*T_0)) ))**dege[i]

    new_strength = (fac_exponent*fac_Q_rot*fac_Q_vib)*strength
    return new_strength



#    if data



##############################################################################
if __name__ == "__main__":



    print ""
    print "Pressure? (mbar)"
    pres = input()
    print ""
    print "Temperature? (K)"
    temp = input()
    print ""
    print "computing cross section of " + MOLNAME + " at P = ", str(pres), " mbar, T = " + str(temp) + " K ...."


    outfile = "xs_" + MOLNAME + "_P-" + str(pres) + "_T-" + str(temp)
    now = datetime.datetime.now()
    with open(outfile, 'w') as f :
        f.write(now.strftime("%Y-%m-%d %H:%M:%S") + '\n')

    WN_lattice = np.linspace(         WN_MIN,          WN_MAX,  num=WN_NUM) # cm-1
    P_lattice  = np.array([pres])
    T_lattice  = np.array([temp])

    xs_line_all = calc_xs(DATAFILE_HEAD+MOLNAME, WN_lattice, T_lattice, P_lattice)
    xs_line     = xs_line_all[:,0,0]

    print "xs_line", xs_line

    vi_out, csh2or_out, cfh2or_out = MT_CKD.get_cntnm(pres, temp)
    vi_min = np.where( vi_out > 0 )[0][0]
    vi_max = np.where( vi_out > 1.0e4 )[0][0]
    with open(outfile+"_tmp", 'w') as f :
        for ii in xrange(vi_min, vi_max):
            f.write(str(vi_out[ii])+"\t"+str(csh2or_out[ii])+"\t"+str(cfh2or_out[ii])+"\n")

    vi     = vi_out[vi_min:vi_max]
    csh2or = csh2or_out[vi_min:vi_max]
    cfh2or = cfh2or_out[vi_min:vi_max]
    ch2or = csh2or + cfh2or
    func_MTCKD  = util_interp.interp_1d_boundary(vi, ch2or, order=3, logx=False)
    xs_cont     = func_MTCKD(WN_lattice)

    print "xs_cont", xs_cont


    with open(outfile, 'a') as f :
        for ii in xrange(WN_NUM):
            f.write(str(WN_lattice[ii]) + '\t' + str(xs_line[ii]) + '\t' + str(xs_cont[ii]) + '\n')

    with open(outfile, 'a') as f :
        f.write(now.strftime("%Y-%m-%d %H:%M:%S") + '\n')
