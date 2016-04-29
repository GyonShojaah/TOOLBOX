import numpy as np
from scipy.optimize import minimize
from scipy import constants
import ConfigParser
import random
import matplotlib.pyplot as plt
import cgs
plt.ion()

MU_ATM = 2.0 # g
GRAV   = cgs.GG
TEMP = 1500.0
P_SURF = 1.01325e5 # bar
P_LIMIT = 0.01
P_NUM = 20
LOWERLIMIT = 1e-12
INVALID = 1e12

SH = cgs.RR*TEMP/(MU_ATM*GRAV) # approximation
##############################################################################
def read_cfgfile(filename):

    dict_mol = {}
    cfgfile = ConfigParser.ConfigParser()
    cfgfile.read(filename)
    molname = cfgfile.sections()

    for ii in xrange(len(molname)):
        dict_mol[ii] = {}
        dict_mol[ii]['name'] = molname[ii]
        list_option = cfgfile.options(molname[ii])

        for jj in xrange(len(list_option)) :
            tmp = cfgfile.get(molname[ii],list_option[jj])
            if ',' in tmp :
                tmp_list = map(float, tmp.split(","))
                if 0.0 in tmp_list:
                    tmp_list.remove(0.0)
                dict_mol[ii][list_option[jj]] = tmp_list
            else :
                dict_mol[ii][list_option[jj]] = float(tmp) 

    return dict_mol


##############################################################################
def get_atmlist(dict_mol):

    set_atm = set([])
    for key in dict_mol:
        set_atm = set_atm | set(dict_mol[key]['atm'])
    return list(set_atm)


##############################################################################
def get_matA(dict_mol, list_atm, dict_solar, PP, TT):

    matA = np.zeros([len(list_atm),len(dict_mol)])
    for jj in xrange(len(dict_mol)) :
        for ii in  xrange(len(dict_mol[jj]['atm'])):
            idx = list_atm.index(dict_mol[jj]['atm'][ii])
            matA[idx][jj] = dict_mol[jj]['num'][ii]

    vecB = np.zeros(len(list_atm))
    for ii in xrange(len(list_atm)):
        vecB[ii] = dict_solar[int(list_atm[ii])]

    matA  = (matA.T/vecB).T
    return matA


##############################################################################
def gibbs_energy(array_n, *arg):

    dict_mol, PP, TT = arg

    print dict_mol[2]['g0']/(constants.R*TT), np.log(PP/1.0e3), np.log(array_n[0])
    gibbs = 0.0
    for ii in xrange(len(dict_mol)):
        if (array_n[ii] < LOWERLIMIT):
            gibbs = INVALID
    if (gibbs < INVALID/10.0):
        for ii in xrange(len(dict_mol)):
            if ( array_n[ii] > LOWERLIMIT ) :
                gibbs = gibbs + array_n[ii]*(dict_mol[ii]['g0']/(constants.R*TT) + np.log(PP/1.0e3) + np.log(array_n[ii])) # 
#    print array_n
#    print "sum:", np.array(sum(array_n)-1)
#    print 'gibbs:', gibbs
    return gibbs


##############################################################################
def gibbs_energy_deriv(array_n, *arg):

    dict_mol, PP, TT = arg
    list_gibbs_energy = []
    for ii in xrange(len(dict_mol)):
        term = (dict_mol[ii]['g0']/(constants.R*TT) + np.log(PP/1.0e5) + array_n[ii]) + 1.0
        list_gibbs_energy.append(term)
    return np.array(list_gibbs_energy)


##############################################################################
def constraints(dict_mol, matA):

    list_cons = []
    list_cons.append({'type': 'eq',
                      'fun' : lambda x: sum(x) - 1.0, 
                      'jac' : lambda x: np.ones(len(x)) })

    for ii in xrange(1,len(matA)):
        list_cons.append({'type': 'eq',
                          'fun' : lambda x, _ii=ii: np.dot(matA[_ii], x) - np.dot(matA[0], x),
                          'jac' : lambda x, _ii=ii: matA[_ii] - matA[0]})

#        list_cons[ii] = ({'type': 'eq',
#                          'fun' : lambda x: np.dot(matA[1], x) - np.dot(matA[0], x),
#                          'jac' : lambda x: matA[1] - matA[0] })

#    for ii in xrange(1,len(matA)):
#        list_cons[ii] = ({'type': 'eq',
#                          'fun' : lambda x, ii_tmp=ii: np.dot(matA[ii_tmp], x) - np.dot(matA[0], x),
#                          'jac' : lambda x, ii_tmp=ii: matA[ii_tmp] - matA[0] })

#    for ii in xrange(1,len(matA)):
#        list_cons[ii] = ({'type': 'eq',
#                          'fun' : (lambda ii_tmp: lambda x: np.dot(matA[ii_tmp], x) - np.dot(matA[0], x))(ii),
#                          'jac' : (lambda ii_tmp: lambda x: matA[ii_tmp] - matA[0])(ii) })

#    for jj in xrange(1,len(matA)):
#        list_cons.append({'type': 'eq',
#                          'fun' : lambda x: np.array(np.dot(matA, x)[jj] - np.dot(matA, x)[0]),
#                          'jac' : lambda x: matA[jj] - matA[0]})

#    for jj in xrange(len(dict_mol)):
#        gibbs_energy_jac = make_gibbs_energy_jac(len(dict_mol), jj)
#        list_cons.append({'type': 'ineq',
#                          'fun' : lambda x: np.array(x[jj]),
#                          'jac' : gibbs_energy_jac})
    return list_cons

#def make_gibbs_energy_jac(num, jj):
#    def returngibbs_energy(x):
#        myarray = np.zeros(num)
#        myarray[jj] = 1.0
#        print "myarray", myarray
#        return myarray
#    return returngibbs_energy


def myplot(dict_mol, matA):

    numpoint = 31
    PRS = 1e5 # mbar
    TEMP = 300.0 # K
    H2vec  = np.logspace(-1, 0, numpoint-1)
    H2Ovec = np.logspace(-7, 0, numpoint)
    nO2 = 1e-3
    nCO2 = 1e-3
    H2A, H2OA = np.meshgrid(H2vec, H2Ovec)
    GA = np.zeros_like(H2OA)
    args=(dict_mol, PRS, TEMP)
    for ii in xrange(len(H2vec)):
        for jj in xrange(len(H2Ovec)):
            nH2  = H2A[jj,ii]
            nH2O = H2OA[jj,ii]
            array_n = np.log(np.array([nH2, nH2O, nO2, nCO2]))
            GA[jj,ii] = gibbs_energy(array_n, *args)

#    print "matA", matA

    plt.figure()    
    plt.contourf(H2OA, H2A, np.log10(-1.0*GA))
    gca = plt.gca()
    plt.xlabel("H2O")
    plt.ylabel("H2")
    gca.set_xscale("log")
    gca.set_yscale("log")
    plt.colorbar()

#    x1 = np.array(H2vec)
#    y1 = 1 - x1 - nO2 - nCO2
#    plt.plot(x1, y1, "k", linewidth=2)

# test
    xH2  = 0.99904689
    xH2O =  0.00020181154
    xO2  =  0.00002797387
    xCO2 = 0.00072332455
    print "-------- test --------"
    print "each term", matA[1,0]*xH2, matA[1,1]*xH2O, matA[1,2]*xO2, matA[1,3]*xCO2, matA[0,0]*xH2, matA[0,1]*xH2O, matA[0,2]*xO2, matA[0,3]*xCO2
    print "term1", (matA[1,0]*xH2 + matA[1,1]*xH2O + matA[1,2]*xO2 + matA[1,3]*xCO2) - (matA[0,0]*xH2 + matA[0,1]*xH2O + matA[0,2]*xO2 + matA[0,3]*xCO2)
    print "term2", (matA[2,0]*xH2 + matA[2,1]*xH2O + matA[2,2]*xO2 + matA[2,3]*xCO2) - (matA[0,0]*xH2 + matA[0,1]*xH2O + matA[0,2]*xO2 + matA[0,3]*xCO2)
    print "term3", xH2 + xH2O + xO2 + xCO2

    xH2, xH2O, xO2, xCO2 = [0.1, 0.949105113338, 0.000759552102057, -0.47441722211]
#[0.107977516233, 0.945112955711, 0.000762437501687, -0.472420629185]
    print xH2, xH2O, xO2, xCO2
    print "-------- test --------"
    print "each term", matA[1,0]*xH2, matA[1,1]*xH2O, matA[1,2]*xO2, matA[1,3]*xCO2, matA[0,0]*xH2, matA[0,1]*xH2O, matA[0,2]*xO2, matA[0,3]*xCO2
    print "term3", xH2 + xH2O + xO2 + xCO2 - 1.0
    print "term1", (matA[1,0]*xH2 + matA[1,1]*xH2O + matA[1,2]*xO2 + matA[1,3]*xCO2) - (matA[0,0]*xH2 + matA[0,1]*xH2O + matA[0,2]*xO2 + matA[0,3]*xCO2)
    print "term2", (matA[2,0]*xH2 + matA[2,1]*xH2O + matA[2,2]*xO2 + matA[2,3]*xCO2) - (matA[0,0]*xH2 + matA[0,1]*xH2O + matA[0,2]*xO2 + matA[0,3]*xCO2)


#    y = nH2O
#    x = nH2
#
#    matA[1,0]*nH2 + matA[1,1]*nH2O + matA[1,2]*nO2 + matA[1,3]*nCO2 = matA[0,0]*nH2 + matA[0,1]*nH2O + matA[0,2]*nO2 + matA[0,3]*nCO2
#    y = ( (matA[1,0]-matA[0,0])*x + (matA[1,2]-matA[0,2])*nO2 + (matA[1,3]-matA[0,3])*nCO2 ) / ( matA[0,1] - matA[1,1] ) 

#    matA[2,0]*nH2 + matA[2,1]*nH2O + matA[2,2]*nO2 + matA[2,3]*nCO2 = matA[0,0]*nH2 + matA[0,1]*nH2O + matA[0,2]*nO2 + matA[0,3]*nCO2
#    y = ( (matA[2,0]-matA[0,0])*x + (matA[2,2]-matA[0,2])*nO2 + (matA[2,3]-matA[0,3])*nCO2 ) / ( matA[0,1] - matA[2,1] ) 

#    EQUATIONS:
#
#    y = 1 - x - nCO2 - nO2         [1]
#    y = c1*x + c2*nO2 + c3*nCO2    [2]
#    y = d1*x + d2*nO2 + d3*nCO2    [3]
#    where:
#      c1 = (matA[1,0]-matA[0,0])/(matA[0,1]-matA[1,1])
#      c2 = (matA[1,2]-matA[0,2])/(matA[0,1]-matA[1,1])
#      c3 = (matA[1,3]-matA[0,3])/(matA[0,1]-matA[1,1])
#      d1 = (matA[2,0]-matA[0,0])/(matA[0,1]-matA[2,1])
#      d2 = (matA[2,2]-matA[0,2])/(matA[0,1]-matA[2,1])
#      d3 = (matA[2,3]-matA[0,3])/(matA[0,1]-matA[2,1])
#
#    From [2] and [3]:
#
#       d2*y = c1*d2*x + c2*d2*nO2 + c3*d2*nCO2
#    -) c2*y = c2*d1*x + c2*d2*nO2 + c2*d3*nCO2
#    -------------------------------------------
#     nCO2 = ((d2-c2)*y-(c1*d2-c2*d1)*x) / (c3*d2 - c2*d3)  [4]
#
#       d3*y = c1*d3*x + c2*d3*nO2 + c3*d3*nCO2
#    -) c3*y = c3*d1*x + c3*d2*nO2 + c3*d3*nCO2
#    -------------------------------------------
#      nO2 = ((d3-c3)*y-(c1*d3-c3*d1)*x)/(c2*d3 - c3*d2)    [5]
#
#    Substitute [4] and [5] to [1]:
#
#       y  = 1 - x - nCO2 - nO2
#          = 1 - x - ((d2-c2)*y-(c1*d2-c2*d1)*x)/(c3*d2-c2*d3) - ((d3-c3)*y-(c1*d3-c3*d1)*x)/(c2*d3-c3*d2)
#          = 1 - x - A1*y + A2*x - B1*y + B2*x
#    <=> y = (1 + (A2+B2-1)*x)/(A1+B1+1)
#    <=> x = (1 - (A1+B1+1)*y)/(1-(A2+B2))
#    where:
#      A1 = (d2-c2)/(c3*d2-c2*d3)
#      A2 = (c1*d2-c2*d1)/(c3*d2-c2*d3)
#      B1 = (d3-c3)/(c2*d3-c3*d2)
#      B2 = (c1*d3-c3*d1)/(c2*d3-c3*d2)
#

    c1 = (matA[1,0]-matA[0,0])/(matA[0,1]-matA[1,1])
    c2 = (matA[1,2]-matA[0,2])/(matA[0,1]-matA[1,1])
    c3 = (matA[1,3]-matA[0,3])/(matA[0,1]-matA[1,1])
    d1 = (matA[2,0]-matA[0,0])/(matA[0,1]-matA[2,1])
    d2 = (matA[2,2]-matA[0,2])/(matA[0,1]-matA[2,1])
    d3 = (matA[2,3]-matA[0,3])/(matA[0,1]-matA[2,1])
    A1 = (d2-c2)/(c3*d2-c2*d3)
    A2 = (c1*d2-c2*d1)/(c3*d2-c2*d3)
    B1 = (d3-c3)/(c2*d3-c3*d2)
    B2 = (c1*d3-c3*d1)/(c2*d3-c3*d2) 
    print "c2*d3-c3*d2", c2*d3-c3*d2

    print "test", (c1*xH2 + c2*xO2 + c3*xCO2- xH2O)
    print "test", (d1*xH2 + d2*xO2 + d3*xCO2- xH2O)


    print "nH2O = ", (A2+B2-1)/(A1+B1+1), " x nH2 + ", 1.0/(A1+B1+1)
    #nO2  = B1*y - B2*x = (B1*(A2+B2-1)/(A1+B1+1)-B2)*x + B1/(A1+B1+1)
    print "nO2  = ", (B1*(A2+B2-1)/(A1+B1+1)-B2), " x nH2 + ", B1/(A1+B1+1)
    #nCO2 = A1*y - A2*x = (A1*(A2+B2-1)/(A1+B1+1)-A2)*x + A1/(A1+B1+1)
    print "nCO2 = ", (A1*(A2+B2-1)/(A1+B1+1)-A2), " x nH2 + ", A1/(A1+B1+1)

    print ""
    print ""

#    x3 = np.logspace(-10, 0, 100)
    y3 = np.array(H2Ovec)
#    y3 = np.logspace(-10, 0, 100)
#    y3 = (1+(A2+B2-1)*x3)/(A1+B1+1)
    x3 = (1.0-(A1+B1+1)*y3)/(1-(A2+B2))
    x4 = (1.0-(1.0-(A1+B1+1)*y3)/(1-(A2+B2)))
#    plt.figure()
    plt.plot(y3,x3,'w', linewidth=2)

    plt.figure()
    plt.xlabel("nH2O")
    plt.ylabel("1-nH2")
    plt.loglog(y3,x4,'k', linewidth=2)

    print 'H2, H2O, O2, CO2'
    for ii in xrange(len(x3)):
        print 1-x3[ii], y3[ii], A1*y3[ii]-A2*x3[ii], B1*y3[ii]-B2*x3[ii]


##############################################################################
if __name__ == "__main__":

    dict_mol = read_cfgfile('data/gibbs_chemicals_5.cfg')
#    dict_mol = read_cfgfile('data/test4_chemicals.cfg')
#    dict_mol = read_cfgfile('chemicals.cfg')
    list_atm  = get_atmlist(dict_mol)


    solarabundance = np.loadtxt("data/solarabundance.txt", usecols=(0,2), dtype={'names':('ID','frac'),'formats':('i2','f8')})
    dict_solar = dict(solarabundance)
    matA = get_matA(dict_mol, list_atm, dict_solar, P_SURF, TEMP)

#    myplot(dict_mol, matA)
#    1/0


#    print "matA", matA
    cons = constraints(dict_mol, matA)
    bnds = [(0,1)]
    for ii in xrange(len(dict_mol)-1):
        bnds = bnds + [(0,1)]
    bnds = tuple(bnds)

    print "#altitude[km]", "pressure[mbar]"
    for ii in xrange(len(dict_mol)):
        print dict_mol[ii]['name']
    grid_P = np.logspace(np.log10(P_LIMIT), np.log10(P_SURF), P_NUM)
    #P = P_surf * exp(-z/H) => z = 
    for zi in xrange(P_NUM):

        array_n0 = np.zeros(len(dict_mol)) + 0.5
        array_n = minimize(gibbs_energy, array_n0, constraints=cons, bounds=bnds, args=(dict_mol, grid_P[zi], TEMP), method='SLSQP', options={'disp': False, 'maxiter':1000})
#    array_n = minimize(gibbs_energy, array_n0, constraints=cons, bounds=((0,1), (0,1), (0,1), (0,1)), args=(dict_mol, P_SURF, TEMP), method='SLSQP', options={'disp': True})
#    array_n = minimize_scalar(gibbs_energy, array_n0, constraints=cons, jac=gibbs_energy_deriv, args=(dict_mol, P_SURF, TEMP), method='L-BFGS-B', options={'disp': True})

        print "----------------------"
        print SH*np.log(P_SURF/grid_P[zi])*1e-5, grid_P[zi], array_n['x']




#    print "----------------------"
#    print "INITIAL:"
#    for ii in xrange(len(dict_mol)):
#        print dict_mol[ii]['name'], '\t', array_n0[ii]
#    print "----------------------"
#    print "FINAL:"
#    for ii in xrange(len(dict_mol)):
#        print "%s\t%e" % (dict_mol[ii]['name'], array_n['x'][ii])
#    print "----------------------"
