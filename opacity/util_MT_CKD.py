import MT_CKD
import numpy as np


#vi = np.array(5050)
#csh2or = np.array(5050)
#cfh2or = np.array(5050)

#MT_CKD.get_cntnm(1013.0, 296.0, vi, csh2or, cfh2or)


print "Pressure? (mbar)"
pres = input()
print ""
print "Temperature? (K)"
temp = input()


vi, csh2or, cfh2or = MT_CKD.get_cntnm(pres, temp)


for ii in xrange(1000):
    print vi[ii], csh2or[ii], cfh2or[ii]
#print vi[1], csh2or[1], cfh2or[1]
