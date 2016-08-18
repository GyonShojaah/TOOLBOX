#INFILE_TAG = "out/transmission/prof_ANN1962-1965aijlE_g10_R10_P10_tl_P100days_cld_trapz-1000_line"
#INFILE_TAG = "out/transmission/profile_USstandard_CH4_N2O_00010-10000_m10001_quad-2000_line"
#INFILE_TAG = "out/transmission/segura2005_Temp_H2O_CH4"
#UNIT   = "wl"
UNIT  = "wn"

#from setting import *
import sys
from copy import deepcopy
import numpy as np

INFILE_TAG = sys.argv[1]
RESOL  = int(sys.argv[2])


#=============================================================================
# Module
#=============================================================================

def lower_resolution(wn_array, sp_array, resolution):

    count = 0

    wn_new_array = np.zeros(1)
    sp_new_array = np.zeros(1)

    wn_tmp = wn_array[0]
    jj = 0
    for ii in range(len(wn_array)):

        if wn_array[ii] > wn_tmp + wn_tmp/resolution :
            wn_new_array[jj] = wn_new_array[jj]/(count*1.0)
            sp_new_array[jj] = sp_new_array[jj]/(count*1.0)
            wn_new_array = np.r_[wn_new_array, np.zeros(1)]
            sp_new_array = np.r_[sp_new_array, np.zeros(1)]
            wn_tmp = wn_tmp + wn_tmp/resolution
            count = 0
            jj += 1

        wn_new_array[jj] += wn_array[ii]
        sp_new_array[jj] += sp_array[ii]
        count += 1

    wn_new_array[jj] = wn_new_array[jj]/(count*1.0)
    sp_new_array[jj] = sp_new_array[jj]/(count*1.0)

    return wn_new_array, sp_new_array

def readfile(filename, key=0) :
    """
    Read cross section look-up table.
    """
    data = np.load(filename)
    if key :
        return data[key]
    else :
        return data


def myplot_func(filename, myfunc, x_min, x_max, xlog=False, ylog=False, str_title="", str_xlabel="x", str_ylabel="y", int_num=10) :
    """
    plot function
    """

    plt.title(str_title)
    plt.xlabel(str_xlabel, fontsize=20)
    plt.ylabel(str_ylabel, fontsize=20)
    plt.grid()
    x = np.linspace(x_min, x_max, num=int_num)
    if xlog : 
        x = np.logspace(np.log10(x_min), np.log10(x_max), num=int_num)
        plt.xscale('log')
    if ylog :
        plt.yscale('log')

    y = np.zeros(int_num)
    for i in range(int_num) :
        y[i] = myfunc(x[i])

    plt.tick_params(axis='both')
    plt.plot(x,y)
    plt.savefig(filename)


def myplot_tbl(filename, x_array, y_array, xlog=False, ylog=False, str_title="", str_xlabel="x", str_ylabel="y") :
    """
    plot tabulated data
    """

    plt.title(str_title)
    plt.xlabel(str_xlabel, fontsize=20)
    plt.ylabel(str_ylabel, fontsize=20)
    plt.grid()
    if xlog : 
        plt.xscale('log')
    if ylog :
        plt.yscale('log')

    plt.tick_params(axis='both')
    plt.plot(x_array,y_array)
    plt.savefig(filename)



#=============================================================================
# main
#=============================================================================

if __name__ == "__main__":

   data = np.loadtxt(INFILE_TAG + ".txt")
   sp_array = data.T[1]
   if (UNIT == "wn"):
       wn_array = data.T[0]
   elif (UNIT == "wl"):
       wn_array = 1.0e4/data.T[0]
   else:
       print "INVALID UNIT"
       sys.exit()


   wn_new_array, sp_new_array = lower_resolution(wn_array, sp_array, RESOL)

   if (UNIT == "wn"):
       data2 = np.dstack([wn_new_array, sp_new_array])
   else:
       data2 = np.dstack([1.0e4/wn_new_array, sp_new_array])

#   print data2
   np.savetxt(INFILE_TAG+"_R"+str(RESOL)+".txt", data2[0])


