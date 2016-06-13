#!/usr/bin/env python
from molecules import molecules
import numpy as np

##############################################################################
# ****************************************************************************
# CAUTION!! This is not originaly. This is based on lbl2od
# ****************************************************************************
##############################################################################

def extract_hitran ( dataFile, wn_limit, molecule, isoNr, strMin=0.0 ):

    # read all lines in interval and return dictionary
    try:             molNr = molecules[molecule]['hitran']
    except KeyError: raise SystemExit, 'ERROR: invalid/unknown molecule ' + repr(Mol)

    try:
        hitran = open ( dataFile )
    except IOError:
        raise SystemExit, 'opening Hitran data file "' + dataFile + '" failed!'

    xBegin, xHigh = wn_limit

    if not (isoNr or strMin):
        lines = extract_Mol (hitran, xBegin,xHigh, molNr)
    elif isoNr==0 and strMin>0.0:
        lines = extract_MolStr (hitran, xBegin,xHigh, molNr, strMin)
    elif strMin<=0.0:
        lines = extract_MolIso (hitran, xBegin,xHigh, molNr, isoNr)
    else:
        lines = extract_MolIsoStr (hitran, xBegin,xHigh, molNr, isoNr, strMin)
        hitran.close()

    linedata = core_parameters ( lines, dataFile )

    print "linedata", linedata
    return linedata
            #return molecules, isotopes, positions, strengths, energies, airWidths, selfWidths, tempExps

##########################################################################################################################################

def bisect_first_line (db, xBegin, xEnd, iw=3,lw=15):
    """ Search and return first line of Hitran or Geisa formatted database. """
    # skip optional header (comment) section
    if 'hit' in db.name.lower():
        # hitran database, skip header section with comment lines (indicated with mol=00)
        mol  = -1
        try:
            while mol<1:  record=db.readline(); mol = int(record[:2])
        except ValueError, msg:
            if len(record.rstrip()) not in (100,160):
                print 'Hitran type database, but length of first data record is not 100 or 160'
            raise SystemExit, str(msg) + '\nERROR reading header section of Hitran-type database (trying to parse molec id number)'
    else:
        # hitran database, skip header section with comment lines (indicated with mol=00)
        record=db.readline()

    # first data record  with spectroscopic line transition in file:  parse it carefully!
    lenRec   = len(record)
    locFirst = db.tell() - lenRec
    try:
        recFirst = record
        xFirst   = float(recFirst[iw:lw])
    except ValueError, msg:
        if 'hit' in db.name.lower() and not len(record.rstrip()) in (100,160):
            print 'Hitran type database, but length of first data records is not 100 or 160'
        elif 'geisa' in db.name.lower() and len(record.rstrip()) in (100,160):
            print 'Geisa type database, but length of first data records is 100 or 160 (looks like Hitran)'
        raise SystemExit, 'ERROR reading first data record of database (trying to parse wavenumber)\n' + str(msg)

    # move to last record in file
    db.seek(-lenRec,2)
    locLast = db.tell()
    recLast = db.readline()
    xLast   = float(recLast[iw:lw])
    # check spectral range
    if  xFirst <= xBegin < xEnd <= xLast:
        pass
    elif  xFirst>xEnd or xLast<xBegin:
        raise SystemExit, 'requested spectral range not in database ' + str(xFirst) + ' --- ' + str(xLast)
    else:
        if xFirst>1.0 or xLast<20000.0:  # apparently a subset of the full hitran database
            print 'WARNING:  requested spectral range only partly in database ' + str(xFirst) + ' --- ' + str(xLast)
        # record number of very first and last lines (locFirst, locLast are byte numbers)
    lineFirst = locFirst/lenRec
    lineLast  = (locLast-locFirst)/lenRec

    # bisecting to desired first line
    while lineLast-lineFirst > 1:
        mid=(lineLast+lineFirst)/2
        db.seek(locFirst+mid*lenRec)
        recMid=db.readline()
        xMid = float(recMid[iw:lw])
        if   xBegin<xMid: lineLast =mid; xLast=xMid
        elif xBegin>xMid: lineFirst=mid; xFirst=xMid
        else:
            # backstep: there are possibly several lines at exactly this position
            while 1:
                db.seek(-2*lenRec,1); rec=db.readline(); xMid=float(rec[iw:lw])
                if xMid<xBegin:
#                    print "# first line in spectral range at record ", mid, "found in", clock(), "sec\n", rec[:67]
                    return        db.readline()
        
    if    xMid<xBegin:      record = db.readline()
    else:                   record = recMid
#    print "# first line in spectral range at record number", mid, "found in", clock(), "sec\n", record[:67]

    return record

##########################################################################################################################################

def extract_All (hitran, xBegin, xHigh):
    """ Read all lines up to an upper wavenumber limit from Hitran formatted database. """
    # proceed to first requested line
    record = bisect_first_line (hitran, xBegin, xHigh)
    # initialize list if lines
    lines = []
    # collect lines
    while record:
        mol    = int(record[:2])
        if mol>0:
            wvn = float(record[3:15])
            if wvn<=xHigh:  lines.append(record)
            else:           break
# read next record
        record = hitran.readline()

#    if len(lines)>0:  print '# last  line     accepted \n', lines[-1][:67]
#    if record:        print '# first line not accepted \n', record[:67]     # empty string returned at end-of-file

    return lines

##########################################################################################################################################

def extract_Mol (hitran, xBegin,xHigh, getMol):
    """ Read lines of a given molecule up to an upper wavenumber limit from Hitran formatted database. """
    # proceed to first requested line
    record = bisect_first_line (hitran, xBegin, xHigh)
    # initialize list if lines
    lines = []
    # collect lines
    while record:
        mol = int(record[:2])
        wvn = float(record[3:15])
        if wvn>xHigh: break
        if mol==getMol: lines.append(record)
# read next record
        record = hitran.readline()

    if len(lines)>0:  print '# last  line     accepted \n', lines[-1][:67]
    if record:        print '# first line not accepted \n', record[:67]     # empty string returned at end-of-file

    return lines

##########################################################################################################################################

def extract_MolIso (hitran, xBegin,xHigh, getMol, getIso):
    """ Read lines of a given molecule/isotope up to an upper wavenumber limit from Hitran formatted database. """
    # proceed to first requested line
    record = bisect_first_line (hitran, xBegin, xHigh)
    # initialize list if lines
    lines = []
    # collect lines
    while record:
        mol = int(record[:2])
        iso = int(record[2:3])
        wvn = float(record[3:15])
        if wvn>xHigh: break
        if mol==getMol and iso==getIso: lines.append(record)
        # read next record
        record = hitran.readline()

#    if len(lines)>0:  print '# last  line     accepted \n', lines[-1][:67]
#    if record:        print '# first line not accepted \n', record[:67]     # empty string returned at end-of-file

    return lines

##########################################################################################################################################

def extract_MolStr (hitran, xBegin,xHigh, getMol, strMin):
    """ Read strong lines of a given molecule up to an upper wavenumber limit from Hitran formatted database. """
    # proceed to first requested line
    record = bisect_first_line (hitran, xBegin, xHigh)
    # initialize list if lines
    lines = []
    # collect lines
    while record:
        mol = int(record[:2])
        wvn = float(record[3:15])
        str = float(record[15:25])
        if wvn>xHigh: break
        if mol==getMol and str>=strMin: lines.append(record)
# read next record
        record = hitran.readline()

#    if len(lines)>0:  print '# last  line     accepted \n', lines[-1][:67]
#    if record:        print '# first line not accepted \n', record[:67]     # empty string returned at end-of-file

    return lines

##########################################################################################################################################

def extract_MolIsoStr (hitran, xBegin,xHigh, getMol, getIso, strMin):
    """ Read strong lines of a given molecule/isotope up to an upper wavenumber limit from Hitran formatted database. """
    # proceed to first requested line
    record = bisect_first_line (hitran, xBegin, xHigh)
    # initialize list if lines
    lines = []
    # collect lines
    while record:
        mol = int(record[:2])
        iso = int(record[2:3])
        str = float(record[15:25])
        if wvn>xHigh: break
        if mol==getMol and iso==getIso and str>=strMin:  lines.append(record)
# read next record
        record = hitran.readline()

#    if len(lines)>0:  print '# last  line     accepted \n', lines[-1][:67]
#    if record:        print '# first line not accepted \n', record[:67]     # empty string returned at end-of-file

    return lines




def core_parameters (lines, dataFile):
    """ Given a list of data base records (one entry per transition) return python lists of the most important spectrocopic line parameters. """
    # column start/stop for all types of parameters
    if 'hit' in dataFile.lower() or 'sao' in dataFile.lower():
        iw,lw, iS,lS, iE, lE, iA,lA, isw, lsw, iT,lT, iI,lI = 3,15, 15,25, 45,55, 35,40, 40,45, 55,59, 2,3
    else:
        from geisa import set_geisa_fields
        fields = set_geisa_fields (dataFile)
        for key,val in fields.items(): exec key + '=' + repr(val)

    # now extract columns, convert to appropriate type (int/float)
    linedata = {}
    linedata['position'] = np.array([float(line[iw:lw]) for line in lines])
#    linedata['strength'] = [float(replace(line[iS:lS],'D','e')) for line in lines]
    linedata['strength'] = np.array([float(line[iS:lS].replace('D','e')) for line in lines])
    linedata['energy']   = np.array([float(line[iE:lE]) for line in lines])
    linedata['airWidth'] = np.array([float(line[iA:lA]) for line in lines])
    linedata['Tdep']     = np.array([float(line[iT:lT]) for line in lines])

#    isoNr     = [int(line[iI:lI]) for line in lines]
    if 'hit' in dataFile.lower() or 'geisa' in dataFile.lower():
        linedata['selfWidth'] = np.array([float(line[isw:lsw]) for line in lines])
    else:
        linedata['selfWidth'] = len(lines)*np.array([0]) # a list with nLines zeros
    return linedata
