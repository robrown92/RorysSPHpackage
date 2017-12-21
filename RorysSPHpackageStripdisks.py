# -*- coding: utf-8 -*-
"""
Created on Wed Feb 17 14:33:00 2016

@author: rob1g10
"""
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import pylab as plt
import os
from scipy.optimize import leastsq

"""----------------------------------------------------------------------------

                    Strip disk for info

----------------------------------------------------------------------------"""

def stripdisk(filename, beforePhrases, numOfValues, timefreq, tol=0.01):
    """
    Goes through a disk file and extracts data.

    INPUT VARIABLES

    filename (str) - name of disk file to extract data from.

    beforePhrases (array of str) - the phrases that come before the desired
                                   data (are obviously best to be as unique as
                                   possible)

    numOfValues (int) - number of values that the user wishes to extract

    timefreq (float) - how often data is to be taken (use timefreq='orbit' for
                        every orbit).

    tol (float) - tolerance for code to allow for matching of desired orbits

    OUTPUT VARIABLES

    values - an array giving the extracted data

    """

    " Open file "
    f = open(filename, 'r')

    timePhraseCount = 0
    findTimeCount = 0

    """ Want to find timeStart from disk file """

    findTimePhrases = ["TIME", ":"]

    for line in f:
        strippedLine = line.strip()
        columns = strippedLine.split()
        for i in range(len(columns)):
            if(timePhraseCount == 2):
                timeStart = float(columns[i])
                timePhraseCount = 0
                findTimeCount += 1
            elif (columns[i] == findTimePhrases[timePhraseCount]):
                timePhraseCount += 1
        if(findTimeCount == 1):
            break

    f.close()

    """ Want to find unit of time for simulation """

    f = open(filename, 'r')

    utimePhraseCount = 0
    finduTimeCount = 0

    """ Want to find timeStart from disk file """

    finduTimePhrases = ["time", ":"]

    for line in f:
        strippedLine = line.strip()
        columns = strippedLine.split()
        for i in range(len(columns)):
            if(utimePhraseCount == 2):
                utime = float(columns[i])
                utimePhraseCount = 0
                finduTimeCount += 1
            elif (columns[i] == finduTimePhrases[utimePhraseCount]):
                utimePhraseCount += 1
        if(finduTimeCount == 1):
            break

    f.close()

    f = open(filename, 'r')

    values = np.array([])
    times = np.array([])

    beforePhraseslen = len(beforePhrases)

    count = 0
    phraseCount = 0
    orbitCount = 0
    timeCount = 0
    sigfig = len(str(timeStart))+1
#    """ Without orbit checking """
#    if(timeStart == None):
#        for line in f:
#            strippedLine = line.strip()
#            columns = strippedLine.split()
#            for i in range(len(columns)):
#                if(phraseCount == beforePhraseslen):
#                    values = np.append(values, float(columns[i]))
#                    phraseCount = 0
#                    orbitCount = 0
#                    count += 1
#                elif (columns[i] == beforePhrases[phraseCount]):
#                    phraseCount += 1
#            if(count == numOfValues):
#                break

    if(timefreq == 'orbit'):
        """ With orbit checking """
        for line in f:
            strippedLine = line.strip()
            columns = strippedLine.split()
            for i in range(len(columns)):
                orbitCheck = orbitCount + timeStart
                orbitCheck = "{:.10e}".format(orbitCheck)
                if(columns[i] == 'orbital'
                   and columns[i+1] == 'period'
                   and columns[i+2] == ':'
                   and abs(float(columns[i+3][0:sigfig]) - float(orbitCheck[0:sigfig])) < tol):
                    orbitCount +=1
                if(count < orbitCount):
                    if(phraseCount == beforePhraseslen):
                        values = np.append(values, float(columns[i]))
                        phraseCount = 0
                        count += 1
                    elif (columns[i] == beforePhrases[phraseCount]):
                        phraseCount += 1

            if(count == numOfValues):
                break

    elif(timefreq == 'day'):
        """ With orbit checking """
        timefreq = (24 * 3600) / utime
        for line in f:
            strippedLine = line.strip()
            columns = strippedLine.split()
            for i in range(len(columns)):
                timeCheck = timeCount*timefreq + timeStart
                timeCheck = "{:.10e}".format(timeCheck)
                if(columns[i] == 'TIME'
                   and columns[i+1] == ':'
                   and abs(float(columns[i+2][0:sigfig]) - float(timeCheck[0:sigfig])) < tol):
                    timeCount +=1
                if(count < timeCount):
                    if(phraseCount == beforePhraseslen):
                        values = np.append(values, float(columns[i]))
                        times = np.append(times, timeCheck)
                        phraseCount = 0
                        count += 1
                    elif (columns[i] == beforePhrases[phraseCount]):
                        phraseCount += 1

            if(count == numOfValues):
                break

    else:
        """ With timefreq checking """
        for line in f:
            strippedLine = line.strip()
            columns = strippedLine.split()
            for i in range(len(columns)):
                timeCheck = timeCount*timefreq + timeStart
                timeCheck = "{:.10e}".format(timeCheck)
                if(columns[i] == 'TIME'
                   and columns[i+1] == ':'
                   and abs(float(columns[i+2][0:sigfig]) - float(timeCheck[0:sigfig])) < tol):
                    timeCount +=1
                if(count < timeCount):
                    if(phraseCount == beforePhraseslen):
                        values = np.append(values, float(columns[i]))
                        phraseCount = 0
                        count += 1
                    elif (columns[i] == beforePhrases[phraseCount]):
                        phraseCount += 1

            if(count == numOfValues):
                break

    f.close()

    return times, values
    
    
def stripdiskForDISKs(filename, beforePhrases, DISKstart, DISKend):
    """
    Goes through a disk file and extracts data.

    INPUT VARIABLES

    filename (str) - name of disk file to extract data from.

    beforePhrases (array of str) - the phrases that come before the desired
                                   data (are obviously best to be as unique as
                                   possible)

    DISKstart (int) - the DISK file to start extraction at.
    
    DISKend (int) - the DISK file to end extraction at.

    tol (float) - tolerance for code to allow for matching of desired orbits

    OUTPUT VARIABLES

    values - an array giving the extracted data

    """

    f = open(filename, 'r')

    values = np.array([])

    beforePhraseslen = len(beforePhrases)

    count = 0
    started = 0
    phraseCount = 0
    onDISK = 0
    BREAK = False

    for line in f:
        if(BREAK):
            break
        strippedLine = line.strip()
        columns = strippedLine.split()
        for i in range(len(columns)):
            if(columns[i][0:4] == 'DISK' and int(columns[i][4:]) >= DISKstart):
                started = 1
            if(columns[i] == 'W' and columns[i+1] == 'R' and columns[i+2] == 'I' 
               and columns[i+3] == 'T' and columns[i+4] == 'T' and columns[i+5] == 'E'
               and columns[i+6] == 'N' and started == 1):
                onDISK = 1
            if(onDISK == 1):
                if(phraseCount == beforePhraseslen):
                    values = np.append(values, float(columns[i]))
                    phraseCount = 0
                    onDISK = 0
                elif (columns[i] == beforePhrases[phraseCount]):
                    phraseCount += 1
#            print "DISK{0:03}".format(DISKend)
            if(columns[i][0:4] == 'DISK' and int(columns[i][4:]) >= DISKend):
                BREAK = True
                

    f.close()

    return values

"""----------------------------------------------------------------------------

                               Strip disks

----------------------------------------------------------------------------"""

""" ACTIVE PARTS -----------------------------------------------------------"""

#numActiveParts_a1 = stripdisk('/data/rob1g10/SPH/DATA/DISKS/injection_off_a1.5e11/disk006',
#                           ['number', 'of', 'active', 'particles', ':'], None, orbitStart=108.8, tol=0.01)
#
#numActiveParts_justBe1 = stripdisk('/data/rob1g10/SPH/DATA/DISKS/injection_off_justBe/disk009',
#                           ['number', 'of', 'active', 'particles', ':'], None, orbitStart=60.58, tol=0.01)
#numActiveParts_justBe2 = stripdisk('/data/rob1g10/SPH/DATA/DISKS/injection_off_justBe/disk010',
#                           ['number', 'of', 'active', 'particles', ':'], None, orbitStart=85.58, tol=0.01)
#numActiveParts_justBe = np.append(numActiveParts_justBe1, numActiveParts_justBe2)
#
#numActiveParts_JJ1 = stripdisk('/data/rob1g10/SPH/DATA/DISKS/injection_off_JJrun/disk003',
#                           ['number', 'of', 'active', 'particles', ':'], None, orbitStart=30.02, tol=0.01)
#numActiveParts_JJ2 = stripdisk('/data/rob1g10/SPH/DATA/DISKS/injection_off_JJrun/disk004',
#                           ['number', 'of', 'active', 'particles', ':'], None, orbitStart=37.02, tol=0.01)
#numActiveParts_JJ3 = stripdisk('/data/rob1g10/SPH/DATA/DISKS/injection_off_JJrun/disk005',
#                           ['number', 'of', 'active', 'particles', ':'], None, orbitStart=39.02, tol=0.01)
#numActiveParts_JJ = np.append(numActiveParts_JJ1, numActiveParts_JJ2)
#numActiveParts_JJ = np.append(numActiveParts_JJ, numActiveParts_JJ3)
#
#plt.plot(range(0, len(numActiveParts_a1)), numActiveParts_a1 / max(numActiveParts_a1), label='~80 day orbit 0.1e')
#plt.plot(np.append(range(0, len(numActiveParts_justBe1)), range(23, len(numActiveParts_justBe2)+23)),
#         numActiveParts_justBe / max(numActiveParts_justBe), label='just Be star')
#plt.plot(range(0, len(numActiveParts_JJ)), numActiveParts_JJ / max(numActiveParts_JJ), label='~24 day orbit 0.3e')
#plt.legend()

""" ACCRETED PARTS ---------------------------------------------------------"""


def createAccretionArrayForDISKs(filenames, SOURCE, DISKstart, DISKend):
    """
    Creates an array of accreted particles over each unit of time.

    INPUT VARIABLES

    filenames (list of str) - name of disk files to extract data from.

    SOURCE (str) - can be NS, Be or KILL - defines what source of loss of
                   particles from system the function will search for.
                   NEW - numparts finds number of active particles in
                   simulation.

    timefreq (float) - how often data is to be taken(use timefreq='orbit' for
                       every orbit).

    DISKstart (int) - the DISK file to start extraction at.
    
    DISKend (int) - the DISK file to end extraction at.

    tol (float) - tolerance for code to allow for matching of desired orbits

    OUTPUT VARIABLES

    values - an array giving the extracted data

    """

    """ Set up phrases for search """

    if(SOURCE == "NS"):
        PHRASES = ['point', 'masses', ':', '1', ':', '2', ':']
    elif(SOURCE == "Be"):
        PHRASES = ['point', 'masses', ':', '1', ':']
    elif(SOURCE == "KILL"):
        PHRASES = ['number', 'of', 'killed', 'particles', ':']
    elif(SOURCE == 'numparts'):
        PHRASES = ['number', 'of', 'active', 'particles', ':']
    elif(SOURCE == 'reassign'):
        PHRASES = ['number', 'of', 'reassigned', 'part.', ':']

    """ Strip disks """

    accreteArray = np.array([])
    breaks = np.array([0])

    for i in range(len(filenames)):
        #timeStart += breaks[i] * timefreq
        times, accreteArray = np.append(accreteArray, stripdiskForDISKs(filenames[i],
                             PHRASES, DISKstart, DISKend))
        breaks = np.append(breaks, len(accreteArray))

    """ Change from cumulative accretion to per unit time """
    orbitalAccreteArray = np.zeros_like(accreteArray)
    #orbitalAccreteArray[-1] = accreteArray[-1] - accreteArray[-2]

    for j in range(1, len(breaks)-1):
        orbitalAccreteArray[breaks[j]] = accreteArray[breaks[j]]

    for j in range(len(breaks)-1, 0, -1):
        for i in range(breaks[-j]-1, breaks[-j-1], -1):
            orbitalAccreteArray[i] = int(accreteArray[i]) - int(accreteArray[i-1])

    return times, orbitalAccreteArray
    

def createAccretionArray(filenames, SOURCE, numOfValues, timefreq, timeStart=None, tol=0.01):
    """
    Creates an array of accreted particles over each unit of time.

    INPUT VARIABLES

    filenames (list of str) - name of disk files to extract data from.

    SOURCE (str) - can be NS, Be or KILL - defines what source of loss of
                   particles from system the function will search for.
                   NEW - numparts finds number of active particles in
                   simulation.

    numOfValues (int) - number of values that the user wishes to extract

    timefreq (float) - how often data is to be taken(use timefreq='orbit' for
                       every orbit).

    timeStart (array of float) - used to extract data only every period, this
                                 chooses at what point during the period it
                                 will start. (if not specified, it will not
                                 use this feature) should be an array of
                                 length filenames.

    tol (float) - tolerance for code to allow for matching of desired orbits

    OUTPUT VARIABLES

    values - an array giving the extracted data

    """

    """ Set up phrases for search """

    if(SOURCE == "NS"):
        PHRASES = ['point', 'masses', ':', '1', ':', '2', ':']
    elif(SOURCE == "NScapture"):
        PHRASES = ['number', 'of', 'part.', 'captured', 'by', 'point', 'mass', '2', ':']
    elif(SOURCE == "Be"):
        PHRASES = ['point', 'masses', ':', '1', ':']
    elif(SOURCE == "KILL"):
        PHRASES = ['number', 'of', 'killed', 'particles', ':']
    elif(SOURCE == 'numparts'):
        PHRASES = ['number', 'of', 'active', 'particles', ':']
    elif(SOURCE == 'reassign'):
        PHRASES = ['number', 'of', 'reassigned', 'part.', ':']

    """ Strip disks """

    accreteArray = np.array([])
    times = np.array([])
    breaks = np.array([0])

    for i in range(len(filenames)):
        #timeStart += breaks[i] * timefreq
        timesAppend, accreteArrayAppend = stripdisk(filenames[i],
                             PHRASES, numOfValues, timefreq, tol)
        accreteArray = np.append(accreteArray, accreteArrayAppend)
        times = np.append(times, timesAppend)
        breaks = np.append(breaks, len(accreteArray))

    """ Change from cumulative accretion to per unit time """
    orbitalAccreteArray = np.zeros_like(accreteArray)
    #orbitalAccreteArray[-1] = accreteArray[-1] - accreteArray[-2]

    for j in range(1, len(breaks)-1):
        orbitalAccreteArray[breaks[j]] = accreteArray[breaks[j]]

    for j in range(len(breaks)-1, 0, -1):
        for i in range(breaks[-j]-1, breaks[-j-1], -1):
            orbitalAccreteArray[i] = int(accreteArray[i]) - int(accreteArray[i-1])

    return times, orbitalAccreteArray

def findCumulativeNumParts(filenames, numOfValues, timefreq, timeStart=None, tol=0.01):
    """
    Creates a cumulative array of particles.

    INPUT VARIABLES

    filenames (list of str) - name of disk files to extract data from.

    numOfValues (int) - number of values that the user wishes to extract

    timefreq (float) - how often data is to be taken(use timefreq='orbit' for
                       every orbit).

    timeStart (array of float) - used to extract data only every period, this
                                 chooses at what point during the period it
                                 will start. (if not specified, it will not
                                 use this feature) should be an array of
                                 length filenames.

    tol (float) - tolerance for code to allow for matching of desired orbits

    OUTPUT VARIABLES

    values - an array giving the extracted data

    """

    """ Set up phrases for search """

    PHRASES = ['number', 'of', 'active', 'particles', ':']

    """ Strip disks """

    accreteArray = np.array([])
    breaks = np.array([0])

    for i in range(len(filenames)):
        accreteArray = np.append(accreteArray, stripdisk(filenames[i],
                             PHRASES, numOfValues, timefreq, tol))
        breaks = np.append(breaks, len(accreteArray))

    return accreteArray

"""----------------24.3d,83.04d---------------------------------------------"""
"""------------------------------NSaccretion--------------------------------"""

#accretionArray = createAccretionArray(['/data/rob1g10/SPH/DATA/DISKS/injection_off_JJrun/disk003',
#                                       '/data/rob1g10/SPH/DATA/DISKS/injection_off_JJrun/disk004',
#                                       '/data/rob1g10/SPH/DATA/DISKS/injection_off_JJrun/disk005'],
#                             "NS", None, 'day', tol=0.01)
#accretionArray2 = createAccretionArray(['/data/rob1g10/SPH/DATA/DISKS/injection_off_a1.5e11/disk006'],
#                             "NS", None, 'day', tol=0.01)
##plt.plot(range(500), accretionArray2[0:500], label='83d, 0.1e')
#plt.plot(range(len(accretionArray)), accretionArray, label='24.3d, 0.3e')
#plt.legend()
#plt.title('NS accretion over ~1800 days after injection is shut off - 24.3d')
#plt.xlabel('Days')
#plt.ylabel('Number of particles accreted by neutron star')
#plt.show()

"""-----------------------numActiveParts------------------------------------"""

#numpartsArrayjustBe = findCumulativeNumParts(['/data/rob1g10/SPH/DATA/DISKS/injection_off_justBe/disk006',
#                                       '/data/rob1g10/SPH/DATA/DISKS/injection_off_justBe/disk007',
#                                       '/data/rob1g10/SPH/DATA/DISKS/injection_off_justBe/disk008'],
#                                       None, 24 * 3600 / 2.8344e5, timeStart=[254.05, 382.05, 668.05], tol=0.01)
#numpartsArray24d = findCumulativeNumParts(['/data/rob1g10/SPH/DATA/DISKS/injection_off_JJrun/disk004',
#                                       '/data/rob1g10/SPH/DATA/DISKS/injection_off_JJrun/disk005'],
#                                       None, 24 * 3600 / 2.7502e5, timeStart=[226.5, 242.225], tol=0.01)
#numpartsArray83d = findCumulativeNumParts(['/data/rob1g10/SPH/DATA/DISKS/injection_off_a1.5e11/disk006'],
#                                       None, 24 * 3600 / 1.142e6, timeStart=[684.05], tol=0.01)
#plt.plot(range(1400), numpartsArray83d[0:1400], label='83d, 0.1e')
#plt.plot(range(1400), numpartsArray24d[0:1400]/max(numpartsArray24d), label='24.3d, 0.3e')
#plt.plot(range(1400), numpartsArrayjustBe[0:1400], label='just Be')
#plt.legend()
#plt.title('Fraction of particles remaining after injection is shut off - 24.3d')
#plt.xlabel('Days')
#plt.ylabel('Number of particles divided by number of particles at shut off')
#plt.show()

"""--------------------------------40d--------------------------------------"""
"""------------------------------NSaccretion--------------------------------"""
"""-----------------------injection_off-------------------------------------"""

#NSArray0e = createAccretionArray(['/data/rob1g10/SPH/DATA/DISKS/40d_0.0e/disk005',
#                                       '/data/rob1g10/SPH/DATA/DISKS/40d_0.0e/disk006'],
#                                      "NS", None, 24 * 3600 / 5.5004e5, timeStart=[286.25, 494.85], tol=0.01)
#NSArray2e = createAccretionArray(['/data/rob1g10/SPH/DATA/DISKS/40d_0.2e/disk005',
#                                       '/data/rob1g10/SPH/DATA/DISKS/40d_0.2e/disk006'],
#                                      "NS",None, 24 * 3600 / 5.5004e5, timeStart=[290.25, 524.05], tol=0.01)
#NSArray4e = createAccretionArray(['/data/rob1g10/SPH/DATA/DISKS/40d_0.4e/disk005',
#                                       '/data/rob1g10/SPH/DATA/DISKS/40d_0.4e/disk006'],
#                                      "NS", None, 24 * 3600 / 5.5004e5, timeStart=[296.25, 610.25], tol=0.01)
#plt.plot(range(3800), NSArray0e[0:3800], label='40d 0.0e')
#plt.plot(range(3800), NSArray2e[0:3800], label='40d 0.2e')
#plt.plot(range(3800), NSArray4e[0:3800], label='40d 0.4e')
#plt.legend()
#plt.title('NS accretion for 40d orbits after injection is shut off')
#plt.xlabel('Days')
#plt.ylabel('Number of particles after injection is shut off')
#plt.show()

"""-----------------------fromStart-----------------------------------------"""

#NSArray0e = createAccretionArray(['/data/rob1g10/SPH/DATA/DISKS/40d_0.0e/disk001',
#                                       '/data/rob1g10/SPH/DATA/DISKS/40d_0.0e/disk002',
#                                       '/data/rob1g10/SPH/DATA/DISKS/40d_0.0e/disk003'],
#                                      "NS", None, 24 * 3600 / 5.5004e5, timeStart=[0.05, 120.65, 188.05], tol=0.01)
#NSArray2e = createAccretionArray(['/data/rob1g10/SPH/DATA/DISKS/40d_0.2e/disk001',
#                                       '/data/rob1g10/SPH/DATA/DISKS/40d_0.2e/disk002',
#                                       '/data/rob1g10/SPH/DATA/DISKS/40d_0.2e/disk003'],
#                                      "NS",None, 24 * 3600 / 5.5004e5, timeStart=[0.05, 118.05, 187.05], tol=0.01)
#NSArray4e = createAccretionArray(['/data/rob1g10/SPH/DATA/DISKS/40d_0.4e/disk001',
#                                       '/data/rob1g10/SPH/DATA/DISKS/40d_0.4e/disk002',
#                                       '/data/rob1g10/SPH/DATA/DISKS/40d_0.4e/disk003'],
#                                      "NS", None, 24 * 3600 / 5.5004e5, timeStart=[0.05, 120.45, 191.45], tol=0.01)
#plt.plot(range(700), NSArray0e[0:700], label='40d 0.0e')
#plt.plot(range(700), NSArray2e[0:700], label='40d 0.2e')
#plt.plot(range(700), NSArray4e[0:700], label='40d 0.4e')
#plt.legend()
#plt.title('NS accretion for 40d orbits after injection is shut off')
#plt.xlabel('Days')
#plt.ylabel('Number of particles after injection is shut off')
#plt.show()

"""-----------------------numActiveParts------------------------------------"""
"""-----------------------injection_off-------------------------------------"""

#numpartsArray0e = findCumulativeNumParts(['/data/rob1g10/SPH/DATA/DISKS/40d_0.0e/disk005',
#                                       '/data/rob1g10/SPH/DATA/DISKS/40d_0.0e/disk006'],
#                                       None, 24 * 3600 / 5.5004e5, timeStart=[286.25, 494.85], tol=0.01)
#numpartsArray2e = findCumulativeNumParts(['/data/rob1g10/SPH/DATA/DISKS/40d_0.2e/disk005',
#                                       '/data/rob1g10/SPH/DATA/DISKS/40d_0.2e/disk006'],
#                                       None, 24 * 3600 / 5.5004e5, timeStart=[290.25, 524.05], tol=0.01)
#numpartsArray4e = findCumulativeNumParts(['/data/rob1g10/SPH/DATA/DISKS/40d_0.4e/disk005',
#                                       '/data/rob1g10/SPH/DATA/DISKS/40d_0.4e/disk006'],
#                                       None, 24 * 3600 / 5.5004e5, timeStart=[296.25, 610.25], tol=0.01)
#plt.plot(range(3800), numpartsArray0e[0:3800]/max(numpartsArray0e), label='40d 0.0e')
#plt.plot(range(3800), numpartsArray2e[0:3800]/max(numpartsArray2e), label='40d 0.2e')
#plt.plot(range(3800), numpartsArray4e[0:3800]/max(numpartsArray4e), label='40d 0.4e')
#plt.legend()
#plt.title('Fraction of particles remaining after injection is shut off')
#plt.xlabel('Days')
#plt.ylabel('Number of particles divided by number of particles at shut off')
#plt.show()

"""-----------------------fromStart-----------------------------------------"""

#numpartsArray0e = findCumulativeNumParts(['/data/rob1g10/SPH/DATA/DISKS/40d_0.0e/disk001',
#                                       '/data/rob1g10/SPH/DATA/DISKS/40d_0.0e/disk002'],
#                                       None, 24 * 3600 / 5.5004e5, timeStart=[0.05, 120.65], tol=0.01)
#numpartsArray2e = findCumulativeNumParts(['/data/rob1g10/SPH/DATA/DISKS/40d_0.2e/disk001',
#                                       '/data/rob1g10/SPH/DATA/DISKS/40d_0.2e/disk002'],
#                                       None, 24 * 3600 / 5.5004e5, timeStart=[0.05, 118.05], tol=0.01)
#numpartsArray4e = findCumulativeNumParts(['/data/rob1g10/SPH/DATA/DISKS/40d_0.4e/disk001',
#                                       '/data/rob1g10/SPH/DATA/DISKS/40d_0.4e/disk002'],
#                                       None, 24 * 3600 / 5.5004e5, timeStart=[0.05, 120.45], tol=0.01)
#plt.plot(range(400), numpartsArray0e[0:400], label='40d 0.0e')
#plt.plot(range(400), numpartsArray2e[0:400], label='40d 0.2e')
#plt.plot(range(400), numpartsArray4e[0:400], label='40d 0.4e')
#plt.legend()
#plt.title('Fraction of particles remaining after injection is shut off')
#plt.xlabel('Days')
#plt.ylabel('Number of particles divided by number of particles at shut off')
#plt.show()

"""----------------------------20d,40d,60d----------------------------------"""
"""------------------------------NSaccretion--------------------------------"""
"""-----------------------injection_off-------------------------------------"""


#NSarray20d_0e = createAccretionArray(['/data/rob1g10/SPH/DATA/DISKS/20d_0.0e/disk005',
#                                        '/data/rob1g10/SPH/DATA/DISKS/20d_0.0e/disk006'],
#                             "NS", None, 'day', tol=0.01)
#NSarray40d_0e = createAccretionArray(['/data/rob1g10/SPH/DATA/DISKS/40d_0.0e/disk005',
#                                        '/data/rob1g10/SPH/DATA/DISKS/40d_0.0e/disk006'],
#                             "NS", None, 'day', tol=0.01)
#NSarray60d_0e = createAccretionArray(['/data/rob1g10/SPH/DATA/DISKS/60d_0.0e/disk006',
#                                        '/data/rob1g10/SPH/DATA/DISKS/60d_0.0e/disk007'],
#                             "NS", None, 'day', tol=0.01)
#NSarray20d_2e = createAccretionArray(['/data/rob1g10/SPH/DATA/DISKS/20d_0.2e/disk005',
#                                        '/data/rob1g10/SPH/DATA/DISKS/20d_0.2e/disk006'],
#                             "NS", None, 'day', tol=0.01)
#NSarray40d_2e = createAccretionArray(['/data/rob1g10/SPH/DATA/DISKS/40d_0.2e/disk005',
#                                        '/data/rob1g10/SPH/DATA/DISKS/40d_0.2e/disk006'],
#                             "NS", None, 'day', tol=0.01)
#NSarray60d_2e = createAccretionArray(['/data/rob1g10/SPH/DATA/DISKS/60d_0.2e/disk006',
#                                        '/data/rob1g10/SPH/DATA/DISKS/60d_0.2e/disk007'],
#                             "NS", None, 'day', tol=0.01)
#NSarray20d_4e = createAccretionArray(['/data/rob1g10/SPH/DATA/DISKS/20d_0.4e/disk005'],
#                             "NS", None, 'day', tol=0.01)
#NSarray40d_4e = createAccretionArray(['/data/rob1g10/SPH/DATA/DISKS/40d_0.4e/disk005',
#                                        '/data/rob1g10/SPH/DATA/DISKS/40d_0.4e/disk006'],
#                             "NS", None, 'day', tol=0.01)
#NSarray60d_4e = createAccretionArray(['/data/rob1g10/SPH/DATA/DISKS/60d_0.4e/disk006',
#                                        '/data/rob1g10/SPH/DATA/DISKS/60d_0.4e/disk007'],
#                             "NS", None, 'day', tol=0.01)

#NUM = 1000
#plt.plot(range(NUM), NSarray20d_0e[0:NUM], label='20d 0.0e')
#plt.plot(range(NUM), NSarray40d_0e[0:NUM], label='40d 0.0e')
#plt.plot(range(NUM), NSarray60d_0e[0:NUM], label='60d 0.0e')
#plt.plot(range(NUM), NSarray20d_2e[0:NUM], label='20d 0.2e')
#plt.plot(range(NUM), NSarray40d_2e[0:NUM], label='40d 0.2e')
#plt.plot(range(NUM), NSarray60d_2e[0:NUM], label='60d 0.2e')
#plt.plot(range(NUM), NSarray20d_4e[0:NUM], label='20d 0.4e')
#plt.plot(range(NUM), NSarray40d_4e[0:NUM], label='40d 0.4e')
#plt.plot(range(NUM), NSarray60d_4e[0:NUM], label='60d 0.4e')
#plt.legend()
#plt.title('NS accretion over ~{} days after injection is shut off'.format(NUM))
#plt.xlabel('Days')
#plt.ylabel('Number of particles accreted by neutron star')
#plt.show()

"""------------------------------Beaccretion--------------------------------"""
"""-----------------------injection_off-------------------------------------"""


#BEarray20d_0e = createAccretionArray(['/data/rob1g10/SPH/DATA/DISKS/20d_0.0e/disk005',
#                                        '/data/rob1g10/SPH/DATA/DISKS/20d_0.0e/disk006'],
#                             "Be", None, 'day', tol=0.01)
#BEarray40d_0e = createAccretionArray(['/data/rob1g10/SPH/DATA/DISKS/40d_0.0e/disk005',
#                                        '/data/rob1g10/SPH/DATA/DISKS/40d_0.0e/disk006'],
#                             "Be", None, 'day', tol=0.01)
#BEarray60d_0e = createAccretionArray(['/data/rob1g10/SPH/DATA/DISKS/60d_0.0e/disk006',
#                                        '/data/rob1g10/SPH/DATA/DISKS/60d_0.0e/disk007'],
#                             "Be", None, 'day', tol=0.01)
#BEarray20d_2e = createAccretionArray(['/data/rob1g10/SPH/DATA/DISKS/20d_0.2e/disk005',
#                                        '/data/rob1g10/SPH/DATA/DISKS/20d_0.2e/disk006'],
#                             "Be", None, 'day', tol=0.01)
#BEarray40d_2e = createAccretionArray(['/data/rob1g10/SPH/DATA/DISKS/40d_0.2e/disk005',
#                                        '/data/rob1g10/SPH/DATA/DISKS/40d_0.2e/disk006'],
#                             "Be", None, 'day', tol=0.01)
#BEarray60d_2e = createAccretionArray(['/data/rob1g10/SPH/DATA/DISKS/60d_0.2e/disk006',
#                                        '/data/rob1g10/SPH/DATA/DISKS/60d_0.2e/disk007'],
#                             "Be", None, 'day', tol=0.01)
#BEarray20d_4e = createAccretionArray(['/data/rob1g10/SPH/DATA/DISKS/20d_0.4e/disk005'],
#                             "Be", None, 'day', tol=0.01)
#BEarray40d_4e = createAccretionArray(['/data/rob1g10/SPH/DATA/DISKS/40d_0.4e/disk005',
#                                        '/data/rob1g10/SPH/DATA/DISKS/40d_0.4e/disk006'],
#                             "Be", None, 'day', tol=0.01)
#BEarray60d_4e = createAccretionArray(['/data/rob1g10/SPH/DATA/DISKS/60d_0.4e/disk006',
#                                        '/data/rob1g10/SPH/DATA/DISKS/60d_0.4e/disk007'],
#                             "Be", None, 'day', tol=0.01)

#NUM = 1000
#plt.plot(range(NUM), BEarray20d_0e[0:NUM], label='20d 0.0e')
#plt.plot(range(NUM), BEarray40d_0e[0:NUM], label='40d 0.0e')
#plt.plot(range(NUM), BEarray60d_0e[0:NUM], label='60d 0.0e')
#plt.plot(range(NUM), BEarray20d_2e[0:NUM], label='20d 0.2e')
#plt.plot(range(NUM), BEarray40d_2e[0:NUM], label='40d 0.2e')
#plt.plot(range(NUM), BEarray60d_2e[0:NUM], label='60d 0.2e')
#plt.plot(range(NUM), BEarray20d_4e[0:NUM], label='20d 0.4e')
#plt.plot(range(NUM), BEarray40d_4e[0:NUM], label='40d 0.4e')
#plt.plot(range(NUM), BEarray60d_4e[0:NUM], label='60d 0.4e')
#plt.legend()
#plt.title('Be accretion over ~{} days after injection is shut off'.format(NUM))
#plt.xlabel('Days')
#plt.ylabel('Number of particles accreted by Be star')

"""-----------------------fromStart-----------------------------------------"""



"""-----------------------numActiveParts------------------------------------"""
"""-----------------------injection_off-------------------------------------"""


#partarray20d_0e = findCumulativeNumParts(['/data/rob1g10/SPH/DATA/DISKS/20d_0.0e/disk005',
#                                        '/data/rob1g10/SPH/DATA/DISKS/20d_0.0e/disk006'],
#                                        None, 'day', tol=0.01)
#partarray40d_0e = findCumulativeNumParts(['/data/rob1g10/SPH/DATA/DISKS/40d_0.0e/disk005',
#                                        '/data/rob1g10/SPH/DATA/DISKS/40d_0.0e/disk006'],
#                                        None, 'day', tol=0.01)
#partarray60d_0e = findCumulativeNumParts(['/data/rob1g10/SPH/DATA/DISKS/60d_0.0e/disk006',
#                                        '/data/rob1g10/SPH/DATA/DISKS/60d_0.0e/disk007'],
#                                        None, 'day', tol=0.01)
#partarray20d_2e = findCumulativeNumParts(['/data/rob1g10/SPH/DATA/DISKS/20d_0.2e/disk005',
#                                        '/data/rob1g10/SPH/DATA/DISKS/20d_0.2e/disk006'],
#                                        None, 'day', tol=0.01)
#partarray40d_2e = findCumulativeNumParts(['/data/rob1g10/SPH/DATA/DISKS/40d_0.2e/disk005',
#                                        '/data/rob1g10/SPH/DATA/DISKS/40d_0.2e/disk006'],
#                                        None, 'day', tol=0.01)
#partarray60d_2e = findCumulativeNumParts(['/data/rob1g10/SPH/DATA/DISKS/60d_0.2e/disk006',
#                                        '/data/rob1g10/SPH/DATA/DISKS/60d_0.2e/disk007'],
#                                        None, 'day', tol=0.01)
#partarray20d_4e = findCumulativeNumParts(['/data/rob1g10/SPH/DATA/DISKS/20d_0.4e/disk005'],
#                                       None, 'day', tol=0.01)
#partarray40d_4e = findCumulativeNumParts(['/data/rob1g10/SPH/DATA/DISKS/40d_0.4e/disk005',
#                                        '/data/rob1g10/SPH/DATA/DISKS/40d_0.4e/disk006'],
#                                        None, 'day', tol=0.01)
#partarray60d_4e = findCumulativeNumParts(['/data/rob1g10/SPH/DATA/DISKS/60d_0.4e/disk006',
#                                        '/data/rob1g10/SPH/DATA/DISKS/60d_0.4e/disk007'],
#                                        None, 'day', tol=0.01)
#NUM = 1400
#plt.plot(range(NUM), partarray20d_0e[0:NUM]/max(partarray20d_0e), r'-', label='20d 0.0e')
#plt.plot(range(NUM), partarray40d_0e[0:NUM]/max(partarray40d_0e), r'-', label='40d 0.0e')
#plt.plot(range(NUM), partarray60d_0e[0:NUM]/max(partarray60d_0e), r'-', label='60d 0.0e')
#plt.plot(range(NUM), partarray20d_2e[0:NUM]/max(partarray20d_2e), r'x', label='20d 0.2e')
#plt.plot(range(NUM), partarray40d_2e[0:NUM]/max(partarray40d_2e), r'x', label='40d 0.2e')
#plt.plot(range(NUM), partarray60d_2e[0:NUM]/max(partarray60d_2e), r'x', label='60d 0.2e')
#plt.plot(range(NUM), partarray20d_4e[0:NUM]/max(partarray20d_4e), r'o', label='20d 0.4e')
#plt.plot(range(NUM), partarray40d_4e[0:NUM]/max(partarray40d_4e), r'o', label='40d 0.4e')
#plt.plot(range(NUM), partarray60d_4e[0:NUM]/max(partarray60d_4e), r'o', label='60d 0.4e')
#plt.legend()
#plt.title('Fraction of active particles over ~300 days after injection is shut off')
#plt.xlabel('Days')
#plt.ylabel('Number of particles / Number of active particles at injection off')
#plt.show()


"""-----------------------fromStart-----------------------------------------"""

#numpartsArray0e = findCumulativeNumParts(['/data/rob1g10/SPH/DATA/DISKS/40d_0.0e/disk001',
#                                       '/data/rob1g10/SPH/DATA/DISKS/40d_0.0e/disk002'],
#                                       None, 24 * 3600 / 5.5004e5, timeStart=[0.05, 120.65], tol=0.01)
#numpartsArray2e = findCumulativeNumParts(['/data/rob1g10/SPH/DATA/DISKS/40d_0.2e/disk001',
#                                       '/data/rob1g10/SPH/DATA/DISKS/40d_0.2e/disk002'],
#                                       None, 24 * 3600 / 5.5004e5, timeStart=[0.05, 118.05], tol=0.01)
#numpartsArray4e = findCumulativeNumParts(['/data/rob1g10/SPH/DATA/DISKS/40d_0.4e/disk001',
#                                       '/data/rob1g10/SPH/DATA/DISKS/40d_0.4e/disk002'],
#                                       None, 24 * 3600 / 5.5004e5, timeStart=[0.05, 120.45], tol=0.01)
#plt.plot(range(400), numpartsArray0e[0:400], label='40d 0.0e')
#plt.plot(range(400), numpartsArray2e[0:400], label='40d 0.2e')
#plt.plot(range(400), numpartsArray4e[0:400], label='40d 0.4e')
#plt.legend()
#plt.title('Fraction of particles remaining after injection is shut off')
#plt.xlabel('Days')
#plt.ylabel('Number of particles divided by number of particles at shut off')
#plt.show()

"""-----------------------fromStart 80d-----------------------------------------"""

#numpartsArray0e = findCumulativeNumParts(['/data/rob1g10/SPH/DATA/DISKS/80d_0.0e/disk001',
#                                       '/data/rob1g10/SPH/DATA/DISKS/80d_0.0e/disk002',
#                                       '/data/rob1g10/SPH/DATA/DISKS/80d_0.0e/disk003',
#                                       '/data/rob1g10/SPH/DATA/DISKS/80d_0.0e/disk004',
#                                       '/data/rob1g10/SPH/DATA/DISKS/80d_0.0e/disk005',
#                                       '/data/rob1g10/SPH/DATA/DISKS/80d_0.0e/disk006',
#                                       '/data/rob1g10/SPH/DATA/DISKS/80d_0.0e/disk007',
#                                       '/data/rob1g10/SPH/DATA/DISKS/80d_0.0e/disk008',
#                                       '/data/rob1g10/SPH/DATA/DISKS/80d_0.0e/disk009',
#                                       '/data/rob1g10/SPH/DATA/DISKS/80d_0.0e/disk010'],
#                                       None, 'day', tol=0.01)
#numpartsArray2e = findCumulativeNumParts(['/data/rob1g10/SPH/DATA/DISKS/80d_0.2e/disk001',
#                                       '/data/rob1g10/SPH/DATA/DISKS/80d_0.2e/disk002',
#                                       '/data/rob1g10/SPH/DATA/DISKS/80d_0.2e/disk003',
#                                       '/data/rob1g10/SPH/DATA/DISKS/80d_0.2e/disk004',
#                                       '/data/rob1g10/SPH/DATA/DISKS/80d_0.2e/disk005',
#                                       '/data/rob1g10/SPH/DATA/DISKS/80d_0.2e/disk006',
#                                       '/data/rob1g10/SPH/DATA/DISKS/80d_0.2e/disk007',
#                                       '/data/rob1g10/SPH/DATA/DISKS/80d_0.2e/disk008',
#                                       '/data/rob1g10/SPH/DATA/DISKS/80d_0.2e/disk009',
#                                       '/data/rob1g10/SPH/DATA/DISKS/80d_0.2e/disk010'],
#                                       None, 'day', tol=0.01)
#numpartsArray4e = findCumulativeNumParts(['/data/rob1g10/SPH/DATA/DISKS/80d_0.4e/disk001',
#                                       '/data/rob1g10/SPH/DATA/DISKS/80d_0.4e/disk002',
#                                       '/data/rob1g10/SPH/DATA/DISKS/80d_0.4e/disk003',
#                                       '/data/rob1g10/SPH/DATA/DISKS/80d_0.4e/disk004',
#                                       '/data/rob1g10/SPH/DATA/DISKS/80d_0.4e/disk005',
#                                       '/data/rob1g10/SPH/DATA/DISKS/80d_0.4e/disk006',
#                                       '/data/rob1g10/SPH/DATA/DISKS/80d_0.4e/disk007',
#                                       '/data/rob1g10/SPH/DATA/DISKS/80d_0.4e/disk008',
#                                       '/data/rob1g10/SPH/DATA/DISKS/80d_0.4e/disk009',
#                                       '/data/rob1g10/SPH/DATA/DISKS/80d_0.4e/disk010'],
#                                       None, 'day', tol=0.01)

#numpartsArray4e = findCumulativeNumParts(['/data/rob1g10/SPH/DATA/DISKS/80d_0.4e_injFixed/disk001',
#                                          '/data/rob1g10/SPH/DATA/DISKS/80d_0.4e_injFixed/disk002',
#                                          '/data/rob1g10/SPH/DATA/DISKS/80d_0.4e_injFixed/disk003',
#                                          '/data/rob1g10/SPH/DATA/DISKS/80d_0.4e_injFixed/disk004',
#                                          '/data/rob1g10/SPH/DATA/DISKS/80d_0.4e_injFixed/disk005'],
#                                       None, 'day', tol=0.01)
#plt.plot(range(len(numpartsArray4e)), numpartsArray4e, label='Num parts')
                                       
#print len(numpartsArray0e)
#plt.plot(range(len(numpartsArray0e)), (numpartsArray0e - 1e5) / 10, label='80d 0.0e - 1e5 / 10')
#plt.plot(range(len(numpartsArray2e)), (numpartsArray2e - 1e5) / 10, label='80d 0.2e - 1e5 / 10')
#plt.plot(range(len(numpartsArray4e)), (numpartsArray4e - 1e5) / 10, label='80d 0.4e - 1e5 / 10')
#plt.legend()
#plt.title('Number of particles')
#plt.xlabel('Days')
#plt.ylabel('Number of particles')
#plt.show()

#array0e = createAccretionArray(['/data/rob1g10/SPH/DATA/DISKS/80d_0.0e/disk001',
#                                       '/data/rob1g10/SPH/DATA/DISKS/80d_0.0e/disk002',
#                                       '/data/rob1g10/SPH/DATA/DISKS/80d_0.0e/disk003',
#                                       '/data/rob1g10/SPH/DATA/DISKS/80d_0.0e/disk004',
#                                       '/data/rob1g10/SPH/DATA/DISKS/80d_0.0e/disk005',
#                                       '/data/rob1g10/SPH/DATA/DISKS/80d_0.0e/disk006',
#                                       '/data/rob1g10/SPH/DATA/DISKS/80d_0.0e/disk007',
#                                       '/data/rob1g10/SPH/DATA/DISKS/80d_0.0e/disk008',
#                                       '/data/rob1g10/SPH/DATA/DISKS/80d_0.0e/disk009',
#                                       '/data/rob1g10/SPH/DATA/DISKS/80d_0.0e/disk010'],
#                                       "Be", None, 'day', tol=0.01)
#array2e = createAccretionArray(['/data/rob1g10/SPH/DATA/DISKS/80d_0.2e/disk001',
#                                       '/data/rob1g10/SPH/DATA/DISKS/80d_0.2e/disk002',
#                                       '/data/rob1g10/SPH/DATA/DISKS/80d_0.2e/disk003',
#                                       '/data/rob1g10/SPH/DATA/DISKS/80d_0.2e/disk004',
#                                       '/data/rob1g10/SPH/DATA/DISKS/80d_0.2e/disk005',
#                                       '/data/rob1g10/SPH/DATA/DISKS/80d_0.2e/disk006',
#                                       '/data/rob1g10/SPH/DATA/DISKS/80d_0.2e/disk007',
#                                       '/data/rob1g10/SPH/DATA/DISKS/80d_0.2e/disk008',
#                                       '/data/rob1g10/SPH/DATA/DISKS/80d_0.2e/disk009',
#                                       '/data/rob1g10/SPH/DATA/DISKS/80d_0.2e/disk010'],
#                                       "Be", None, 'day', tol=0.01)
#array4e = createAccretionArray(['/data/rob1g10/SPH/DATA/DISKS/80d_0.4e/disk001',
#                                       '/data/rob1g10/SPH/DATA/DISKS/80d_0.4e/disk002',
#                                       '/data/rob1g10/SPH/DATA/DISKS/80d_0.4e/disk003',
#                                       '/data/rob1g10/SPH/DATA/DISKS/80d_0.4e/disk004',
#                                       '/data/rob1g10/SPH/DATA/DISKS/80d_0.4e/disk005',
#                                       '/data/rob1g10/SPH/DATA/DISKS/80d_0.4e/disk006',
#                                       '/data/rob1g10/SPH/DATA/DISKS/80d_0.4e/disk007',
#                                       '/data/rob1g10/SPH/DATA/DISKS/80d_0.4e/disk008',
#                                       '/data/rob1g10/SPH/DATA/DISKS/80d_0.4e/disk009',
#                                       '/data/rob1g10/SPH/DATA/DISKS/80d_0.4e/disk010'],
#                                       "Be", None, 'day', tol=0.01)

#array4e = createAccretionArray(['/data/rob1g10/SPH/DATA/DISKS/80d_0.4e_injFixed/disk001',
#                                '/data/rob1g10/SPH/DATA/DISKS/80d_0.4e_injFixed/disk002',
#                                '/data/rob1g10/SPH/DATA/DISKS/80d_0.4e_injFixed/disk003',
#                                '/data/rob1g10/SPH/DATA/DISKS/80d_0.4e_injFixed/disk004',
#                                '/data/rob1g10/SPH/DATA/DISKS/80d_0.4e_injFixed/disk005'],
#                               "Be", None, 'day', tol=0.01)

#plt.plot(range(len(array4e)), array4e, label='Be accretion')
#plt.plot(range(len(array0e)), array0e, label='80d 0.0e accretion')
#plt.plot(range(len(array2e)), array2e, label='80d 0.2e accretion')
#plt.plot(range(len(array4e)), array4e, label='80d 0.4e accretion')
#plt.legend()
#plt.title('Be accretion')
#plt.xlabel('Days')
#plt.ylabel('Number of particles accreted by Be star')

"""-----------------------fromStart isoBe-----------------------------------------"""

#numpartsArray = findCumulativeNumParts(['/data/rob1g10/SPH/DATA/DISKS/isoBe/disk001',
#                                       '/data/rob1g10/SPH/DATA/DISKS/isoBe/disk002',
#                                       '/data/rob1g10/SPH/DATA/DISKS/isoBe/disk003',
#                                       '/data/rob1g10/SPH/DATA/DISKS/isoBe/disk004',
#                                       '/data/rob1g10/SPH/DATA/DISKS/isoBe/disk005',
#                                       '/data/rob1g10/SPH/DATA/DISKS/isoBe/disk006',
#                                       '/data/rob1g10/SPH/DATA/DISKS/isoBe/disk007',
#                                       '/data/rob1g10/SPH/DATA/DISKS/isoBe/disk008',
#                                        '/data/rob1g10/SPH/DATA/DISKS/isoBe/disk009'],
#                                       None, 'day', timeStart=[0.05, 120.65], tol=0.01)

#numpartsArray = findCumulativeNumParts(['/data/rob1g10/SPH/DATA/DISKS/isoBe_fromBeginningTenth/disk001',
#                                        '/data/rob1g10/SPH/DATA/DISKS/isoBe_fromBeginningTenth/disk002'],
#                                       None, 'day', tol=0.01)
#numpartsArray = findCumulativeNumParts(['/data/rob1g10/SPH/DATA/DISKS/20d_0.0e_eqTest/disk027',
#                                        '/data/rob1g10/SPH/DATA/DISKS/20d_0.0e_eqTest/disk028',
#                                        '/data/rob1g10/SPH/DATA/DISKS/20d_0.0e_eqTest/disk029'],
#                                       None, 'day', tol=0.01)
#numpartsArray = findCumulativeNumParts(['/data/rob1g10/SPH/DATA/DISKS/40d_0.2e_eqTest/disk050',
#                                        '/data/rob1g10/SPH/DATA/DISKS/40d_0.2e_eqTest/disk051',
#                                        '/data/rob1g10/SPH/DATA/DISKS/40d_0.2e_eqTest/disk052'],
#                                       None, 'day', tol=0.01)
#numpartsArray = findCumulativeNumParts(['/data/rob1g10/SPH/DATA/DISKS/60d_0.4e_eqTest/disk011',
#                                        '/data/rob1g10/SPH/DATA/DISKS/60d_0.4e_eqTest/disk012',
#                                        '/data/rob1g10/SPH/DATA/DISKS/60d_0.4e_eqTest/disk013'],
#                                       None, 'day', tol=0.01)
#numpartsArray = findCumulativeNumParts(['/data/rob1g10/SPH/DATA/DISKS/80d_0.4e_eqTest/disk004',
#                                        '/data/rob1g10/SPH/DATA/DISKS/80d_0.4e_eqTest/disk005',
#                                        '/data/rob1g10/SPH/DATA/DISKS/80d_0.4e_eqTest/disk006'],
#                                       None, 'day', tol=0.01)
#numpartsArray = findCumulativeNumParts(['disk005'],
#                                       None, 'day', tol=0.01)

#plt.plot(range(len(numpartsArray)), (numpartsArray), label='Number of particles')
#plt.title('Number of particles')
#plt.xlabel('Days')
#plt.ylabel('Number of particles')
#plt.ylim(0, 2e5)
#plt.show()

#plt.plot( [0, 30, 60, 90], [5.4, 25.01, 28.15, 28.67])
#p = np.polyfit([0, 30, 60, 90], [5.4, 25.01, 28.15, 28.67],2)
#r = np.arange(0,90)
#plt.plot(range(0,90), 5.4 + 23.27*np.sin(np.radians(range(0,90))))
#plt.plot(r, p[0]*r*r + p[1]*r + p[0])
#plt.xlabel('Inclination Angle')
#plt.ylabel('Mean negative EW')
#plt.ylim(0, 2e5)
#plt.show()

#BEarray = createAccretionArray(['/data/rob1g10/SPH/DATA/DISKS/isoBe/disk001',
#                               '/data/rob1g10/SPH/DATA/DISKS/isoBe/disk002',
#                               '/data/rob1g10/SPH/DATA/DISKS/isoBe/disk003',
#                               '/data/rob1g10/SPH/DATA/DISKS/isoBe/disk004',
#                               '/data/rob1g10/SPH/DATA/DISKS/isoBe/disk005',
#                               '/data/rob1g10/SPH/DATA/DISKS/isoBe/disk006',
#                               '/data/rob1g10/SPH/DATA/DISKS/isoBe/disk007',
#                               '/data/rob1g10/SPH/DATA/DISKS/isoBe/disk008',
#                                '/data/rob1g10/SPH/DATA/DISKS/isoBe/disk009'],
#                             "Be", None, 'day', tol=0.01)
#
#BEarray2 = createAccretionArray(['/data/rob1g10/SPH/DATA/DISKS/isoBe_injFixed/disk001','/data/rob1g10/SPH/DATA/DISKS/isoBe_injFixed/disk002'],
#                             "Be", None, 'day', tol=0.01)
#
#plt.plot(range(len(BEarray)), BEarray, label='Be accretion')
#plt.plot(range(len(BEarray2)), BEarray2, label='Be accretion')
#plt.legend()
#plt.title('Be accretion')
#plt.xlabel('Days')
#plt.ylabel('Number of particles accreted by Be star')

"""-----------------------Scratch-------------------------------------------"""


#BEarrayhalf = createAccretionArray([['/data/rob1g10/SPH/DATA/DISKS/justBe/disk007',
#                                    '/data/rob1g10/SPH/DATA/DISKS/justBe/disk008'][1]],
#                             "Be", None, 'day', tol=0.01)
#BEarrayquarter = createAccretionArray([['/data/rob1g10/SPH/DATA/DISKS/justBe_quarterInj/disk005',
#                                       '/data/rob1g10/SPH/DATA/DISKS/justBe_quarterInj/disk006'][1]],
#                             "Be", None, 'day', tol=0.01)
#BEarrayeighth = createAccretionArray([['/data/rob1g10/SPH/DATA/DISKS/justBe_eighthInj/disk004',
#                                      '/data/rob1g10/SPH/DATA/DISKS/justBe_eighthInj2/disk005'][1]],
#                             "Be", None, 'day', tol=0.01)
##numpartsArrayhalf = findCumulativeNumParts([['/data/rob1g10/SPH/DATA/DISKS/justBe_halfInj2/disk011',
##                                    '/data/rob1g10/SPH/DATA/DISKS/justBe_halfInj2/disk012'][1]],
##                                       None, 'day', tol=0.01)
##numpartsArrayquarter = findCumulativeNumParts([['/data/rob1g10/SPH/DATA/DISKS/justBe_quarterInj/disk005',
##                                       '/data/rob1g10/SPH/DATA/DISKS/justBe_quarterInj/disk006'][1]],
##                                       None, 'day', timeStart=[0.05, 118.05], tol=0.01)
##numpartsArrayeighth = findCumulativeNumParts([['/data/rob1g10/SPH/DATA/DISKS/justBe_eighthInj/disk004',
##                                      '/data/rob1g10/SPH/DATA/DISKS/justBe_eighthInj2/disk005'][1]],
##                                       None, 'day', tol=0.01)
#lens = [len(BEarrayhalf), len(BEarrayquarter), len(BEarrayeighth)]
##lens = [len(numpartsArrayhalf), len(numpartsArrayquarter), len(numpartsArrayeighth)]
#NUM = min(lens)
#plt.plot(range(NUM), BEarrayhalf[0:NUM], label='standard')
#plt.plot(range(NUM), BEarrayquarter[0:NUM], label='half')
#plt.plot(range(NUM), BEarrayeighth[0:NUM], label='quarter')
##plt.plot(range(NUM), numpartsArrayhalf[0:NUM], label='standard')
##plt.plot(range(NUM), numpartsArrayquarter[0:NUM], label='half')
##plt.plot(range(NUM), numpartsArrayeighth[0:NUM], label='quarter')
#plt.legend()
#plt.title('Be accretion over ~{} days'.format(NUM))
#plt.xlabel('Days')
#plt.ylabel('Number of particles accreted by Be star')


#array = createAccretionArray(['/data/rob1g10/SPH/DATA/DISKS/20d_0.2e_eqTest/disk001',
#                              '/data/rob1g10/SPH/DATA/DISKS/20d_0.2e_eqTest/disk002',
#                              '/data/rob1g10/SPH/DATA/DISKS/20d_0.2e_eqTest/disk003',
#                              '/data/rob1g10/SPH/DATA/DISKS/20d_0.2e_eqTest/disk004',
#                              '/data/rob1g10/SPH/DATA/DISKS/20d_0.2e_eqTest/disk005',
#                              '/data/rob1g10/SPH/DATA/DISKS/20d_0.2e_eqTest/disk006',
#                              '/data/rob1g10/SPH/DATA/DISKS/20d_0.2e_eqTest/disk007',
#                              '/data/rob1g10/SPH/DATA/DISKS/20d_0.2e_eqTest/disk008',
#                              '/data/rob1g10/SPH/DATA/DISKS/20d_0.2e_eqTest/disk009',
#                              '/data/rob1g10/SPH/DATA/DISKS/20d_0.2e_eqTest/disk010',
#                              '/data/rob1g10/SPH/DATA/DISKS/20d_0.2e_eqTest/disk011',
#                              '/data/rob1g10/SPH/DATA/DISKS/20d_0.2e_eqTest/disk012'],
#                                       "NS", None, 'day', tol=0.01)
#BH0292array = createAccretionArray(['/data/rob1g10/SPH/DATA/DISKS/BH20d_0.292e/disk001',
#                                '/data/rob1g10/SPH/DATA/DISKS/BH20d_0.292e/disk002'],
#                                       "Be", None, 'day', tol=0.01)
#BH167array = createAccretionArray(['/data/rob1g10/SPH/DATA/DISKS/BH16.7d_0.2e/disk001',
#                                '/data/rob1g10/SPH/DATA/DISKS/BH16.7d_0.2e/disk002'],
#                                       "Be", None, 'day', tol=0.01)
#BH143array = createAccretionArray(['/data/rob1g10/SPH/DATA/DISKS/BH14.3d_0.2e/disk001',
#                                '/data/rob1g10/SPH/DATA/DISKS/BH14.3d_0.2e/disk002'],
#                                       "Be", None, 'day', tol=0.01)
#BH28array = createAccretionArray(['/data/rob1g10/SPH/DATA/DISKS/BH28d_0.2e/disk001',
#                                '/data/rob1g10/SPH/DATA/DISKS/BH28d_0.2e/disk002',
#                                '/data/rob1g10/SPH/DATA/DISKS/BH28d_0.2e/disk003'],
#                                       "NS", None, 'day', tol=0.01)
#BH24array = createAccretionArray(['/data/rob1g10/SPH/DATA/DISKS/BH24d_0.2e/disk001',
#                                '/data/rob1g10/SPH/DATA/DISKS/BH24d_0.2e/disk002',
#                                '/data/rob1g10/SPH/DATA/DISKS/BH24d_0.2e/disk003'],
#                                       "NS", None, 'day', tol=0.01)
#BH24_3array = createAccretionArray(['/data/rob1g10/SPH/DATA/DISKS/BH24.3d_0.2e/disk001'],
#                                       "NS", None, 'day', tol=0.01)
#BH25array = createAccretionArray(['/data/rob1g10/SPH/DATA/DISKS/BH25d_0.2e/disk001'],
#                                       "NS", None, 'day', tol=0.01)
#array = findCumulativeNumParts(['/data/rob1g10/SPH/DATA/DISKS/20d_0.2e_eqTest/disk001',
#                                '/data/rob1g10/SPH/DATA/DISKS/20d_0.2e_eqTest/disk002',
#                                '/data/rob1g10/SPH/DATA/DISKS/20d_0.2e_eqTest/disk003',
#                                '/data/rob1g10/SPH/DATA/DISKS/20d_0.2e_eqTest/disk004',
#                                '/data/rob1g10/SPH/DATA/DISKS/20d_0.2e_eqTest/disk005',
#                                '/data/rob1g10/SPH/DATA/DISKS/20d_0.2e_eqTest/disk006',
#                                '/data/rob1g10/SPH/DATA/DISKS/20d_0.2e_eqTest/disk007',
#                                '/data/rob1g10/SPH/DATA/DISKS/20d_0.2e_eqTest/disk008',
#                                '/data/rob1g10/SPH/DATA/DISKS/20d_0.2e_eqTest/disk009',
#                                '/data/rob1g10/SPH/DATA/DISKS/20d_0.2e_eqTest/disk010',
#                                '/data/rob1g10/SPH/DATA/DISKS/20d_0.2e_eqTest/disk011',
#                                '/data/rob1g10/SPH/DATA/DISKS/20d_0.2e_eqTest/disk012'],
#                                       None, 'day', tol=0.01)                                       
#BH0292array = findCumulativeNumParts(['/data/rob1g10/SPH/DATA/DISKS/BH20d_0.292e/disk001',
#                                  '/data/rob1g10/SPH/DATA/DISKS/BH20d_0.292e/disk002'],
#                                       None, 'day', tol=0.01)
#BH167array = findCumulativeNumParts(['/data/rob1g10/SPH/DATA/DISKS/BH16.7d_0.2e/disk001',
#                                '/data/rob1g10/SPH/DATA/DISKS/BH16.7d_0.2e/disk002'],
#                                       None, 'day', tol=0.01)
#BH143array = findCumulativeNumParts(['/data/rob1g10/SPH/DATA/DISKS/BH14.3d_0.2e/disk001',
#                                '/data/rob1g10/SPH/DATA/DISKS/BH14.3d_0.2e/disk002'],
#                                       None, 'day', tol=0.01)
#BH28array = findCumulativeNumParts(['/data/rob1g10/SPH/DATA/DISKS/BH28d_0.2e/disk001',
#                                '/data/rob1g10/SPH/DATA/DISKS/BH28d_0.2e/disk002',
#                                '/data/rob1g10/SPH/DATA/DISKS/BH24d_0.2e/disk003'],
#                                       None, 'day', tol=0.01)
#BH24array = findCumulativeNumParts(['/data/rob1g10/SPH/DATA/DISKS/BH24d_0.2e/disk001',
#                                '/data/rob1g10/SPH/DATA/DISKS/BH24d_0.2e/disk002',
#                                '/data/rob1g10/SPH/DATA/DISKS/BH24d_0.2e/disk003'],
#                                       None, 'day', tol=0.01)
#                                       
#plt.plot(range(len(array)), array, label='20d 0.2e')
#plt.plot(range(len(BH24array)), BH24array, label='BH 24d 0.2e')
#plt.plot(range(len(BH28array)), BH28array, label='BH 28d 0.2e')
#plt.plot(range(len(BH24_3array)), BH24_3array, label='BH 24.3d 0.2e')
#plt.plot(range(len(BH25array)), BH25array, label='BH 25d 0.2e')
#plt.plot(range(len(BH0292array)), BH0292array, label='BH 20d 0.292e')
#plt.plot(range(len(BH167array)), BH167array, label='BH 16.7d 0.2e')
#plt.plot(range(len(BH143array)), BH143array, label='BH 14.3d 0.2e')
#for i in range(10, 2600, 20):
#    plt.plot(i * np.ones(10), np.linspace(np.min(array), np.max(array), 10), color='orange')
#plt.legend(loc='best')
#plt.xlabel('Days')
#plt.ylabel('Number of particles accreted by Be star')  
#plt.xlim(3500, 3600)
#plt.ylim(0, 20000)       
#plt.ylim(150, 350)   
#plt.ylim(0, 15)
#print np.sum(array[3000:3400]), np.sum(BH24array[3000:3480]), np.sum(BH28array[3000:3560])

#------------------------------------------
# NS PLOUGH
#------------------------------------------

#array02 = findCumulativeNumParts(['/data/rob1g10/SPH/DATA/DISKS/20d_0.2e_eqTest/disk001',
#                                '/data/rob1g10/SPH/DATA/DISKS/20d_0.2e_eqTest/disk002',
#                                '/data/rob1g10/SPH/DATA/DISKS/20d_0.2e_eqTest/disk003',
#                                '/data/rob1g10/SPH/DATA/DISKS/20d_0.2e_eqTest/disk004',
#                                '/data/rob1g10/SPH/DATA/DISKS/20d_0.2e_eqTest/disk005',
#                                '/data/rob1g10/SPH/DATA/DISKS/20d_0.2e_eqTest/disk006',
#                                '/data/rob1g10/SPH/DATA/DISKS/20d_0.2e_eqTest/disk007',
#                                '/data/rob1g10/SPH/DATA/DISKS/20d_0.2e_eqTest/disk008',
#                                '/data/rob1g10/SPH/DATA/DISKS/20d_0.2e_eqTest/disk009',
#                                '/data/rob1g10/SPH/DATA/DISKS/20d_0.2e_eqTest/disk010',
#                                '/data/rob1g10/SPH/DATA/DISKS/20d_0.2e_eqTest/disk011',
#                                '/data/rob1g10/SPH/DATA/DISKS/20d_0.2e_eqTest/disk012'],
#                                       None, 'day', tol=0.01)
#array04 = findCumulativeNumParts(['/data/rob1g10/SPH/DATA/DISKS/20d_0.4e_eqTest/disk001',
#                                '/data/rob1g10/SPH/DATA/DISKS/20d_0.4e_eqTest/disk002',
#                                '/data/rob1g10/SPH/DATA/DISKS/20d_0.4e_eqTest/disk003'],
#                                       None, 'day', tol=0.01)
#array09 = findCumulativeNumParts(['/data/rob1g10/SPH/DATA/DISKS/20d_0.9e/disk001',
#                                '/data/rob1g10/SPH/DATA/DISKS/20d_0.9e/disk002'],
#                                       None, 'day', tol=0.01)
#array02 = createAccretionArray(['/data/rob1g10/SPH/DATA/DISKS/20d_0.2e_eqTest/disk001',
#                              '/data/rob1g10/SPH/DATA/DISKS/20d_0.2e_eqTest/disk002',
#                              '/data/rob1g10/SPH/DATA/DISKS/20d_0.2e_eqTest/disk003',
#                              '/data/rob1g10/SPH/DATA/DISKS/20d_0.2e_eqTest/disk004',
#                              '/data/rob1g10/SPH/DATA/DISKS/20d_0.2e_eqTest/disk005',
#                              '/data/rob1g10/SPH/DATA/DISKS/20d_0.2e_eqTest/disk006',
#                              '/data/rob1g10/SPH/DATA/DISKS/20d_0.2e_eqTest/disk007',
#                              '/data/rob1g10/SPH/DATA/DISKS/20d_0.2e_eqTest/disk008',
#                              '/data/rob1g10/SPH/DATA/DISKS/20d_0.2e_eqTest/disk009',
#                              '/data/rob1g10/SPH/DATA/DISKS/20d_0.2e_eqTest/disk010',
#                              '/data/rob1g10/SPH/DATA/DISKS/20d_0.2e_eqTest/disk011',
#                              '/data/rob1g10/SPH/DATA/DISKS/20d_0.2e_eqTest/disk012'],
#                                       "Be", None, 'day', tol=0.01)
#
#array04 = createAccretionArray(['/data/rob1g10/SPH/DATA/DISKS/20d_0.4e_eqTest/disk001',
#                              '/data/rob1g10/SPH/DATA/DISKS/20d_0.4e_eqTest/disk002',
#                              '/data/rob1g10/SPH/DATA/DISKS/20d_0.4e_eqTest/disk003'],
#                                       "Be", None, 'day', tol=0.01)
#array09 = createAccretionArray(['/data/rob1g10/SPH/DATA/DISKS/20d_0.9e/disk001',
#                              '/data/rob1g10/SPH/DATA/DISKS/20d_0.9e/disk002'],
#                                       "Be", None, 'day', tol=0.01)

#plt.plot(range(len(array02)), array02, label='20d 0.2e')                                     
#plt.plot(range(len(array04)), array04, label='20d 0.4e', color='g')
#plt.plot(range(len(array09)), array09, label='20d 0.9e', color='r')
#plt.xlabel('Days')
#plt.ylabel('Number of particles accreted by Be star') 
#plt.legend(loc='best')  
#plt.xlim(1400, 1600)
#plt.ylim(100, 750) 

#------------------------------------------
# BHs
#------------------------------------------


#array14 = findCumulativeNumParts(['/data/rob1g10/SPH/DATA/DISKS/M14_20d_0.0e/disk001',
#                                '/data/rob1g10/SPH/DATA/DISKS/M14_20d_0.0e/disk002',
#                                '/data/rob1g10/SPH/DATA/DISKS/M14_20d_0.0e/disk003'],
#                                       None, 'day', tol=0.01)
#array10 = findCumulativeNumParts(['/data/rob1g10/SPH/DATA/DISKS/M10_20d_0.0e/disk001',
#                               '/data/rob1g10/SPH/DATA/DISKS/M10_20d_0.0e/disk002',
#                               '/data/rob1g10/SPH/DATA/DISKS/M10_20d_0.0e/disk003'],
#                                       None, 'day', tol=0.01)
#array8 = findCumulativeNumParts(['/data/rob1g10/SPH/DATA/DISKS/M8_20d_0.0e/disk001',
#                                '/data/rob1g10/SPH/DATA/DISKS/M8_20d_0.0e/disk002',
#                                '/data/rob1g10/SPH/DATA/DISKS/M8_20d_0.0e/disk003'],
#                                       None, 'day', tol=0.01)
#array6 = findCumulativeNumParts(['/data/rob1g10/SPH/DATA/DISKS/M6_20d_0.0e/disk001',
#                               '/data/rob1g10/SPH/DATA/DISKS/M6_20d_0.0e/disk002',
#                               '/data/rob1g10/SPH/DATA/DISKS/M6_20d_0.0e/disk003'],
#                                       None, 'day', tol=0.01)                                       
#array1_4 = findCumulativeNumParts(['/data/rob1g10/SPH/DATA/DISKS/M14_20d_0.0e/disk001',
#                                '/data/rob1g10/SPH/DATA/DISKS/M1.4_20d_0.0e/disk002',
#                                '/data/rob1g10/SPH/DATA/DISKS/M1.4_20d_0.0e/disk003'],
#                                       None, 'day', tol=0.01)

                                       

#array14 = createAccretionArray(['/data/rob1g10/SPH/DATA/DISKS/M14_20d_0.0e/disk001',
#                               '/data/rob1g10/SPH/DATA/DISKS/M14_20d_0.0e/disk002',
#                               '/data/rob1g10/SPH/DATA/DISKS/M14_20d_0.0e/disk003',
#                               '/data/rob1g10/SPH/DATA/DISKS/M14_20d_0.0e/disk004'],
#                                       "NS", None, 'day', tol=0.01)
#array10 = createAccretionArray(['/data/rob1g10/SPH/DATA/DISKS/M10_20d_0.0e/disk001',
#                               '/data/rob1g10/SPH/DATA/DISKS/M10_20d_0.0e/disk002',
#                               '/data/rob1g10/SPH/DATA/DISKS/M10_20d_0.0e/disk003',
#                               '/data/rob1g10/SPH/DATA/DISKS/M10_20d_0.0e/disk004'],
#                                       "NS", None, 'day', tol=0.01)
#array8 = createAccretionArray(['/data/rob1g10/SPH/DATA/DISKS/M8_20d_0.0e/disk001',
#                               '/data/rob1g10/SPH/DATA/DISKS/M8_20d_0.0e/disk002',
#                               '/data/rob1g10/SPH/DATA/DISKS/M8_20d_0.0e/disk003',
#                               '/data/rob1g10/SPH/DATA/DISKS/M8_20d_0.0e/disk004'],
#                                       "NS", None, 'day', tol=0.01)
#array6 = createAccretionArray(['/data/rob1g10/SPH/DATA/DISKS/M6_20d_0.0e/disk001',
#                               '/data/rob1g10/SPH/DATA/DISKS/M6_20d_0.0e/disk002',
#                               '/data/rob1g10/SPH/DATA/DISKS/M6_20d_0.0e/disk003',
#                               '/data/rob1g10/SPH/DATA/DISKS/M6_20d_0.0e/disk004'],
#                                       "NS", None, 'day', tol=0.01)
#array1_4 = createAccretionArray(['/data/rob1g10/SPH/DATA/DISKS/M1.4_20d_0.0e/disk001',
#                               '/data/rob1g10/SPH/DATA/DISKS/M1.4_20d_0.0e/disk002',
#                               '/data/rob1g10/SPH/DATA/DISKS/M1.4_20d_0.0e/disk003',
#                               '/data/rob1g10/SPH/DATA/DISKS/M1.4_20d_0.0e/disk004'],
#                                       "NS", None, 'day', tol=0.01)
#                                 
#
#
#plt.plot(range(len(array14)), array14, label='M14')                                     
#plt.plot(range(len(array10)), array10, label='M10')#, color='g')
#plt.plot(range(len(array8)), array8, label='M8')#, color='r')
#plt.plot(range(len(array6)), array6, label='M6')#, color='r')
#plt.plot(range(len(array1_4)), array1_4, label='M1.4')#, color='r')
#plt.xlabel('Days')
###plt.ylabel('Number of particles accreted by Be star') 
#plt.ylabel('Number of particles accreted by compact object')
#plt.legend(loc='best')
#plt.xlim(700, 800)
#plt.xlim(1400, 1600)
#plt.ylim(0, 15)
#Ms = [14, 10, 8, 6, 1.4]
#start = 200
#end = 700 
#norm = (end - start) / 20
#AccPerOrbit = np.array([np.sum(array14[start:end]) / norm, np.sum(array10[start:end]) / norm, 
#                        np.sum(array8[start:end]) / norm, np.sum(array6[start:end]) / norm,
#                        np.sum(array1_4[start:end]) / norm])
#print AccPerOrbit
#plt.plot(Ms, AccPerOrbit)      
#plt.xlabel('Mass of compact object in solar masses')
###plt.ylabel('Number of particles accreted by Be star') 
#plt.ylabel('Number of particles accreted by compact object per orbit')
#for i in range(len(AccPerOrbit)):
#    plt.text(Ms[i]-0.6, max(AccPerOrbit), "{:0.1f}".format(AccPerOrbit[i]))

# BLACK HOLES

filenames=np.array([

#                    'M1.4_20d_0.0e',
#                    'M3_19d_0.0e',   
#                    'M4_18.4d_0.0e', 
#                    'M5_17.9d_0.0e', 
#                    'M6_17.4d_0.0e', 
#                    'M7_17d_0.0e',   
#                    'M8_16.6d_0.0e', 
#                    'M9_16.2d_0.0e',
#                    'M10_15.8d_0.0e',
#                    
#                    'M1.4_23.4d_0.1e', 
#                    'M3_22.2d_0.1e',   
#                    'M4_21.6d_0.1e', 
#                    'M5_21d_0.1e', 
#                    'M6_20.4d_0.1e', 
#                    'M7_19.9d_0.1e',   
#                    'M8_19.4d_0.1e',
#                    'M9_19d_0.1e', 
#                    'M10_18.5d_0.1e', 
#                    
#                    'M1.4_28d_0.2e', 
#                    'M3_26.5d_0.2e', 
#                    'M4_25.7d_0.2e', 
#                    'M5_25d_0.2e',   
#                    'M6_24.3d_0.2e', 
#                    'M7_23.7d_0.2e', 
#                    'M8_23.1d_0.2e',
#                    'M9_22.6d_0.2e',  
#                    'M10_22.1d_0.2e',

#                     "M1.4_34.1d_0.3e",
#                     "M3_32.4d_0.3e",
#                     "M4_31.4d_0.3e",
#                     "M5_30.5d_0.3e",
#                     "M6_29.7d_0.3e",
#                     "M7_29d_0.3e",
#                     "M8_28.3d_0.3e",
#                     "M9_27.6d_0.3e",
#                     "M10_27d_0.3e",       
                    
#                    'M1.4_43d_0.4e', 
#                    'M3_40.8d_0.4e', 
#                    'M4_39.6d_0.4e', 
#                    'M5_38.5d_0.4e', 
#                    'M6_37.5d_0.4e', 
#                    'M7_36.5d_0.4e', 
#                    'M8_35.6d_0.4e', 
#                    'M9_34.8d_0.4e',
#                    'M10_34.1d_0.4e',
#                    'M10_34.1d_0.4e_10000',
                    
#                    'M1.4_56.6d_0.5e', 
#                    'M3_53.7d_0.5e', 
#                    'M4_52.1d_0.5e', 
#                    'M5_50.6d_0.5e', 
#                    'M6_49.2d_0.5e', 
#                    'M7_48d_0.5e', 
#                    'M8_46.8d_0.5e', 
#                    'M9_45.8d_0.5e', 
#                    'M10_44.8d_0.5e',
                    
#                    'M1.4_79.1d_0.6e', 
#                    'M3_75d_0.6e', 
#                    'M4_72.8d_0.6e', 
#                    'M5_70.7d_0.6e', 
#                    'M6_68.8d_0.6e', 
#                    'M7_67.1d_0.6e', 
#                    'M8_65.5d_0.6e', 
#                    'M9_64d_0.6e', 
#                    'M10_62.6d_0.6e',
                    
#                    "M20_13.2d_0.0e", 
#                    "M20_18.5d_0.2e", 
#                    "M20_28.4d_0.4e",
                    
#                    "SXP5.05_M-11_alpha0.1",
#                    "SXP5.05_M-10_alpha0.1",
#                    "SXP5.05_M-9_alpha0.1",
#                    
#                    "SXP5.05_M-11_alpha0.3",
#                    "SXP5.05_M-10_alpha0.3",
#                    "SXP5.05_M-9_alpha0.3",
#                    
#                    "SXP5.05_M-11_alpha1",
#                    "SXP5.05_M-10_alpha1",
#                    "SXP5.05_M-9_alpha1",
                    
#                    "MWC656",                    
                    
#                    'isoBe_M15_Mdot9',
#                    'isoBe_M15_Mdot10',
#                    'isoBe_M15_Mdot11',
                    ])
    
discSizeByEye = np.array([7.6, 5.0, 4.2,
                     8.0, 6.5, 6.1,
                     10.2, 6.7, 6.6])
discSizeBy13cutoff = np.array([9.89, 6.86, 6.57,
                               10.28, 8.82, 6.41,
                               10.96, 8.05, 7.71])
                               
#discSizeBy13cutoffAve = np.array([9.34496185141, 6.43417567368, 6.0038530793,
#                               9.87757604912, 6.48824271039, 6.05448715052,
#                               12.9974279454, 6.54423984746, 6.31058968542])
                               
discSizeBy13cutoffAve = np.array([9.46188771845, 7.58371367007, 6.75718735671, 6.41041726437, 6.25431580103, 6.15325731571, 6.26965414598, 5.64591948605,
                                  9.71139226483, 7.71351937642, 7.25491756203, 7.0506396408, 6.43728448136, 6.57165240994, 6.53867961734, 6.02314874814,
                                  13.1382444591, 11.811350156, 8.56327253211, 7.728423661, 10.5713026043, 14.3949363267, 7.17138425154, 5.96694752017])

discSizeBy12cutoffAve = np.array([7.22221409951, 5.97505469443, 5.79970553491, 5.20237481989, 5.21161958273, 5.03472727644, 5.48019605986, 5.06675210707,
                                  7.42376306644, 6.86082955432, 6.59115691822, 6.41280700028, 6.07888515902, 5.7951635589, 5.5921994887, 5.50082631155,
                                  9.59554182492, 8.72651748796, 7.63723231198, 6.72979852311, 9.86385838692, 13.625203094, 6.41633509154, 5.5566450489])

orbit = np.array([])
Mx = np.array([])
#for BN in filenames:
#    underscores = np.array([])
#    index = 0
#    for char in BN:
##        print char, index
#        if(char == 'M'):
#            Mplace = index
#        if(char == 'd'):
#            dplace = index
#        if(char == '_'):
#            underscores = np.append(underscores, index)
#        index += 1
#    underscore1 = int(np.min(underscores))
##    print underscore1
#    orbit = np.append(orbit, float(BN[underscore1+1:dplace]) )
#    Mx = np.append(Mx, float(BN[Mplace+1:underscore1]) )

#plt.scatter(Mx[0:8]/2e30, discSizeBy13cutoffAve[0:8], marker='o', color='r', label="cutoff=10^-13")
#plt.scatter(Mx[0:8]/2e30, discSizeBy12cutoffAve[0:8], marker='o', color='b', label="cutoff=10^-12")
#plt.scatter(Mx[8:16]/2e30, discSizeBy13cutoffAve[8:16], marker='x', color='r', label="cutoff=10^-13")
#plt.scatter(Mx[8:16]/2e30, discSizeBy12cutoffAve[8:16], marker='x', color='b', label="cutoff=10^-12")
#plt.scatter(Mx[16:24]/2e30, discSizeBy13cutoffAve[16:24], marker='^', color='r', label="cutoff=10^-13")
#plt.scatter(Mx[16:24]/2e30, discSizeBy12cutoffAve[16:24], marker='^', color='b', label="cutoff=10^-12")
                    
#Rx = np.array([10e3, 17.8e3, 29.6e3,
#               10e3, 17.8e3, 29.6e3,
#               10e3, 17.8e3, 29.6e3,
#               29.6e3])
               
#Orbitx = np.array([20,   19, 18.4, 17.9, 17.4,   17, 16.6, 15.8,
#                   28, 26.5, 25.7,   25, 24.3, 23.7, 23.1, 22.1,
#                   43, 40.8, 39.6, 38.5, 37.5, 36.5, 35.6, 34.1,
#                   34.1])
Orbitx = orbit                 
#orbit = np.array([20,   19, 18.4, 17.9, 17.4,   17, 16.6, 15.8,
#                  28, 26.5, 25.7,   25, 24.3, 23.7, 23.1, 22.1,
#                  43, 40.8, 39.6, 38.5, 37.5, 36.5, 35.6, 34.1])
bPhases04e = np.array([0.11618648, 0.5976398, 0.38550201, 0.83917929, 0.93665309, 0.75851662, 0.83093399, 0.91999362])
         
arrayMax = np.array([])
mcapturePhase = np.array([])

arrayAverage = np.array([])

starts = np.array([])
ends = np.array([])

sumParts = True

count=0              
for j in range(0, len(filenames)):
    filename = '/data/rob1g10/SPH/DATA/DISKS/{}/disk001'.format(filenames[j])
#    print filename
    # find units

    phrase1 = ["following", "units", "units", "of", ":",  "mass", ":"]
    phrase2 = ["following", "units", "units", "of", ":", "distance", ":"]
    phrase3 = ["following", "units", "units", "of", ":", "time", ":" ]

    """ Strip disk """   
    
    PHRASES = [phrase1, phrase2, phrase3]
    values = np.array([])
    
    for i in range(3):
        beforePhrases = PHRASES[i]
        f = open(filename, 'r')
    
        beforePhraseslen = len(beforePhrases)
    
        count = 0
        phraseCount = 0
    
        """ Without orbit checking """
        for line in f:
            strippedLine = line.strip()
            columns = strippedLine.split()
            for i in range(len(columns)):
                if(phraseCount == beforePhraseslen):
                    values = np.append(values, columns[i])
                    phraseCount = 0
                    count += 1
                elif (columns[i] == beforePhrases[phraseCount]):
                    phraseCount += 1
            if(count == 1):
                break
        
        f.close()
        
    umass = float(values[0][0:6])*10**float(values[0][7:10])
    udist = float(values[1])
    utime = float(values[2])    

    # plot 
    if(sumParts==True):
        filename = '/data/rob1g10/SPH/DATA/DISKS/'+filenames[j]
        
        os.system("ls /data/rob1g10/SPH/DATA/DISKS/{}/disk* > lsdisks.txt".format(filenames[j]))
    
        f = open("lsdisks.txt", 'r')
        lines = f.readlines()
    	
        numdisks = int(lines[-1][-3:])
        
        f.close()

        disks = np.array([])        
        
        for ND in range(1, numdisks+1):
            disks = np.append(disks, '/data/rob1g10/SPH/DATA/DISKS/{}/disk{:03}'.format(filenames[j], ND))
        

        array = findCumulativeNumParts(disks, None, 'day', tol=0.01)
                                    
        plt.plot(range(len(array)), array, label=filenames[j])
        plt.legend(loc='best')
        plt.xlabel("Days")
        plt.ylabel("Number of Particles")
        print filenames[j]
    
    
    
    
    
    
    else:
        filename = '/data/rob1g10/SPH/DATA/DISKS/'+filenames[j]
        
        os.system("ls /data/rob1g10/SPH/DATA/DISKS/{}/disk* > lsdisks.txt".format(filenames[j]))
    
        f = open("lsdisks.txt", 'r')
        lines = f.readlines()
    	
        numdisks = int(lines[-1][-3:])
        
        f.close()

        disks = np.array([])        
        
        for ND in range(1, numdisks+1):
            disks = np.append(disks, '/data/rob1g10/SPH/DATA/DISKS/{}/disk{:03}'.format(filenames[j], ND))
        
        times, array = createAccretionArray(disks, "NS", None, 'day', tol=0.01)#*1e-15 * umass / (24*3600)
        
#        plt.plot(range(len(array)), array)
    #    print(len((np.arange(len(array))-2000)/Orbitx[i]), len(6.67e-11*Mx[i]*array / Rx[i]))
        arrayMax = np.append(arrayMax, np.max(array))
#        periastrons = np.array([10604.8956481,10876.6435185,9724.82893519,19263.4161806,10703.5259722,22310.5965972,
#                                14831.2488889,9698.40474537])
#        mcapturePhase = np.append(mcapturePhase, ((np.argmax(array[-201:-1])-periastrons[j-16]) % Orbitx[j])/Orbitx[j] )
        print len(array)
        # find range of 10 orbits for the system
        kend = int( len(array) / int(Orbitx[j]) )
#        if(kend > 50): 
#            kend = 50
        print kend
        end = -1 #5*int(Orbitx[j]) 
        arrayAve = array[end-int(Orbitx[j]):end]
        for k in range(kend-40, kend):
            start = end - int(Orbitx[j])
            arrayAve += array[start:end]
            end -= int(Orbitx[j])  
        print arrayAve
        
        for index in range(len(arrayAve)):
            if(arrayAve[index] <= 1e13):
                arrayAve[index] = 1e15
        print arrayAve
        Mx = 2e30*np.array([1.4, 3., 4., 5., 6., 7., 8., 10.,
                    1.4, 3., 4., 5., 6., 7., 8., 10.,
                    1.4, 3., 4., 5., 6., 7., 8., 10.,
                    10.])
        print Mx, 
        print j-16, bPhases04e[j-16],#, Rx
        arrayAverage = np.append(arrayAverage, np.mean(arrayAve))
#        xTIME = np.append(np.arange(len(arrayAve)), Orbitx[j]+np.arange(len(arrayAve)))/Orbitx[j]
#        plt.plot( xTIME, 
#                 6.67e-11*np.append(arrayAve, arrayAve), label="M = "+str(Mx[j]/2e30),
#                color=['b', 'g', 'r', 'y', 'black', 'brown', 'turquoise', 'pink'][j-16])
#        plt.plot(np.ones(100)*bPhases04e[int(j-16)], np.linspace(0, 6.67e-11*np.max(arrayAve), 100),
#                 color=['b', 'g', 'r', 'y', 'black', 'brown', 'turquoise', 'pink'][j-16])
#        plt.plot(np.ones(100)*bPhases04e[int(j-16)]+1, np.linspace(0, 6.67e-11*np.max(arrayAve), 100),
#                 color=['b', 'g', 'r', 'y', 'black', 'brown', 'turquoise', 'pink'][j-16])
#        count += 1
#Mx = np.array([1.4, 3., 4., 5., 6., 7., 8., 10.,
#                    1.4, 3., 4., 5., 6., 7., 8., 10.,
#                    1.4, 3., 4., 5., 6., 7., 8., 10.,
#                    10.])

#mcapturePhase = np.array([ 0.1627907 ,  0.17647059,  0.13636364,  0.49350649,  0.22666667, 0.57534247,  0.32022472,  0.46627566])
bPhases04ePre = [1.16186477e-01,   4.88907044e-01,   3.85502011e-01,   1.69191919e-02,
                 3.14814815e-04,   1.31640030e-02,   3.85295724e-01,  1.54868578e-02]
bPhases04ePost = [0.98886305,  0.5976398,   0.51282618,  0.83917929,  0.93665309,  0.75851662,
                  0.83093399,  0.91999362]
Age = np.array([10628.1640046, 10875.3448148, 9727.3499537, 19289.1513889, 10732.17375, 22354.746088,
                14839.18125, 9743.99321759])
periastrons = np.array([10604.8956481,10876.6435185,9724.82893519,19263.4161806,10703.5259722,22310.5965972,
                        14831.2488889,9698.40474537])
#bPhases04e = np.array([ 0.16660476,  0.5976398-0.5,   0.38550201,  0.83917929-0.5,  0.84116049-0.5,  0.75851662, 0.83093399-0.5,  0.82450038-0.5])
#plt.scatter(Mx[16:24], bPhases04ePre, label='Min Relative Velocity Pre Per')
#plt.scatter(Mx[16:24], bPhases04ePost, color='g', label='Min Relative Velocity Post Per')
#plt.scatter(Mx[16:24], mcapturePhase, color='r', label='Peak mass capture')
#plt.legend(loc='best')

#if(sumParts==False):
#arrayMax *= 6.67e-11*Mx[j] / Rx[j]
#arrayAverage *= 6.67e-11*Mx[j] / Rx[j]
#arrayAverage[1] = 10
#arrayAverage[7] = 10
#plt.scatter(Mx[0:9], arrayMax[0:9], marker='o', color='b', label="e = 0.0")
#plt.scatter(Mx[9:18], arrayMax[9:18], marker='x', color='r', label="e = 0.1")
#plt.scatter(Mx[18:27], arrayMax[18:27], marker='^', color='g', label="e = 0.2")
#plt.scatter(Mx[27:36], arrayMax[27:36], marker='H', color='black', label="e = 0.4")
#plt.scatter(Mx[36:45], arrayMax[36:45], marker='<', color='purple', label="e = 0.6")
#plt.legend(loc='best')
#plt.xlabel("Compact object mass in solar masses")
##plt.ylabel("Disc size (by specific cut off) in stellar radii")
#plt.ylabel("log Maximum X-ray Luminosity in ergs / s")
#plt.xlabel("Binary phase")
#plt.ylabel("Average X-ray Luminosity in ergs / s")
###plt.xlim(0, 5)
#plt.ylim(0, 10)
                