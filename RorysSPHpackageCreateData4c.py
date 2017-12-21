# -*- coding: utf-8 -*-
"""
Created on Mon Aug 22 11:15:09 2016

@author: rob1g10
"""

import numpy as np
import os
#from RorysSPHpackageSpectra import plotAndSaveSpectralLine, calcEW

"""----------------------------------------------------------------------------

                PLOT OPTICAL DEPTH

----------------------------------------------------------------------------"""

def createOptDepthData4c(CcodeName, partcoords, partvels, u, rhos, numparts, partmass, BeRadius, BeMass,
                         powerLawIndex, baseDensity, theta, phi):
    
    T = 2e5
    m_e = 9.1e-28 # mass of electron
    e = 4.8e-10 # electron charge
    lambda0 = 656.28e-7 # rest wavelength
    c = 3e10 # speed of light
    nu0 = c / lambda0 # rest frequency
    f = 0.6407 # absorption oscillator strength of the transition. for n=3 to n=2 Balmer line (emitting at 656nm), this is 0.6407
    kb = 1.38e-16 # Boltzmann constant 
    h = 6.626e-27 # plank constant 
    amu = 1.66e-24 # atomic mass unit
    G = 6.67e-8 # gravitational constant
    M_stellar = 18 * 2e33;
    
    """
    inclination angle in degrees
    """
    theta1 = theta
    phi1 = phi
    theta = np.deg2rad(theta)
    phi = np.deg2rad(phi)
    
    xvelsInit = partvels[:,0]
    yvelsInit = partvels[:,1]
    zvelsInit = partvels[:,2]
    
    """ rotate coords and vels """
    
    xcoords = partcoords[:,0] * np.cos(phi) - partcoords[:,1] * np.sin(phi)
    ycoords = (partcoords[:,0] * np.sin(phi) + partcoords[:,1] * np.cos(phi))*np.cos(theta) - partcoords[:,2] * np.sin(theta)
    zcoords = partcoords[:,1] * np.sin(theta) + partcoords[:,2] * np.cos(theta)
    print min(zcoords) / (7 * 6.96e10), max(zcoords) / (7 * 6.96e10)
    xvels = xvelsInit * np.cos(phi) - yvelsInit * np.sin(phi)
    yvels = (xvelsInit * np.sin(phi) + yvelsInit * np.cos(phi)) * np.cos(theta) - zvelsInit * np.sin(theta)
    zvels = yvelsInit * np.sin(theta) + zvelsInit * np.cos(theta)  
    
#    """ sort arrays """
#    
    sortOrder = np.argsort(zcoords)[::-1]
#    print sortOrder
#    
#    xcoords = xcoords[sortOrder]
#    ycoords = ycoords[sortOrder]
#    zcoords = zcoords[sortOrder]
#    xvels = xvels[sortOrder]
#    yvels = yvels[sortOrder]
#    zvels = zvels[sortOrder]
    
    """ calc ionFrac """
    R = np.sqrt(xcoords**2 + ycoords**2)
    T = T / (R / (7 * 6.96e10))
    ionCoeff = (2 * np.pi * m_e * kb * T) / (h*h)
    ionExp = np.exp( (-13.6 * 1.6e-12) / (kb * T))
    ionNumDens = 6.022e23 * rhos
    ionA = 1.
    ionB = (1. / ionNumDens) * ionCoeff**1.5 * ionExp;
    ionC = -ionB;
    ionFrac = (-ionB + np.sqrt(ionB*ionB - 4*ionA*ionC)) / (2*ionA);
    
    f = open("INPUT/optDepthX4c{}.txt".format(CcodeName), "w")
    for i in range(len(partcoords)):
        f.write("{}\n".format(partcoords[i,0]))
    f.close()
    f = open("INPUT/optDepthY4c{}.txt".format(CcodeName), "w")
    for i in range(len(partcoords)):
        f.write("{}\n".format(partcoords[i,1]))
    f.close()
    f = open("INPUT/optDepthZ4c{}.txt".format(CcodeName), "w")
    for i in range(len(partcoords)):
        f.write("{}\n".format(partcoords[i,2]))
    f.close()
    f = open("INPUT/optDepthXrot4c{}.txt".format(CcodeName), "w")
    for i in range(len(partcoords)):
        f.write("{}\n".format(xcoords[i]))
    f.close()
    f = open("INPUT/optDepthYrot4c{}.txt".format(CcodeName), "w")
    for i in range(len(partcoords)):
        f.write("{}\n".format(ycoords[i]))
    f.close()
    f = open("INPUT/optDepthZrot4c{}.txt".format(CcodeName), "w")
    for i in range(len(partcoords)):
        f.write("{}\n".format(zcoords[i]))
    f.close()
    f = open("INPUT/optDepthVX4c{}.txt".format(CcodeName), "w")
    for i in range(len(partcoords)):
        f.write("{}\n".format(xvels[i]))
    f.close()
    f = open("INPUT/optDepthVY4c{}.txt".format(CcodeName), "w")
    for i in range(len(partcoords)):
        f.write("{}\n".format(yvels[i]))
    f.close()
    f = open("INPUT/optDepthVZ4c{}.txt".format(CcodeName), "w")
    for i in range(len(partcoords)):
        f.write("{}\n".format(zvels[i]))
    f.close()
#    f = open("optDepthVels4c.txt", "w")
#    for i in range(len(lineOfSightDists)):
#        f.write("{}\n".format(partvels[i,0]))
#        f.write("{}\n".format(partvels[i,1]))
#        f.write("{}\n".format(partvels[i,2]))
#    f.close()
    f = open("INPUT/optDepthDens4c{}.txt".format(CcodeName), "w")
    for i in range(len(rhos)):
        f.write("{:.10}\n".format(rhos[i]))
    f.close()
    f = open("INPUT/optDepthIntergs4c{}.txt".format(CcodeName), "w")
    for i in range(len(u)):
        f.write("{:.10}\n".format(u[i]))
    f.close()
    f = open("INPUT/optDepthAzimuth4c{}.txt".format(CcodeName), "w")
    for i in range(len(xcoords)):
        f.write("{:.10}\n".format(np.tan(xcoords[i]/ycoords[i])))
    f.close()
    f = open("INPUT/optDepthSort4c{}.txt".format(CcodeName), "w")
    for i in range(len(xcoords)):
        f.write("{}\n".format(sortOrder[i]))
    f.close()
    f = open("INPUT/optDepthIonFrac4c{}.txt".format(CcodeName), "w")
    for i in range(len(xcoords)):
        if(1 - ionFrac[i] < 0):
            f.write("{}\n".format(0.0))
#            print i, 1-ionFrac[i], - ionB[i] + np.sqrt(ionB*ionB - 4*ionA*ionC)[i], rhos[i]
        else:
            f.write("{}\n".format(1 - ionFrac[i]))  
    f.close()
    
    f = open("INPUT/optDepthParams4c{}.txt".format(CcodeName), "w")
    f.write("{}\n".format(numparts))
    f.write("{:.10}\n".format(partmass))
    f.write("{:.10}\n".format(BeRadius))
    f.write("{:.10}\n".format(BeMass))
    f.write("{:.10}\n".format(powerLawIndex))
    f.write("{:.10}\n".format(baseDensity))
    f.write("{:.10}\n".format(theta))
    f.write("{:.10}\n".format(phi))
    f.write("{:.10}\n".format(CcodeName))
    f.close()
    return 1-ionFrac


def readOptDepthData4py(filename):
    optDepth = np.array([])
    " Open file "
    f = open(filename, 'r')
    
    for line in f:
        line = line.strip()
        optDepth = np.append(optDepth, float(line))

    f.close()

    return optDepth

def calcOptDepthsOverTime(CcodeName, inputName, baseName, inputdirectory, startDISKs, 
                          numDISKs, phi, theta, SubQDISK, dropbox='n'):
    """
    Produces an animation using pngs created and placed in a given folder.


    INPUT VARIABLES
    ---------------

    inDirs (list of strings) - directories where the QGDISK and IPDISK files 
                               reside.
    
    utime, orbitPeriod (floats) - unit of time and orbit period for simulation.

    OUTPUT VARIABLES
    ----------------

    An animation of the respective files.

    """   
        
    os.system("rm lsQG.txt lsIP.txt")
    
    QGfilenames = np.array([])
    
    f = open("lsQG.txt", "w")
    for j in range(startDISKs, startDISKs + numDISKs):
        for k in range(SubQDISK[0], SubQDISK[1]):
            NAME = "{}QGDISK{:03}{:02}\n".format(inputdirectory, j, k)
            if(os.system("ls {}".format(NAME)) == 0):
                f.write(NAME)
            else:
                break
    f.close()
    
    f = open("lsIP.txt", "w")
    for j in range(startDISKs, startDISKs + numDISKs):
        NAME = "{}IPDISK{:03}\n".format(inputdirectory, j)
        f.write(NAME)
    f.close()
    
    """ QGDISKs """    
    
    f = open("lsQG.txt", "r")
    
    lines = f.readlines()

    for line in lines:
        f2 = open("{}".format(line[0:-1]))
        temp = f2.readlines()
        f2.close()
        if(len(temp) != 0):
            QGfilenames = np.append(QGfilenames, line[0:-1])
        
    f.close()
        
    for i in range(len(QGfilenames)):
        filename = QGfilenames[i]
        diskname = filename[:-26] + "disk001"
        IPname = filename[:-11] + "IPDISK" + filename[-5:-2]
        Becoords, NScoords, partcoords, partvels, rho, u, num_parts, time, umass, udist, utime = readFILES(filename, diskname, IPname)
        ionFrac = createOptDepthData4c(inputName, #CcodeName,
                             partcoords * udist, partvels * udist/utime, 
                             u * udist**2/utime**2, 
                             rho * umass / udist**3, num_parts, 
                             1e-15 * umass / 2e30, 
                             (7 * 6.96e10), 18 * 2e33, 2.5, 
                             max(rho) * umass / udist**3, theta, phi)
        os.system("./opticalDepth{}".format(CcodeName))
        name = "{}_theta{}_phi{}_time{:0.1f}".format(baseName, theta, phi, time * utime / (24 * 3600))
        if(dropbox == 'y'):        
            os.system("cp optDepthData4py{}.txt /data/rob1g10/SPH/DATA/DISKS/{}/optDepthData{}{}.txt".format(CcodeName, baseName, CcodeName, name))
            os.system("mv optDepthData4py{}.txt DATA/optDepthData{}{}.txt".format(CcodeName, CcodeName, name))
        else:
            os.system("mv optDepthData4py{}.txt /data/rob1g10/SPH/DATA/DISKS/{}/optDepthData{}{}.txt".format(CcodeName, baseName, CcodeName, name))
#        print time * utime / (24 *3600)
        print name, "done"
    os.system("rm lsQG.txt lsIP.txt")
    
    return ionFrac

"""----------------------------------------------------------------------------

                FILE READING

----------------------------------------------------------------------------"""


""" ASCII reading """

def readQDISK(filename):

    " Open file "
    f = open(filename, 'r')

    " calculate number of particles and set arrays to be big enough for the data "
    num_parts = 0
    for line in f:
        #if(num_parts == 2):
        #    print line
        num_parts += 1

    f.close()

    partcoords = np.zeros((num_parts, 3))
    partvels = np.zeros((num_parts, 3))
    xcoords = np.zeros(num_parts)
    ycoords = np.zeros(num_parts)
    zcoords = np.zeros(num_parts)
    xvels = np.zeros(num_parts)
    yvels = np.zeros(num_parts)
    zvels = np.zeros(num_parts)
    rho = np.zeros(num_parts)
    u = np.zeros(num_parts)


    f = open(filename, 'r')
    " Loop over lines and extract variables of interest "
    count = 0
    for line in f:
        line = line.strip()
        columns = line.split()
        partcoords[count] = [float(columns[0]), float(columns[1]), float(columns[2])]
        xcoords[count] = float(columns[0])
        ycoords[count] = float(columns[1])
        zcoords[count] = float(columns[2])
        partvels[count] = [float(columns[3]), float(columns[4]), float(columns[5])]
        xvels[count] = float(columns[3])
        yvels[count] = float(columns[4])
        zvels[count] = float(columns[5])
        rho[count] = float(columns[6])
        u[count] = float(columns[9])
        count += 1

    f.close()

    Becoords = np.array([xcoords[0], ycoords[0], zcoords[0]])
    NScoords = np.array([xcoords[1], ycoords[1], zcoords[1]])

    partcoords -= Becoords
    Becoords -= Becoords
    NScoords -= Becoords

    return Becoords, NScoords, partcoords[2:], partvels[2:], rho[2:], u[2:], num_parts-2

def readIPDISK(filename):# utime, orbitPeriod):

    time = np.array([])
    
    " Open file "
    f = open(filename, 'r')
    
    for line in f:
        line = line.strip()
        columns = line.split()
        time = np.append(time, float(columns[0]))# * utime / (24 * 3600 * orbitPeriod))

    f.close()

    return time

def readUnits(filename):
    """
    Goes through a disk file and extracts the units of the simulation.

    INPUT VARIABLES

    filename (str) - name of disk file to extract data from.

    tol (float) - tolerance for code to allow for matching of desired orbits

    OUTPUT VARIABLES

    values - an array giving the extracted data

    """
    
    """ Set up phrases for search """
    
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
    
    return umass, udist, utime
    
def readBeRadius(filename):
    return
    
def readFILES(filename, diskname, IPname):
    Becoords, NScoords, partcoords, partvels, rho, u, num_parts = readQDISK(filename)
    
    umass, udist, utime = readUnits(diskname)
    
    #timeArray = readIPDISK(filename[0:-11]+"IPDISK"+filename[-5:-2])
    timeArray = readIPDISK(IPname)
#    try:
    time = timeArray[int(filename[-2:]) - 1]
    
    
    return Becoords, NScoords, partcoords, partvels, rho, u, num_parts, time, umass, udist, utime
#    except:
#        return Becoords, NScoords, partcoords, partvels, rho, u, num_parts, umass, udist, utime


"""----------------------------------------------------------------------------

                Running

----------------------------------------------------------------------------"""

filename = 'INPUT/QGDISK41101'
diskname = 'INPUT/disk001'
IPname = 'INPUT/IPDISK411'
Becoords, NScoords, partcoords, partvels, rho, u, num_parts, time, umass, udist, utime = readFILES(filename, diskname, IPname)
phi = 0
#for theta in range(0, 91, 10):
createOptDepthData4c("BASIC",
                     partcoords * udist, partvels * udist/utime, 
                     u * udist**2/utime**2, 
                     rho * umass / udist**3, num_parts, 
                     1e-15 * umass / 2e30, 
                     (7 * 6.96e10), 18 * 2e33, 2.5, 
                     max(rho) * umass / udist**3, theta, phi)