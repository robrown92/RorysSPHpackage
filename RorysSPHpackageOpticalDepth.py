# -*- coding: utf-8 -*-
"""
Created on Mon Aug 22 11:15:09 2016

@author: rob1g10
"""

import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import pylab as plt
import os
from scipy.optimize import leastsq, curve_fit
from RorysSPHpackage import readQDISK, readUnits, readFILES, plotAndSaveFig
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
    theta = np.deg2rad(theta)
    phi = np.deg2rad(phi)
    
    xvelsInit = partvels[:,0]
    yvelsInit = partvels[:,1]
    zvelsInit = partvels[:,2]
    
    """ rotate coords and vels """
    
    xcoords = partcoords[:,0] * np.cos(phi) - partcoords[:,1] * np.sin(phi)
    ycoords = (partcoords[:,0] * np.sin(phi) + partcoords[:,1] * np.cos(phi))*np.cos(theta) - partcoords[:,2] * np.sin(theta)
    zcoords = partcoords[:,1] * np.sin(theta) + partcoords[:,2] * np.cos(theta)
    print min(zcoords), max(zcoords)
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
    
    f = open("optDepthX4c{}.txt".format(CcodeName), "w")
    for i in range(len(partcoords)):
        f.write("{}\n".format(partcoords[i,0]))
    f.close()
    f = open("optDepthY4c{}.txt".format(CcodeName), "w")
    for i in range(len(partcoords)):
        f.write("{}\n".format(partcoords[i,1]))
    f.close()
    f = open("optDepthZ4c{}.txt".format(CcodeName), "w")
    for i in range(len(partcoords)):
        f.write("{}\n".format(partcoords[i,2]))
    f.close()
    f = open("optDepthXrot4c{}.txt".format(CcodeName), "w")
    for i in range(len(partcoords)):
        f.write("{}\n".format(xcoords[i]))
    f.close()
    f = open("optDepthYrot4c{}.txt".format(CcodeName), "w")
    for i in range(len(partcoords)):
        f.write("{}\n".format(ycoords[i]))
    f.close()
    f = open("optDepthZrot4c{}.txt".format(CcodeName), "w")
    for i in range(len(partcoords)):
        f.write("{}\n".format(zcoords[i]))
    f.close()
    f = open("optDepthVX4c{}.txt".format(CcodeName), "w")
    for i in range(len(partcoords)):
        f.write("{}\n".format(xvels[i]))
    f.close()
    f = open("optDepthVY4c{}.txt".format(CcodeName), "w")
    for i in range(len(partcoords)):
        f.write("{}\n".format(yvels[i]))
    f.close()
    f = open("optDepthVZ4c{}.txt".format(CcodeName), "w")
    for i in range(len(partcoords)):
        f.write("{}\n".format(zvels[i]))
    f.close()
#    f = open("optDepthVels4c.txt", "w")
#    for i in range(len(lineOfSightDists)):
#        f.write("{}\n".format(partvels[i,0]))
#        f.write("{}\n".format(partvels[i,1]))
#        f.write("{}\n".format(partvels[i,2]))
#    f.close()
    f = open("optDepthDens4c{}.txt".format(CcodeName), "w")
    for i in range(len(rhos)):
        f.write("{:.10}\n".format(rhos[i]))
    f.close()
    f = open("optDepthIntergs4c{}.txt".format(CcodeName), "w")
    for i in range(len(u)):
        f.write("{:.10}\n".format(u[i]))
    f.close()
    f = open("optDepthAzimuth4c{}.txt".format(CcodeName), "w")
    for i in range(len(xcoords)):
        f.write("{:.10}\n".format(np.tan(xcoords[i]/ycoords[i])))
    f.close()
    f = open("optDepthSort4c{}.txt".format(CcodeName), "w")
    for i in range(len(xcoords)):
        f.write("{}\n".format(sortOrder[i]))
    f.close()
    f = open("optDepthIonFrac4c{}.txt".format(CcodeName), "w")
    for i in range(len(xcoords)):
        if(1 - ionFrac[i] < 0):
            f.write("{}\n".format(0.0))
#            print i, 1-ionFrac[i], - ionB[i] + np.sqrt(ionB*ionB - 4*ionA*ionC)[i], rhos[i]
        else:
            f.write("{}\n".format(1 - ionFrac[i]))  
    f.close()
    
    f = open("optDepthParams4c{}.txt".format(CcodeName), "w")
    f.write("{}\n".format(numparts))
    f.write("{:.10}\n".format(partmass))
    f.write("{:.10}\n".format(BeRadius))
    f.write("{:.10}\n".format(BeMass))
    f.write("{:.10}\n".format(powerLawIndex))
    f.write("{:.10}\n".format(baseDensity))
    f.write("{:.10}\n".format(theta))
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
    

def OPTICALDEPTHMAIN(partcoords, patvels, u, num_parts, udist, utime, umass, outFilename):
    " create input "
    createOptDepthData4c(partcoords * (udist*0.01), partvels * (udist * 0.01)/utime, u * (udist * 0.01)**2/utime**2, num_parts, 
                     1e-15 * 0.001 * umass / 2e30, 0, (7 * 6.96e8), 2.5)
    " run calculation "
    os.system("./opticalDepth")
    " move output "
    os.system("mv optDepthData4py.txt {}.txt".format(outFilename))
    

def calcRotatingOptDepths(filenames, baseName, thetaRange, phiRange):
    """
    Produces a scatter plots of varying viewing angles of a system
    in 2D of the particles and stars.

    GEOMETRY
    --------

    z axis is height of disk, y axis is width and x axis is depth.

    INPUT VARIABLES
    ---------------
    
    filename (list of 3 strings) - name of the filename with the particle 
                                   data on it and name of disk file for 
                                   units and also IP disk file name.

    pngdirectory (string) - directory where png files are found.

    baseName (string) - base of the name of the file that the plot will be saved as.

    OUTPUT VARIABLES
    ----------------

    Optical Depth values for all particles in given file from all angles..

    """
    
#    for theta in range(0, 181, 15):
#        for phi in range(0, 181, 15):
    for theta in range(thetaRange[0], thetaRange[1], thetaRange[2]):
        for phi in range(phiRange[0], phiRange[1], phiRange[2]):
            Becoords, NScoords, partcoords, partvels, rho, u, num_parts, time, umass, udist, utime = readFILES(filenames[0], filenames[1], filenames[2])
            createOptDepthData4c(partcoords * (udist*0.01), partvels * (udist * 0.01)/utime, u * (udist * 0.01)**2/utime**2, num_parts, 
                                 1e-15 * 0.001 * umass / 2e30, (7 * 6.96e10), 2.5, 
                                 max(rho) * umass * 1e6 / udist**3, theta, phi)
            print max(rho) * umass * 1e6 / udist**3
            os.system("./opticalDepth")
            name = "{}_theta{}_phi{}_time{:0.1f}".format(baseName, theta, phi, time * utime / (24 * 3600))
            
#            os.system("cp optDepthData4py.txt /home/rob1g10/Dropbox/PhD/DISKs/optDepthData{}.txt".format(name))
            os.system("mv optDepthData4py.txt /data/rob1g10/SPH/DATA/DISKS/{}/optDepthData{}.txt".format(baseName, name))
            print(theta, phi)
    return
    
#def calcOptDepthsOverTime(inputdirectory, baseName, theta, phi, optDepthCode="C"):
#    """
#    Produces a scatter plots of varying viewing angles of a system
#    in 2D of the particles and stars.
#
#    GEOMETRY
#    --------
#
#    z axis is height of disk, y axis is width and x axis is depth.
#
#    INPUT VARIABLES
#    ---------------
#    
#    
#
#    OUTPUT VARIABLES
#    ----------------
#
#
#    """
#    
#    QGfilenames = np.array([])
#    
#    os.system("ls {}QG* > lsQG.txt".format(inputdirectory))
#    
#    umass, udist, utime = readUnits("{}disk001".format(inputdirectory))
#    
#    """ QGDISKs """    
#    
#    f = open("lsQG.txt", "r")
#    
#    lines = f.readlines()
#
#    for line in lines:
#        f2 = open("{}".format(line[0:-1]))
#        temp = f2.readlines()
#        f2.close()
#        if(len(temp) != 0):
#            QGfilenames = np.append(QGfilenames, line[0:-1])
#        
#    f.close() 
#    
#    os.system("rm lsQG.txt".format(inputdirectory))
#    
#    if(optDepthCode=="C"):    
#        for i in range(len(QGfilenames)):
#            try:
#                Becoords, NScoords, partcoords, partvels, rho, u, num_parts, time, umass, udist, utime = readFILES(QGfilenames[i], 
#                                                                                                         "{}disk001".format(inputdirectory),
#                                                                                                         "{}IPDISK{}".format(inputdirectory, QGfilenames[i][-5:-2]))
#                createOptDepthData4c(partcoords * (udist*0.01), partvels * (udist * 0.01)/utime, u * (udist * 0.01)**2/utime**2, num_parts, 
#                                     1e-15 * 0.001 * umass / 2e30, (7 * 6.96e10), 2.5, theta, phi)
#                print(time * utime / (24 * 3600))
#                os.system("./opticalDepthTime")
#                name = "{}_theta{}_phi{}_time{:0.1f}".format(baseName, theta, phi, time * utime / (24 * 3600))
#                os.system("mv optDepthData4pyTime.txt /data/rob1g10/SPH/DATA/DISKS/{}/optDepthData{}.txt".format(baseName, name))
#            except:
#                pass 
#    elif(optDepthCode=="Python"):
#        for i in range(len(QGfilenames)):
#            Becoords, NScoords, partcoords, partvels, rho, u, num_parts, time, umass, udist, utime = readFILES(QGfilenames[i], 
#                                                                                                     "{}disk001".format(inputdirectory),
#                                                                                                     "{}IPDISK{}".format(inputdirectory, QGfilenames[i][-5:-2]))
#            opticalDepthPython((partcoords - Becoords) * (udist*0.01), partvels * (udist * 0.01)/utime, u * (udist * 0.01)**2/utime**2, num_parts, 
#                                 1e-15 * 0.001 * umass / 2e30, (7 * 6.96e10), 2.5, 90, 90)
#            print(time * utime / (24 * 3600))
#            name = "{}_theta{}_phi{}_time{:0.1f}".format(baseName, theta, phi, time * utime / (24 * 3600))
#            os.system("mv optDepthsFromPy.txt /data/rob1g10/SPH/DATA/DISKS/{}/optDepthData{}.txt".format(baseName, name))
#    else:
#        print("optDepthCode must be 'C' or 'Python'.")    
#        
#    return

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
#        diskname = 'INPUT/disk001'
#        IPname = 'INPUT/IPDISK411'
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
        print filename[-26:], name, "done"
    os.system("rm lsQG.txt lsIP.txt")
    
    return ionFrac
    

def plotRotatingOptDepthSpectraFigs(CcodeName, filenames, pngdirectory, baseName, thetaRange, phiRange):
    """
    Produces a scatter plots of varying viewing angles of a system
    in 2D of the particles and stars.

    GEOMETRY
    --------

    z axis is height of disk, y axis is width and x axis is depth.

    INPUT VARIABLES
    --------------- c 
    
    filename (list of 3 strings) - name of the filename with the particle 
                                   data on it and name of disk file for 
                                   units and also IP disk file name.

    pngdirectory (string) - directory where png files are found.

    baseName (string) - base of the name of the file that the plot will be saved as.

    OUTPUT VARIABLES
    ----------------

    Various plots of the particles.

    """
    
#    for theta in range(0, 181, 15):
#        for phi in range(0, 181, 15):
    for theta in range(thetaRange[0], thetaRange[1], thetaRange[2]):
        for phi in range(phiRange[0], phiRange[1], phiRange[2]):
                Becoords, NScoords, partcoords, partvels, rho, u, num_parts, time, umass, udist, utime = readFILES(filenames[0], filenames[1], filenames[2])
                nameZoom = "{}{}_theta{}_phi{}_time{:0.1f}_zoomY.png".format(pngdirectory, baseName, theta, phi, time * utime / (24 * 3600))
                name = "{}{}_theta{}_phi{}_time{:0.1f}_zoomN.png".format(pngdirectory, baseName, theta, phi, time * utime / (24 * 3600))
                oD = readOptDepthData4py("/data/rob1g10/SPH/DATA/DISKS/{}/optDepthData{}{}_theta{}_phi{}_time{:0.1f}.txt".format(baseName, CcodeName, baseName, theta, phi, time * utime / (24 * 3600)))
                T = np.exp(-oD) / 1e-3
                print(theta, phi)
                print(len(T))
                plotAndSaveSpectralLine(partvels, T, [phi, theta], 200, udist, 
                        utime, "uniform", nameZoom, save='y')
                plotAndSaveSpectralLine(partvels, T, [phi, theta], 200, udist, 
                        utime, "uniform", name, save='y', xlims="y")
    return


def plotRotatingOptDepthFigs(CcodeName, filenames, pngdirectory, baseName, thetaRange, phiRange, obsAngles, zoom):
    """
    Produces a scatter plots of varying viewing angles of a system
    in 2D of the particles and stars.

    GEOMETRY
    --------

    z axis is height of disk, y axis is width and x axis is depth.

    INPUT VARIABLES
    ---------------
    
    filename (list of 3 strings) - name of the filename with the particle 
                                   data on it and name of disk file for 
                                   units and also IP disk file name.

    pngdirectory (string) - directory where png files are found.

    baseName (string) - base of the name of the file that the plot will be saved as.

    OUTPUT VARIABLES
    ----------------

    Various plots of the particles.

    """
    
    for theta in range(thetaRange[0], thetaRange[1], thetaRange[2]):
        for phi in range(phiRange[0], phiRange[1], phiRange[2]):
            for zoom in [zoom]:#, 10, 5]:
                Becoords, NScoords, partcoords, partvels, rho, u, num_parts, time, umass, udist, utime = readFILES(filenames[0], filenames[1], filenames[2])
                name = "{}{}_theta{}_phi{}_zoom{}_obsAngle{}_time{:0.1f}.png".format(pngdirectory, baseName, theta, phi, zoom, obsAngles[0], time * utime / (24 * 3600))
                OD = readOptDepthData4py("/data/rob1g10/SPH/DATA/DISKS/{}/optDepthData{}{}_theta{}_phi{}_time{:0.1f}.txt".format(baseName, CcodeName, baseName, obsAngles[0], obsAngles[1], time * utime / (24 * 3600)))
                plotAndSaveFig([filenames[0], filenames[1]] , time * utime / (24 * 3600), name, 'OptDepth', orientation=[0 + phi, 0 + theta], 
                    show='n', save='y', xylims=[-zoom, zoom, -zoom, zoom], T=np.exp(-OD))
                print(theta, phi)
    
    Tsum = 0
    xcoords = partcoords[:,0] * udist / (7 * 6.96e10)
    ycoords = partcoords[:,1] * udist / (7 * 6.96e10)
    R = np.sqrt(xcoords**2 + ycoords**2)
    for i in range(len(OD)):
        if(R[i] < 10):
            Tsum += np.exp(-OD[i])

    return Tsum, partcoords

def plotRotatingOptDepthFigsDropbox(ODName, filenames, pngdirectory, baseName, plotType, thetaRange, phiRange, obsAngles, zoom):
    """
    Produces a scatter plots of varying viewing angles of a system
    in 2D of the particles and stars.

    GEOMETRY
    --------

    z axis is height of disk, y axis is width and x axis is depth.

    INPUT VARIABLES
    ---------------
    
    filename (list of 3 strings) - name of the filename with the particle 
                                   data on it and name of disk file for 
                                   units and also IP disk file name.

    pngdirectory (string) - directory where png files are found.

    baseName (string) - base of the name of the file that the plot will be saved as.

    OUTPUT VARIABLES
    ----------------

    Various plots of the particles.

    """
    
    for theta in range(thetaRange[0], thetaRange[1], thetaRange[2]):
        for phi in range(phiRange[0], phiRange[1], phiRange[2]):
            for zoom in [zoom]:#, 10, 5]:
                Becoords, NScoords, partcoords, partvels, rho, u, num_parts, time, umass, udist, utime = readFILES(filenames[0], filenames[1], filenames[2])
                name = "{}{}_theta{}_phi{}_zoom{}_obsAngle{}_time{:0.1f}.png".format(pngdirectory, baseName, theta, phi, zoom, obsAngles[0], time * utime / (24 * 3600))
                OD = readOptDepthData4py("{}".format(ODName)) #* 1e2
                print "Max OD:", max(OD)
                plotAndSaveFig([filenames[0], filenames[1]] , time * utime / (24 * 3600), name, plotType, orientation=[0 + phi, 0 + theta], 
                    show='n', save='y', xylims=[-zoom, zoom, -zoom, zoom], T=np.exp(-OD))
                print(theta, phi)
    
    Tsum = 0
    Tsum = np.sum(np.exp(-OD))
    ODplot = np.array([])
    Zplot = np.array([])
    xcoords = partcoords[:,0] * udist / (7 * 6.96e10)
    ycoords = partcoords[:,1] * udist / (7 * 6.96e10)
    R = np.sqrt(xcoords**2 + ycoords**2)
    #for i in range(len(OD)):
    #    if(OD[i] > 80):
    #        #Tsum += np.exp(-OD[i])
    #        ODplot = np.append(ODplot, OD[i])
    #        Zplot = np.append(Zplot, partcoords[i,2] * udist / (7 * 6.96e10))

    
    #print(len(OD), len(partcoords[:,2]))
    #plt.scatter(ODplot, Zplot, s=1)

    return Tsum, partcoords


def plotRotatingOptDepthHistsDropbox(ODName, filenames, thetaRange, colour, show="y"):
    """
    Produces a scatter plots of varying viewing angles of a system
    in 2D of the particles and stars.

    GEOMETRY
    --------

    z axis is height of disk, y axis is width and x axis is depth.

    INPUT VARIABLES
    ---------------
    
    filename (list of 3 strings) - name of the filename with the particle 
                                   data on it and name of disk file for 
                                   units and also IP disk file name.

    pngdirectory (string) - directory where png files are found.

    baseName (string) - base of the name of the file that the plot will be saved as.

    OUTPUT VARIABLES
    ----------------

    Various plots of the particles.

    """
    colourCount = 0
    for theta in range(thetaRange[0], thetaRange[1], thetaRange[2]):
        #for phi in range(phiRange[0], phiRange[1], phiRange[2]):
            #for zoom in [zoom]:#, 10, 5]:
        Becoords, NScoords, partcoords, partvels, rho, u, num_parts, time, umass, udist, utime = readFILES(filenames[0], filenames[1], filenames[2])
        OD = 0.48 * readOptDepthData4py("{}{}.txt".format(ODName, theta)) #* 1e2
        count = 0
        for i in range(len(OD)):
            if(OD[i] < 1e-9):
                count += 1
        print count
        for i in range(len(OD)):
            if(OD[i] < 1e-9):
                OD[i] = 1e-9
        counts, bins = np.histogram(OD, bins=np.logspace(-9, 2, 20))
        #plt.close()
        xplot = np.array([bins[0]])
        yplot = np.array([counts[0]])
        for i in range(1, len(bins)-1):
            xplot = np.append(xplot, bins[i])
            xplot = np.append(xplot, bins[i])
            yplot = np.append(yplot, counts[i-1])
            yplot = np.append(yplot, counts[i])
        xplot = np.append(xplot, bins[-1])
        yplot = np.append(yplot, counts[-1])
        #plt.plot(bins[:-1], counts, label="{}".format(theta), color=colour[colourCount]) 
        plt.plot(xplot, yplot, label="{}".format(theta), color=colour[colourCount])       
        colourCount += 1
        print(theta, phi)

    if(show == "y"):
        plt.xlabel("Optical Depth")
        plt.ylabel("Number of particles")
        plt.legend(loc="best")
        plt.ylim(1, 1.3e4)
        #plt.xlim(0, 5)
        #plt.yscale('log', nonposy='clip')
        plt.xscale('log')#, nonposy='clip')
        plt.show()

    return

def opticalDepthPython(filenames, outfilename, BeRadius, powerLawIndex, theta, phi, intSteps=10):
                             
    Becoords, NScoords, partcoords, partvels, rho, u, numparts, time, umass, udist, utime = readFILES(filenames[0], filenames[1], filenames[2])

    partcoords *= udist
    partvels *= udist/utime 
    u *= udist**2/utime**2 
    rho *= umass / udist**3
#    partmass = 1e-15 * umass / 2e30, 
#    BeRadius = (7 * 6.96e10)
#    BeMass = 18 * 2e33
#    powerLawIndex = 2.5, 
                           
    opticalDepth = np.zeros(numparts)
    
    CUTOFF = 15
    
    percentCount = 0
    
    m_e = 9.1e-28 #"mass of electron"
    e = 4.8e-10 #"electron charge"
    lambda0 = 656.28e-7 #"rest wavelength"
    c = 3e10 #"speed of light"
    f = 0.6407 #"absorption oscillator strength of the transition. for n=3 to n=2 Balmer line (emitting at 656nm), this is 3.03x10^-19"
    n_rs = 6.022e18 #"population density of level 2 electrons at stellar radius"
    kb = 1.38e-23 #"Boltzmann constant"
    idealGasConstant = 8.314 #"ideal gas constant"
    sigma = 6.3e-18 # absorption cross section at the Lyman limit for hydrogen
    nuFrac = (91.2 / 656.28)**3.5 # fraction for Lyman limit (91.2nm)
    " calc coefficient "
    coefficient = (np.sqrt(np.pi) * e*e * lambda0 * f) / (m_e * c * np.sqrt(2));

    xcoords = partcoords[:,0]
    ycoords = partcoords[:,1]
    zcoords = partcoords[:,2]
    vx = partvels[:,0]
    vy = partvels[:,1]
    vz = partvels[:,2]
    
    theta = np.deg2rad(theta)
    phi = np.deg2rad(phi)
 
    " Rotate coordinates and velocities "

    xcoords = xcoords * np.cos(phi) - ycoords * np.sin(phi)
    ycoords = (xcoords * np.sin(phi) + ycoords * np.cos(phi)) * np.cos(theta) - zcoords * np.sin(theta)
    zcoords = ycoords * np.sin(theta) + zcoords * np.cos(theta)
    vx = vx * np.cos(phi) - vy * np.sin(phi)
    vy = (vx * np.sin(phi) + vy * np.cos(phi)) * np.cos(theta) - vz * np.sin(theta)
    vz = vy * np.sin(theta) + vz * np.cos(theta)     
        
    "calculate radius for all particles"
    R = np.sqrt(xcoords**2 + ycoords**2 + zcoords**2)

    "calculate speed"
    v = np.sqrt(vx**2 + vy**2 + vz**2)

    "calculate scale Height"
#    H = (R / BeRadius) ** 1.5
    
    """
    calculate thermal velocity for all particles
    temperature -> U = 0.5 * (3 * kb * T) _-> 2 * U = 3 * kb * T
    thermal velocity = sqrt (3 * kb * T / m) = sqrt (2 * U / m)
    """
#    v_th = np.sqrt(2 * u / partmass)
    
    """
    calculate ionised fraction of hydrogen for all particles
    """
    ionCoeff = (2 * np.pi * (9.1e-31 * 1.38e-23 * 1e4) / 6.626e-34**2)
    ionExp = np.exp( (-13.6 * 1.6e-19) / (1.38e-23 * 1e4) )
    ionNumDens = 6.022 * pow(10, 23) * rho * 1e6
    ionA = 1
    ionB = (1 / ionNumDens) * ionCoeff**1.5 * ionExp
    ionC = -ionB
    ionFrac = (-ionB + np.sqrt(ionB*ionB - 4*ionA*ionC)) / (2*ionA)

    "calculate line of sight velocity and distance for all particles"
    "inclinationAngle = 0 corresponds to birdseye view of disk"
    v_D = vz;
    z_D = zcoords;
    maxz_D = max(z_D)
#    behindStarCount = 0
    for i in range(numparts):
        intSum = 0
        stepSize = (maxz_D - z_D[i]) / intSteps
        
        for j in range(intSteps):
            " IF PARTICLE IS BEHIND BE STAR, GIVE IT INFINITE OPTICAL DEPTH "
#            if(np.sqrt(xcoords[i]**2 + ycoords[i]**2) < BeRadius and zcoords[i] < 0
#                or intSum > CUTOFF):
#                intSum = CUTOFF
##                behindStarCount += 1
#                break
            a = z_D[i] + (j*stepSize) # integration lower limit
            b = a + stepSize # integration upper limit
            
            "find nearest particle to end a and get values"
            partA = np.argmin( np.sqrt( (xcoords - xcoords[i])**2 + (ycoords - ycoords[i])**2 + (zcoords - a)**2))
#            print i, partA
#            if(i%10000 == 0):
#                print "{}%".format(i/ 1000)
            "calc modified population density"
            if(1-ionFrac[partA] > 0):
                n = (1-ionFrac[partA]) * 6.022e22 * rho[partA]
            else:
                n = 0
            "add this to the sum"
#            intSum += (b-a) * coefficient * n * np.exp(-0.5 * (v[partA] - v_D[partA]) ** 2  / (v_th[partA] ** 2)) / v_th[partA]  
            intSum += (b-a) * n * sigma * nuFrac
    
        opticalDepth[i] = intSum
        percentage = 100. * i / numparts
        if(percentage > 10 and percentCount == 0):
            print "10% done"
            percentCount +=1
        if(percentage > 20 and percentCount == 1):
            print "20% done"
            percentCount +=1
        if(percentage > 30 and percentCount == 2):
            print "30% done"
            percentCount +=1
        if(percentage > 40 and percentCount == 3):
            print "40% done"
            percentCount +=1
        if(percentage > 50 and percentCount == 4):
            print "50% done"
            percentCount +=1
        if(percentage > 60 and percentCount == 5):
            print "60% done"
            percentCount +=1
        if(percentage > 70 and percentCount == 6):
            print "70% done"
            percentCount +=1
        if(percentage > 80 and percentCount == 7):
            print "80% done"
            percentCount +=1
        if(percentage > 90 and percentCount == 8):
            print "90% done"
            percentCount +=1
#    print behindStarCount
    f = open(outfilename, "w")
    for i in range(numparts):
        f.write("{}\n".format(opticalDepth[i]))
    f.close()

    print "{} done".format(outfilename)

    return opticalDepth

"""----------------------------------------------------------------------------

                ANALYTIC TESTS

----------------------------------------------------------------------------"""

def analyticData(xsteps, ysteps, zsteps):
    """
    Creates analytic data for a set number of particles in an isolated Be
    disc.

    GEOMETRY
    --------

    z axis is height of disk, y axis is width and x axis is depth.

    INPUT VARIABLES
    ---------------
    
    xsteps - number of particle steps in the x direction

    ysteps - number of particle steps in the x direction

    zsteps - number of particle steps in the z direction

    OUTPUT VARIABLES
    ----------------

    The coordinates, velocities, densities and internal energies
    of the particles.

    """

    num_parts = xsteps * ysteps * zsteps

    Becoords = np.array([0,0,0])
    NScoords = np.array([0,0,0])
    partcoords = np.zeros((num_parts, 3))
    partvels = np.zeros((num_parts, 3))
    rho = np.zeros(num_parts)
    u = np.zeros(num_parts)
    #for i in range(xsteps):
    #    for j in range(ysteps):
    #        for k in range(zsteps):
    #            partcoords[]

    umass = udist = utime = 1
    return Becoords, NScoords, partcoords, partvels, rho, u, num_parts, time, umass, udist, utime


def analyticEW(RHO, discRadius, discFlare, xSteps, ySteps, zSteps, thetaRange, phi, CUTOFF):

    m_e = 9.1e-28 #"mass of electron"
    e = 4.8e-10 #"electron charge"
    lambda0 = 656.28e-7 #"rest wavelength"
    c = 3e10 #"speed of light"
    f = 0.6407 #"absorption oscillator strength of the transition. for n=3 to n=2 Balmer line (emitting at 656nm), this is 3.03x10^-19"
    n_rs = 6.022e18 #"population density of level 2 electrons at stellar radius"
    kb = 1.38e-16 #"Boltzmann constant"
    h = 6.626e-27 # Planck constant
    idealGasConstant = 8.314 #"ideal gas constant"
    sigma = 6.3e-18 # absorption cross section at the Lyman limit for hydrogen
    nuFrac = (91.2 / 656.28)**3.5 # fraction for Lyman limit (91.2nm)
    rho0 = 2e-11 # base gas density
    alpha = 3.8 # power law exponent
    G = 6.67e-8 # gravitational constant
    M_stellar = 18 * 1.989e33 # stellar mass
    amu = 1.66e-24 # atomic mass unit

    #xcoords = xcoords * np.cos(phi) - ycoords * np.sin(phi)
    #ycoords = (xcoords * np.sin(phi) + ycoords * np.cos(phi)) * np.cos(theta) - zcoords * np.sin(theta)
    #zcoords = ycoords * np.sin(theta) + zcoords * np.cos(theta)

    filename = 'DATA/QGDISK41101'
    diskname = 'INPUT/disk001'
    IPname = 'INPUT/IPDISK411'

    Becoords, NScoords, partcoords, partvels, rhoparts, u, num_parts, time, umass, udist, utime = readFILES(filename, diskname, IPname)

    partcoords *= udist / (7 * 6.96e10)
    Rparts = np.sqrt(partcoords[:,0]**2 + partcoords[:,1]**2)
    Zparts = partcoords[:,2]
    rhoparts *= umass / udist**3

    #print min(Rparts), max(Rparts)
    #print min(Zparts), max(Zparts)

    # ionFrac Calculation
    #ionCoeff = (2 * pi * m_e * kb * T) / (h*h);
    #ionExp = exp( (-13.6 * 1.6e-12)) / (kb * T))
    #ionNumDens = 6.022e23 * rhos;
    #ionA = 1;
    #ionB = (1 / ionNumDens) * ionCoeff**1.5 * ionExp;
    #ionC = -ionB;
    #ionFrac = (-ionB + sqrt(ionB*ionB - 4*ionA*ionC)) / (2*ionA);

    rhoInDiscs = np.array([])
    EWs = np.array([])

    Xset = np.array([])
    Yset = np.array([])
    Zset = np.array([])

    for thetaDegs in range(thetaRange[0], thetaRange[1], thetaRange[2]):

        XsetTemp = np.array(Xset)
        YsetTemp = np.array(Yset)
        ZsetTemp = np.array(Zset)

        theta = np.deg2rad(thetaDegs)
        discFlare = np.deg2rad(discFlare)

        Rmin = 1
        zmin = 0.03
        zmax = discRadius #* np.tan(discFlare)
        #print discFlare

        EWsum = 0
        xmax = discRadius 
        ymax = discRadius

        xmaxrot = xmax
        ymaxrot = ymax #np.abs(ymax * np.cos(theta) - zmax * np.sin(theta))
        zmaxrot = zmax #np.abs(ymax * np.sin(theta) + zmax * np.cos(theta))

        #print xmaxrot, ymaxrot, zmaxrot

        "lines of line of sight"
        xLOS = np.linspace(xmaxrot, -xmaxrot, xSteps)
        yLOS = np.linspace(ymaxrot, -xmaxrot, ySteps)
        zLOS = np.linspace(zmaxrot, -xmaxrot, zSteps)

        Rs = np.array([])
        Zs = np.array([])
        ns = np.array([])
        Rplot = np.array([])
        Zplot = np.array([])
        nplot = np.array([])
        tauSums = np.array([])
        ionFracs = np.array([])
        ionFracs2 = np.array([])

        pointsInDisc = 0
        rhoInDisc = 0
        nInDisc = 0

        for i in xLOS:
            for j in yLOS:
                tauSum = 0
                #zold = -zmaxrot
                kold = zLOS[0]
                for k in zLOS:
                    # calculate radius and z of position
                    R = np.sqrt(i**2 + (j*np.cos(-theta) - k*np.sin(-theta))**2)
                    Z = j*np.sin(-theta) + k*np.cos(-theta)

                    Zmax = 1.00070310e-05*R**4 - 1.06179392e-03*R**3 + 2.70965278e-02*R**2 - 2.42317601e-02*R + 1.34718187e-01
                    if(thetaDegs == 0 and R < discRadius and R > 1 and np.abs(Z) < Zmax):
                        Xset = np.append(Xset, i)
                        Yset = np.append(Yset, (j*np.cos(-theta) - k*np.sin(-theta)))
                        Zset = np.append(Zset, j*np.sin(-theta) + k*np.cos(-theta))
                    elif(R < discRadius and R > 1 and np.abs(Z) < Zmax):
                        try:
                            cGP = np.argmin((XsetTemp - i)**2 + (YsetTemp - (j*np.cos(-theta) - k*np.sin(-theta)))**2 + 
                                            (ZsetTemp - j*np.sin(-theta) + k*np.cos(-theta))**2) # closest grid point from previous calculations
                            R = np.sqrt(XsetTemp[cGP]**2 + YsetTemp[cGP]**2)
                            Z = ZsetTemp[cGP]

                            #delete used particle
                            XsetTemp = np.delete(XsetTemp, cGP)
                            YsetTemp = np.delete(YsetTemp, cGP)
                            ZsetTemp = np.delete(ZsetTemp, cGP)
                        except:
                            print i, j, k
                            

                    # FOR COMPARISON WITH REAL DISC

                    #closestPart = np.argmin(np.sqrt((R - Rparts)**2 + (Z - Zparts)**2))
                    #xrot = i
                    #yrot = j*np.cos(-theta) - k*np.sin(-theta)
                    #zrot = j*np.sin(-theta) + k*np.cos(-theta)
                    #closestPart = np.argmin(np.sqrt((xrot - partcoords[:,0])**2 + (yrot - partcoords[:,1])**2 + (zrot - partcoords[:,2])**2))
                    #rho = rhoparts[closestPart]

                    #Z = 1
                    Rs = np.append(Rs, R)
                    Zs = np.append(Zs, Z)
                    # calcualte H
                    #v_kep = np.sqrt(G * M_stellar / R)
                    # calculate Temperature at given radius
                    #T = twoPartBeTemp2([R], 1.6e4, 0.9e4, 1.6e4, 10, 15)[0]
                    T = 1.3e4
                    # speed of sound = sqrt ( gamma * P / rho )
                    c_s = np.sqrt(1.4 * kb * T / amu)
                    # scale factor
                    H = c_s * R**1.5 / np.sqrt( G * M_stellar )
                    # calculate density at given radius and z
                    #if(R < 3.): 
                    #    alpha = 0.7
                    #elif(R < 5): 
                    #    alpha = 2.5
                    #elif(R < 7.5):
                    #    alpha = 3.1
                    #elif(R < 10):
                    #    alpha = 3.3
                    #else:
                    #    alpha = 4.
                    #rho = rho0 * (R ** -alpha) * np.exp(-0.5 * Z**2 / R**1.5)
                    pfit = np.array([-0.09690832,  -0.45785672,  -0.53858783,  -2.14967644, -10.75693424])
                    if( R < 1):
                        xfit = np.log10(1.0)
                    else:
                        xfit = np.log10(R)
                    #if(R < 1.1):
                    #    rho = (8.58735734e-11*xfit - 8.00917462e-11) * np.exp(-0.5 * Z**2 / R**1.5)
                    #else:
                    rho = 10.**(pfit[0]*xfit**4. + pfit[1]*xfit**3. + pfit[2]*xfit**2. + pfit[3]*xfit + pfit[4]) * np.exp(-0.5 * Z**2.) #/ R**3)#H**2.)
                    if( np.isnan(rho) ):
                        print rho, (i, j, k)
                        rho = 0
                    #rho = 1e-15
                    # calculate ion fraction
                    ionCoeff = (2. * np.pi * m_e * kb * T) / (h*h)
                    ionExp = np.exp( (-13.6 * 1.6e-12) / (kb * T))
                    ionNumDens = 6.022e23 * rho
                    ionA = 1.
                    ionB = (1. / ionNumDens) * ionCoeff**1.5 * ionExp
                    ionC = -ionB
                    ionFrac = (-ionB + np.sqrt(ionB*ionB - 4.*ionA*ionC)) / (2.*ionA)

                    ionFracs = np.append(ionFracs, ionFrac)
                    ionFracs2 = np.append(ionFracs2, 1. - ionFrac)
                    # calculate modified population density
                    if(1.-ionFrac > 0):
                        n = (1.-ionFrac) * 6.022e23  * rho
                    else:
                        n = 0.
                    ns = np.append(ns, n)
                    # check point is within disc
    #                if(R > discRadius or np.abs(Z) > R * np.tan(discFlare)):
    #                    zold = k
                    Zmax = 1.00070310e-05*R**4 - 1.06179392e-03*R**3 + 2.70965278e-02*R**2 - 2.42317601e-02*R + 1.34718187e-01
                    if(R > discRadius or R < 1 or np.abs(Z) > Zmax):
                        kold = k
                        #zold = k
                    else:
                        # calculate  optical depth
                        tauSum += np.abs(k - kold) * (7. * 6.96e10) * n * sigma * nuFrac #(2.*zmaxrot / zSteps) 8e3
                        kold = k
                        #zold = k
                        # EWsum
                        EWsum += np.exp(-tauSum) * rho
                        pointsInDisc += 1
                        rhoInDisc += rho
                        nInDisc += n
                    #if(k == zLOS[19] or k == zLOS[20]):
                    #if(tauSum > CUTOFF):
                    #    break
                    #    if(np.abs(Z) < 1):
                    #        #if(j<0):
                        Rplot = np.append(Rplot, R)#-np.sqrt(i**2 + j**2))
                    #        #else:
                    #        Rplot = np.append(Rplot, j*np.cos(-theta) - k*np.sin(-theta))#np.sqrt(i**2 + j**2))
                    #        Zplot = np.append(Zplot, i)#k)
                    #        nplot = np.append(nplot, n)
                    #        tauSums = np.append(tauSums, tauSum)
            #print Zmax
            #try:
            #    print tauSum #max(ionFracs), max(ionFracs2) #rho#, np.exp(-0.5 * Z**2 / R **1.5), (R ** -alpha)
            #except:
            #    pass
        #plt.scatter(Rs, Zs, s=0.5)
        #R = np.log10(1)
        #print 10.**(pfit[0]*R**4. + pfit[1]*R**3. + pfit[2]*R**2. + pfit[3]*R + pfit[4])
        #tauSums = np.exp(-tauSums)
        #cm = plt.cm.get_cmap('rainbow')
        #sc = plt.scatter(Rplot, Zplot, c=tauSums, vmin=np.min(tauSums), vmax=np.max(tauSums), cmap=cm, edgecolor='')
        #cb = plt.colorbar(sc)
        #plt.xlabel("R / stellar radii")
        #plt.ylabel("z / stellar radii")
        #plt.ylim(-25, 25)
        #plt.xlim(-25, 25)
        #plt.scatter(np.log10(ionFracs), np.log10(ionFracs2))
        #plt.hist(Rplot, bins=25)
        print EWsum, pointsInDisc, rhoInDisc, thetaDegs
        EWs = np.append(EWs, EWsum)
        rhoInDiscs = np.append(rhoInDiscs, rhoInDisc)

    #EWs = np.append(EWs[1:5], np.append(EWs[0], EWs[5:]))
    I1 = np.array([ 91.92051094,  90.46173838,  88.29814729,  88.47741181,  89.37006442, 89.96901982,  89.9781312,   89.2090977,   87.90792666,  80.54662987])
    #EWs *= (np.max(I1) - np.min(I1)) / (np.max(EWs) - np.min(EWs))
    #EWs += np.max(I1) - np.max(EWs) 
    EWs *= (I1[5] - np.min(I1)) / (EWs[5] - np.min(EWs))
    EWs += I1[5] - EWs[5] 

    plt.plot(range(0,91,10), EWs, label="Analytic {}".format(gridSteps))
    #plt.plot(range(0,91,10), rhoInDiscs / max(rhoInDiscs))
    #plt.plot(range(0,91,10), I1)
    return EWsum, pointsInDisc, rhoInDisc, nInDisc\


def analyticEW2(RHO, discRadius, discFlare, xSteps, ySteps, zSteps, intSteps, thetaRange, phi, CUTOFF):

    m_e = 9.1e-28 #"mass of electron"
    e = 4.8e-10 #"electron charge"
    lambda0 = 656.28e-7 #"rest wavelength"
    c = 3e10 #"speed of light"
    f = 0.6407 #"absorption oscillator strength of the transition. for n=3 to n=2 Balmer line (emitting at 656nm), this is 3.03x10^-19"
    n_rs = 6.022e18 #"population density of level 2 electrons at stellar radius"
    kb = 1.38e-16 #"Boltzmann constant"
    h = 6.626e-27 # Planck constant
    idealGasConstant = 8.314 #"ideal gas constant"
    sigma = 6.3e-18 # absorption cross section at the Lyman limit for hydrogen
    nuFrac = (91.2 / 656.28)**3.5 # fraction for Lyman limit (91.2nm)
    rho0 = 2e-11 # base gas density
    alpha = 3.8 # power law exponent
    G = 6.67e-8 # gravitational constant
    M_stellar = 18 * 1.989e33 # stellar mass
    amu = 1.66e-24 # atomic mass unit

    #xcoords = xcoords * np.cos(phi) - ycoords * np.sin(phi)
    #ycoords = (xcoords * np.sin(phi) + ycoords * np.cos(phi)) * np.cos(theta) - zcoords * np.sin(theta)
    #zcoords = ycoords * np.sin(theta) + zcoords * np.cos(theta)

    filename = 'DATA/QGDISK41101'
    diskname = 'INPUT/disk001'
    IPname = 'INPUT/IPDISK411'

    Becoords, NScoords, partcoords, partvels, rhoparts, u, num_parts, time, umass, udist, utime = readFILES(filename, diskname, IPname)

    partcoords *= udist / (7 * 6.96e10)
    Rparts = np.sqrt(partcoords[:,0]**2 + partcoords[:,1]**2)
    Zparts = partcoords[:,2]
    rhoparts *= umass / udist**3

    #print min(Rparts), max(Rparts)
    #print min(Zparts), max(Zparts)

    # ionFrac Calculation
    #ionCoeff = (2 * pi * m_e * kb * T) / (h*h);
    #ionExp = exp( (-13.6 * 1.6e-12)) / (kb * T))
    #ionNumDens = 6.022e23 * rhos;
    #ionA = 1;
    #ionB = (1 / ionNumDens) * ionCoeff**1.5 * ionExp;
    #ionC = -ionB;
    #ionFrac = (-ionB + sqrt(ionB*ionB - 4*ionA*ionC)) / (2*ionA);

    rhoInDiscs = np.array([])
    EWs = np.array([])

    for thetaDegs in range(thetaRange[0], thetaRange[1], thetaRange[2]):

        theta = np.deg2rad(thetaDegs)
        discFlare = np.deg2rad(discFlare)

        Rmin = 1
        zmin = 0.03
        zmax = discRadius #* np.tan(discFlare)
        #print discFlare

        EWsum = 0
        xmax = discRadius 
        ymax = discRadius

        xmaxrot = xmax
        ymaxrot = ymax #np.abs(ymax * np.cos(theta) - zmax * np.sin(theta))
        zmaxrot = zmax #np.abs(ymax * np.sin(theta) + zmax * np.cos(theta))

        #print xmaxrot, ymaxrot, zmaxrot

        "lines of line of sight"
        xLOS = np.linspace(xmaxrot, -xmaxrot, xSteps)
        yLOS = np.linspace(ymaxrot, -xmaxrot, ySteps)
        zLOS = np.linspace(zmaxrot, -xmaxrot, zSteps)

        Rs = np.array([])
        Zs = np.array([])
        ns = np.array([])
        Rplot = np.array([])
        Zplot = np.array([])
        nplot = np.array([])
        tauSums = np.array([])
        ionFracs = np.array([])
        ionFracs2 = np.array([])

        pointsInDisc = 0
        rhoInDisc = 0
        nInDisc = 0

        for i in xLOS:
            for j in yLOS:
                for k in zLOS:
                    tauSum = 0
                    for l in range(0, intSteps):
                        # calculate radius and z of position
                        #R = np.sqrt(i**2 + (j*np.cos(-theta) - k*np.sin(-theta))**2)
                        #Z = j*np.sin(-theta) + k*np.cos(-theta)
                        x = i
                        y = 0*np.cos(-theta) - discRadius*np.sin(-theta)
                        z = 0*np.sin(-theta) + discRadius*np.cos(-theta)
                        ystepsize = (y - j) / intSteps
                        zstepsize = (z - k) / intSteps
                        R = np.sqrt(x**2 + (y - l*ystepsize)**2)
                        Z = z - l*zstepsize
                        #if(y < 1):
                        #    print R, y, ystepsize
                        # FOR COMPARISON WITH REAL DISC

                        #closestPart = np.argmin(np.sqrt((R - Rparts)**2 + (Z - Zparts)**2))
                        #xrot = i
                        #yrot = j*np.cos(-theta) - k*np.sin(-theta)
                        #zrot = j*np.sin(-theta) + k*np.cos(-theta)
                        #closestPart = np.argmin(np.sqrt((xrot - partcoords[:,0])**2 + (yrot - partcoords[:,1])**2 + (zrot - partcoords[:,2])**2))
                        #rho = rhoparts[closestPart]

                        #Z = 1
                        Rs = np.append(Rs, R)
                        Zs = np.append(Zs, Z)
                        # calcualte H
                        #v_kep = np.sqrt(G * M_stellar / R)
                        # calculate Temperature at given radius
                        #T = twoPartBeTemp2([R], 1.6e4, 0.9e4, 1.6e4, 10, 15)[0]
                        T = 1.3e4
                        # speed of sound = sqrt ( gamma * P / rho )
                        c_s = np.sqrt(1.4 * kb * T / amu)
                        # scale factor
                        H = c_s * (R*(7. * 6.96e10))**1.5 / np.sqrt( G * M_stellar )
                        # calculate density at given radius and z
                        #if(R < 3.): 
                        #    alpha = 0.7
                        #elif(R < 5): 
                        #    alpha = 2.5
                        #elif(R < 7.5):
                        #    alpha = 3.1
                        #elif(R < 10):
                        #    alpha = 3.3
                        #else:
                        #    alpha = 4.
                        #rho = rho0 * (R ** -alpha) * np.exp(-0.5 * Z**2 / R**1.5)
                        pfit = np.array([-0.09690832,  -0.45785672,  -0.53858783,  -2.14967644, -10.75693424])
                        if( R < 1):
                            xfit = np.log10(1.0)
                        else:
                            xfit = np.log10(R)
                        #if(R < 1.1):
                        #    rho = (8.58735734e-11*xfit - 8.00917462e-11) * np.exp(-0.5 * Z**2 / R**1.5)
                        #else:
                        rho = 10.**(pfit[0]*xfit**4. + pfit[1]*xfit**3. + pfit[2]*xfit**2. + pfit[3]*xfit + pfit[4]) * np.exp(-0.5*Z**2.)# / R**3.)
                        if( np.isnan(rho) ):
                            print rho, (i, j, k)
                            rho = 0
                        #rho = 1e-15
                        # calculate ion fraction
                        ionCoeff = (2. * np.pi * m_e * kb * T) / (h*h)
                        ionExp = np.exp( (-13.6 * 1.6e-12) / (kb * T))
                        ionNumDens = 6.022e23 * rho
                        ionA = 1.
                        ionB = (1. / ionNumDens) * ionCoeff**1.5 * ionExp
                        ionC = -ionB
                        ionFrac = (-ionB + np.sqrt(ionB*ionB - 4.*ionA*ionC)) / (2.*ionA)

                        ionFracs = np.append(ionFracs, ionFrac)
                        ionFracs2 = np.append(ionFracs2, 1. - ionFrac)
                        # calculate modified population density
                        if(1.-ionFrac > 0):
                            n = (1.-ionFrac) * 6.022e23  * rho
                        else:
                            n = 0.
                        ns = np.append(ns, n)
                        # check point is within disc
        #                if(R > discRadius or np.abs(Z) > R * np.tan(discFlare)):
        #                    zold = k
                        Zmax = 1.00070310e-05*R**4 - 1.06179392e-03*R**3 + 2.70965278e-02*R**2 - 2.42317601e-02*R + 1.34718187e-01
                        if(R < discRadius or R > 1 or np.abs(Z) < Zmax):
                            # calculate  optical depth
                            tauSum += np.sqrt(ystepsize**2 + zstepsize**2) * (7. * 6.96e10) * n * sigma * nuFrac #(2.*zmaxrot / zSteps) 8e3
                            if(l==intSteps-1):
                                rhoInDisc += rho
                    
                    # EWsum
                    EWsum += np.exp(-tauSum) * rho
                    #pointsInDisc += 1
                    #rhoInDisc += rho
                    #nInDisc += n

        print EWsum, rhoInDisc, thetaDegs #, pointsInDisc, rhoInDisc, thetaDegs
        EWs = np.append(EWs, EWsum)
        rhoInDiscs = np.append(rhoInDiscs, rhoInDisc)

    #EWs = np.append(EWs[1:5], np.append(EWs[0], EWs[5:]))
    I1 = np.array([ 91.92051094,  90.46173838,  88.29814729,  88.47741181,  89.37006442, 89.96901982,  89.9781312,   89.2090977,   87.90792666,  80.54662987])
    EWs *= (np.max(I1) - np.min(I1)) / (np.max(EWs) - np.min(EWs))
    EWs += np.max(I1) - np.max(EWs) 
    #EWs *= (I1[5] - np.min(I1)) / (EWs[5] - np.min(EWs))
    #EWs += I1[5] - EWs[5] 

    plt.plot(range(0,91,10), EWs, label="Analytic {}, {}".format(gridSteps, intSteps))
    #plt.plot(range(0,91,10), rhoInDiscs / max(rhoInDiscs))
    #plt.plot(range(0,91,10), I1)
    return EWsum, pointsInDisc, rhoInDisc, nInDisc

"""----------------------------------------------------------------------------

                RUNNING FUNCTIONS

----------------------------------------------------------------------------"""

#filename = 'QGDISK09101'
#filename = '/data/rob1g10/SPH/OriginalOkazaki/images/imagep3DSP/MyRead/infoPartsDISK001'
#filename = '/data/rob1g10/SPH/RoryTest2/sph3D/RBDISK1.31000E+02'
#filename = '/data/rob1g10/SPH/RoryTest2/images/imagep3DSP/QGDISK09001'
#filename = '/data/rob1g10/SPH/RoryTest/sph3D/infoPartsDISK  1'
#filename = '/data/rob1g10/SPH/DATA/NSorbit_DISKsize_rel/a1.5e11_isotherm/processedDISKs/QGDISK10925'
#filename = '/data/rob1g10/SPH/DATA/DISKS/justBe_halfInj2/processedDISKs/QGDISK50301'
#filename = '/data/rob1g10/SPH/DATA/DISKS/justBe_quarterInj/processedDISKs/QGDISK59301'
#filename = '/data/rob1g10/SPH/DATA/DISKS/isoBe_fromBeginningTenth/processedDISKs/QGDISK41101'
filename = 'DATA/QGDISK41101'
diskname = 'INPUT/disk001'
IPname = 'INPUT/IPDISK411'
#filename = '/data/rob1g10/SPH/DATA/DISKS/20d_0.0e_eqTest/processedDISKs/QGDISK38701'
#diskname = filename[:-26] + "disk001"
#IPname = filename[:-11] + "IPDISK" + filename[-5:-2]
#filename = '/data/rob1g10/SPH/DATA/DISKS/20d_0.4e/processedDISKs/QGDISK28301'
#filename = '/data/rob1g10/SPH/DATA/DISKS/20d_0.4e/processedDISKs/QGDISK26501'
#filename = '/data/rob1g10/SPH/DATA/DISKS/40d_0.4e/processedDISKs/QGDISK13401'
#filename = '/data/rob1g10/SPH/DATA/DISKS/injection_off_JJrun/processedDISKs/QGDISK13601'
#filename = '/data/rob1g10/SPH/DATA/DISKS/injection_off_a1.5e11/processedDISKs/QGDISK14701'
#filename = '/data/rob1g10/SPH/DATA/DISKS/injection_off_justBe/processedDISKs/QGDISK17001'
#filename = '/data/rob1g10/SPH/injFiddle/sph3D/processedDISKs/QGDISK13201'
#filename = '/data/rob1g10/SPH/DATA/DISKS/injFiddle_largebox/processedDISKs/QGDISK32701'
#filename = '/data/rob1g10/SPH/DATA/DISKS/inclined/processedDISKs/QGDISK04819'
#filename = '/data/rob1g10/SPH/DATA/DISKS/inclined/processedDISKs/QGDISK08322'
# filename = '/data/rob1g10/SPH/DATA/DISKS/sph3D_f90_radialTemp/processedDISKs/QGDISK04501'
# diskname = '/data/rob1g10/SPH/DATA/DISKS/sph3D_f90_radialTemp/disk001'
# IPname = '/data/rob1g10/SPH/DATA/DISKS/sph3D_f90_radialTemp/processedDISKs/IPDISK045'
#filename = raw_input('Type the name of the ASCII file to be read') 
""" Read File """

Becoords, NScoords, partcoords, partvels, rho, u, num_parts, time, umass, udist, utime = readFILES(filename, diskname, IPname)
#print max(rho) * umass / udist**3
#
#sortOrder = np.argsort(partcoords[:,2])[::-1]
#print sortOrder
#
#partcoords = partcoords[sortOrder]
#partvels = partvels[sortOrder]
#rho = rho[sortOrder]
theta = 0
phi = 0
#ionFrac = createOptDepthData4c("TESTNEW", 
#                             partcoords * udist, partvels * udist/utime, 
#                             u * udist**2/utime**2, 
#                             rho * umass / udist**3, num_parts, 
#                             1e-15 * umass / 2e30, 
#                             (7 * 6.96e10), 18 * 2e33, 2.5, 
#                             max(rho) * umass / udist**3, theta, phi)
#print ionFrac
#plt.scatter((np.sqrt(partcoords[:,0]**2 + partcoords[:,1]**2) *udist) / (7 * 6.96e10), ionFrac, s=0.5)
#plt.ylim(-1e-6, 1.1*max(ionFrac))
#plt.scatter((np.sqrt(partcoords[:,0]**2 + partcoords[:,1]**2) *udist) / (7 * 6.96e10), 
#            1.3e4 / (np.sqrt(partcoords[:,0]**2 + partcoords[:,1]**2) *udist / (7 * 6.96e10)), s=0.5)
#plt.scatter(ionFrac, 
#            1 - 1 / np.sqrt(np.sqrt(partcoords[:,0]**2 + partcoords[:,1]**2) *udist), s=0.5)
#plt.ylim(-1e-6, 1.1*max(ionFrac))

#for theta in range(0, 91, 10):
#    opticalDepthPython([filename, diskname, IPname], "OUTPUT/PythonOD{}.txt".format(theta), 
#                       (7 * 6.96e10), 2.5, theta, 0)

"Temp func-------------------------------------------------------------------"
def beTemp(R, rho0, T0, minT, minTx):
    """ function calculate temperature for particles in a Be star disk base on radius."""
    #minTx -= 5
    T = (5000 * (4- (4*(R) + 6) / (0.1*(R)**3 + 6)))
    minTx -= R[np.argmin(T)]
    T = (5000 * (4- (4*(R - minTx) + 6) / (0.1*(R - minTx)**3 + 6)))
    return T

"2 part Temp func-------------------------------------------------------------------"
def twoPartBeTemp(R, rho0, T0, minT, minTxstart, minTxend, levelledT):
    """ function calculate temperature for particles in a Be star disk base on radius."""
    T = np.zeros_like(R)
    for i in range(len(R)):
        midMin = 0.5*(minTxend - minTxstart) + minTxstart
        if(R[i] < minTxstart):
            stretch = (T0 - minT) / (midMin - 1)**2
            T[i] = stretch * (R[i] - midMin)**2 + minT
        else:
            stretch = levelledT
            setMin = minT - (levelledT - minT)*(1 / (1+np.exp(-midMin+3+(minTxend-midMin))))
            T[i] = (levelledT - minT)*(1 / (1+np.exp(-R[i]+3+(minTxend-midMin)))) + setMin
    return T

def twoPartBeTemp2(R, T0, minT, levelledT, minTx, levelledx):
    """ function calculate temperature for particles in a Be star disk base on radius."""
    T = np.zeros_like(R)
    for i in range(len(R)):
        xstretch = 1 / minTx
        if(R[i] <= minTx):
            ystretch = (T0 - minT) / ((minTx - 1)**4)
            T[i] = ystretch * (R[i] - minTx)**4 + minT
        elif(R[i] > minTx and R[i] < levelledx):
            ystretch = (levelledT - minT) / ((minTx - 1)**2)
            T[i] = ystretch * (R[i] - minTx)**2 + minT
        else:
            T[i] = levelledT
            #T[i] = - 1000 / (R[i] + 2) + levelledT
    return T

" calcs by hand -------------------------------------------------------------"

theta = 0
# discFlare = 11.3 or 10.2
#for discFlare in range(10, 31, 5):
gridSteps = 30
##for rho in [1e-11, 1e-12, 1e-13, 1e-14, 1e-15]:
count = 0
#for gridSteps in [10, 20, 30, 40]:#, 14, 18]:#, 4, 6, 8, 10, 12]:
#    #for intSteps in [6]:#, 12, 16]:
#    #    intSteps = [8, 10, 20, 30][count]
#    count += 1
#    EWs = np.array([])
#    pointsInDiscs = np.array([])
#    rhoInDiscs = np.array([])
#    nInDiscs = np.array([])
#
#    #EW, pointsInDisc, rhoInDisc, nInDisc = analyticEW2(rho, 25, 10, gridSteps, gridSteps, gridSteps, intSteps, [0,91,10], 0, 100)
#    EW, pointsInDisc, rhoInDisc, nInDisc = analyticEW(rho, 25, 10, gridSteps, gridSteps, gridSteps, [0,91,10], 0, 100)
#    EWs = np.append(EWs, EW)# / 0.75e-12)
#    pointsInDiscs = np.append(pointsInDiscs, pointsInDisc)
#    rhoInDiscs = np.append(rhoInDiscs, rhoInDisc)
#    nInDiscs = np.append(nInDiscs, nInDisc)
    #theta = np.deg2rad(theta)
    #print 25, 25 * np.cos(theta) - 5 * np.sin(theta), 25 * np.sin(theta) + 5 * np.cos(theta)
    #print EWs[-1], theta

#    EWs = np.array([])
#    count = 0
#    for theta in range(0, 91, 10):
#        EW, pointsInDisc, rhoInDisc, nInDisc = analyticEW(rhoInDiscs[count], 25, 10, gridSteps, gridSteps, gridSteps, theta, 0, 100)
#        EWs = np.append(EWs, EW)# / 0.75e-12)
#        count += 1
#        #theta = np.deg2rad(theta)
#        #print 25, 25 * np.cos(theta) - 5 * np.sin(theta), 25 * np.sin(theta) + 5 * np.cos(theta)
#        print EWs[-1], theta
    #EWs /= pointsInDiscs
    #I1 = np.array([42, 43, 42, 42.5, 44, 46, 48, 47.5, 39, 24])
#    I1 = np.array([ 91.92051094,  90.46173838,  88.29814729,  88.47741181,  89.37006442, 89.96901982,  89.9781312,   89.2090977,   87.90792666,  80.54662987])
#    EWs *= (np.max(I1) - np.min(I1)) / (np.max(EWs) - np.min(EWs))
#    EWs += np.max(I1) - np.max(EWs)
    #EWs = np.append(EWs[0], EWs[0:-1])
    #plt.plot(range(0,91,10), EWs, label="Analytic {}".format(gridSteps))#"{} degrees".format(discFlare))
    #plt.plot(range(0,91,10), rhoInDiscs / max(rhoInDiscs), label="Analytic {}".format(gridSteps))
#plt.plot(range(0,91,10), rhoInDiscs1 / max(rhoInDiscs1), label="Analytic {}".format(gridSteps))
    #plt.plot(range(0,91,10), nInDiscs, label="Analytic {}".format(gridSteps))
#plt.legend(loc="best")
#plt.show()

# iso Be SPH Data
#steps48 = np.array([104.637709337, 106.201670945, 105.830646857, 100.217708517, 91.5878249126, 
#                    94.2286300503, 90.6869841091, 97.1217019663, 78.4684559926, 18.5450552703])
#plt.plot(range(0,91,10), steps48 / 2.2, label="Analytic using SPH data")#"{} degrees".format(discFlare))
#steps48T1.6_0.1_1.6 = np.array([49.6996724838, 49.1277233829, 46.7037548846, 42.817798982, 37.8021711116, 
#                                51.9697044416, 38.3547555424, 27.0281987133, 23.5945739068, 18.5470200817])

#steps48T1.6_0.9_1.6_R3-5 = np.array([65.0456214129, 62.5667849493, 62.1817446167, 56.4766124872, 52.6127786043,
#                                     70.2251125592, 64.5652873367, 43.1413846517, 38.3687814433, 18.5471233802])
#steps48T1.6_0.5_1.6_R3-5 = np.array([50.1161524716, 49.708375863, 48.9171333504, 45.5383685236, 41.3570529421, 
#                                     57.7275059966, 39.8060439539, 29.6417327485, 25.6581145064, 18.5470874757])
#steps48T1.6_0.9_1.6_R5p5-10 = np.array([83.218465969, 78.2212803861, 79.8121253709, 70.8157549373, 63.1996732716,
#                                        78.3963415373, 90.5974237596, 62.2756345148, 29.050018719, 17.9596553761])
#steps48T1.6_0.5_1.6_R5p5-10 = np.array([28.6063071389, 23.017095145, 23.5984767018, 17.8509536534, 14.0815117135,
#                                        32.3122137374, 40.9880408101, 12.22649362, 7.07841050979, 6.67365914454])
#steps48T1.6_0.9_1.6_R10-15 = np.array([103.08268336, 105.513456134, 104.70295183, 99.2166816502, 89.033017091,
#                                       105.650610784, 124.582666814, 102.526615773, 64.0623354442, 18.0963422623])
#steps48T1.6_0.5_1.6_R10-15 = np.array([34.1422687268, 30.7896973528, 32.6243843594, 26.7087098191, 22.5990272123,
#                                       41.3477682995, 62.5140144185, 39.8454527107, 5.56122437066, 1.24185378406])

"constants"
theta = np.deg2rad(0)
T = 1.3e4
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
print min(rho * umass / udist**3)

R = np.sqrt(partcoords[:,0]**2 + partcoords[:,1]**2) * udist / (7 * 6.96e10)
#T = (5000 * (4- (4*(R-2.5) + 6) / ((R-2.5)**3 + 6)))
#T1 = twoPartBeTemp(R, 1e-11, 1.5e4, 1e4, 2, 3, 1.5e4)
#T2 = twoPartBeTemp(R, 1e-11, 1.6e4, 0.9e4, 2, 4, 1.8e4)
#T3 = twoPartBeTemp(R, 1e-11, 2e4, 0.8e4, 2, 8, 1.4e4)
#T4 = twoPartBeTemp(R, 1e-11, 1.5e4, 1e4, 5, 8, 1.5e4)
"Fits"
"----"
T1 = twoPartBeTemp2(R, 1.5e4, 1e4, 1.5e4, 3, 5)
T2 = twoPartBeTemp2(R, 1.6e4, 0.9e4, 1.8e4, 5.5, 10)
T3 = twoPartBeTemp2(R, 2e4, 0.8e4, 1.4e4, 6.5, 12)
T4 = twoPartBeTemp2(R, 2.3e4, 0.8e4, 1.4e4, 10, 20)
T5 = twoPartBeTemp2(R, 2.3e4, 0.3e4, 1.4e4, 10, 20)
#T1 = twoPartBeTemp2(R, 1.5e4, 1e4, 1.5e4, 3, 5)
#T2 = twoPartBeTemp2(R, 1.6e4, 0.9e4, 1.8e4, 5.5, 10)
#T3 = twoPartBeTemp2(R, 2e4, 0.8e4, 1.4e4, 6.5, 12)
#T4 = twoPartBeTemp2(R, 2.3e4, 0.8e4, 1.4e4, 10, 20)
#T5 = twoPartBeTemp2(R, 2.3e4, 0.3e4, 1.4e4, 10, 20)
#T1 = twoPartBeTemp2(R, 1.6e4, 0.9e4, 1.6e4, 3, 5)
#T2 = twoPartBeTemp2(R, 1.6e4, 0.7e4, 1.6e4, 3, 5)
#T3 = twoPartBeTemp2(R, 1.6e4, 0.5e4, 1.6e4, 3, 5)  
#T4 = twoPartBeTemp2(R, 1.6e4, 0.9e4, 1.6e4, 5.5, 10)
#T5 = twoPartBeTemp2(R, 1.6e4, 0.7e4, 1.6e4, 5.5, 10)
#T6 = twoPartBeTemp2(R, 1.6e4, 0.5e4, 1.6e4, 5.5, 10)
#T7 = twoPartBeTemp2(R, 1.6e4, 0.9e4, 1.6e4, 8, 15)
#T8 = twoPartBeTemp2(R, 1.6e4, 0.7e4, 1.6e4, 8, 15)
#T9 = twoPartBeTemp2(R, 1.6e4, 0.5e4, 1.6e4, 8, 15)
#print max(R), min(R)
#print max(T), min(T)

#plt.scatter(R, T1, s=0.5, color='blue', label='5e-12')
#plt.scatter(R, T2, s=0.5, color='green', label='1e-11')
#plt.scatter(R, T3, s=0.5, color='black', label='5e-11')
#plt.scatter(R, T4, s=0.5, color='red', label='1e-10')
#plt.scatter(R, T5, s=0.5, color='yellow', label='BREAK')
#plt.scatter(R, T1, s=0.5, color='blue')
#plt.scatter(R, T2, s=0.5, color='blue')
#plt.scatter(R, T3, s=0.5, color='blue')
#plt.scatter(R, T4, s=0.5, color='red')
#plt.scatter(R, T5, s=0.5, color='red')
#plt.scatter(R, T6, s=0.5, color='red')
#plt.scatter(R, T7, s=0.5, color='green')
#plt.scatter(R, T8, s=0.5, color='green')
#plt.scatter(R, T9, s=0.5, color='green')
##plt.plot(range(0, 10), 2e4*np.ones(10))
##plt.plot(np.range(0))
#plt.xlim(1, 20)
#plt.xlabel("R / stellar radii")
#plt.ylabel("Temperature / K")
#plt.legend(loc="best")
#plt.show()

#coefficient = (np.sqrt(np.pi) * e**2 * lambda0 * f) / (m_e * c * np.sqrt(2))

#azimuthalAngle = np.tan(partcoords[:,0]/partcoords[:,1])
#c_s = np.sqrt(2 * (u* udist**2 / utime**2) / 3)
#R = np.sqrt(partcoords[:,0]**2 + partcoords[:,1]**2) * udist
#v_kep = np.sqrt(G * M_stellar / R)
#H = c_s * R**1.5 / (G * M_stellar)

#v_th = np.sqrt(3 * kb * T / amu)
#nu_th = nu0 * (v_th / c)
#Q =  (3 / 4) * (v_kep / v_th) * (H / R)
#delta_nu = nu_th * np.sqrt(1 + Q**2 * np.sin(theta)**2 * np.sin(azimuthalAngle)**2 * np.tan(theta)**2 * np.cos(azimuthalAngle)**2)

#v_D = v_kep * np.sin(theta) * np.sin(azimuthalAngle) * (1 - (2/3) * (partcoords[:,2] / R) 
#                                                        * np.tan(theta) * np.cos(azimuthalAngle) )
#nu_D = nu0 * (1 - (v_D / c))

"number density of ionised hydrogen"
colours = ["blue", "green", "black", "red", "yellow", "purple"]
#Ts = [T1]
#for i in range(len(Ts)):
#    T = Ts[i]
#tau1s = np.array([])
#Tplot = np.array([])
#rhos = [1e-11, 1e-12, 1e-13, 1e-14, 1e-15, 1e-16, 1e-17]
#Rs =   [1, 1e1**(1./3), 1e2**(1./3), 1e3**(1./3), 1e4**(1./3), 1e5**(1./3), 1e6**(1./3)]
#plt.scatter(R, np.log10(rho* umass / udist**3), s=0.5)
#plt.ylim(0, max(np.log(rho* umass / udist**3)))
#for i in range(len(rhos)):
#    rho = rhos[i]
#    R = Rs[i]
#    T = np.linspace(1, 15000, 1000)
#    ionCoeff = (2 * np.pi * m_e * kb * T) / (h*h)
#    ionExp = np.exp( (-13.6 * 1.6e-12) / (kb * T))
#    ionNumDens = 6.022e23 * rho * umass / udist**3
#    ionA = 1.
#    ionB = (1. / ionNumDens) * ionCoeff**1.5 * ionExp
#    ionC = -ionB
#    ionFrac = (-ionB + np.sqrt(ionB*ionB - 4*ionA*ionC)) / (2*ionA)
#    n = (1-ionFrac) * 6.022 * pow(10, 23) * rho #* umass / udist**3
#    #    #if(rho[i] < 1e-13):
#    #    #    print i, 1-ionFrac, - ionB + np.sqrt(ionB*ionB - 4*ionA*ionC), rho[i]
#    #    #print max(1-ionFrac)
#    #plt.scatter(R, 1-ionFrac, s=0.5, color=colours[i])
#    #plt.plot(T, n*(R * 7 * 6.96e10)*6.3e-18*(91.2/656.28)**3.5, label="{}".format(rho))
#    #plt.ylim(0, 1)#)max(1-ionFrac)*1.1)
#    tau = n*(R * 7 * 6.96e10)*6.3e-18*(91.2/656.28)**3.5
#    arg = np.argmin(np.abs(tau - 1))
#    tau1s = np.append(tau1s, tau[arg])
#    Tplot = np.append(Tplot, T[arg])
#plt.plot(np.log10(rhos), Tplot)
##plt.plot(Rs, Tplot)
##plt.ylim(4000, 9000)
#plt.xlim(-11, -15)
##plt.legend(loc="best")
#plt.show()
cm = plt.cm.get_cmap('rainbow')
#T = T9
#T = 1.3e4 * np.ones(len(rho))
#ionCoeff = (2 * np.pi * m_e * kb * T) / (h*h)
#ionExp = np.exp( (-13.6 * 1.6e-12) / (kb * T))
#ionNumDens = 6.022e23 * rho * umass / udist**3
#ionA = 1.
#ionB = (1. / ionNumDens) * ionCoeff**1.5 * ionExp
#ionC = -ionB
#ionFrac = (-ionB + np.sqrt(ionB*ionB - 4*ionA*ionC)) / (2*ionA)
#n = (1-ionFrac) * 6.022 * pow(10, 23) * rho #* umass / udist**3
#sc = plt.scatter(R, partcoords[:,2], c=n, s = 1, vmin=np.min(n), vmax=np.max(n), cmap=cm, edgecolor='')
#cb = plt.colorbar(sc)
#plt.xlim(1, 3)
#plt.ylim(-0.02, 0.02)
#plt.show()

#Ts = [5e3, 6e3, 7e3, 8e3, 9e3, 1e4]
#rhos = [1e-10, 1e-12, 1e-14, 1e-16, 1e-18]
#for i in range(len(rhos)):
#    rho = rhos[i]
#    #T = Ts[i]
#    T = np.linspace(2e3, 9e3, 100)
#    ionCoeff = (2 * np.pi * m_e * kb * T) / (h*h)
#    ionExp = np.exp( (-13.6 * 1.6e-12) / (kb * T))
#    ionNumDens = 6.022e23 * rho * umass / udist**3
#    ionA = 1.
#    ionB = (1. / ionNumDens) * ionCoeff**1.5 * ionExp
#    ionC = -ionB
#    ionFrac = (-ionB + np.sqrt(ionB*ionB - 4*ionA*ionC)) / (2*ionA)
#    #n = (1-ionFrac) * 6.022 * pow(10, 23) * rho * umass / udist**3
#    #tau = n*(R * 7 * 6.96e10)*6.3e-18*(91.2/656.28)**3.5
#    plt.plot(T, 1-ionFrac, label="{}".format(rho))#, color=colours[i])
##plt.xlim(0, 10)
##plt.ylim(0, 5)
#plt.xlabel("Temperature / K")
#plt.ylabel("Fraction of neutral hydrogen")
#plt.legend(loc="best")
#plt.show()


"ans"
#ans = (7 * 6.96e10) * (1/np.cos(theta)) * coefficient * n * np.exp(-0.5 * (nu0 - nu_D) * (nu0 - nu_D) / (nu_th * nu_th)) / nu_th
#tau0 = (7 * 6.96e10) * ((np.pi * e**2 * lambda0 * f) / (m_e * c * np.sqrt(2))) * n * H / v_th
#print max(tau0), min(tau0)

#coeff = (2 * np.pi * 3.1e-31 * 1.38e-23 * 1.3e4) / (6.626e-34**2)
#expTerm = np.exp(-(13.6 * 1.6e-19) / (1.38e-23 * 1.3e4))
#n = 6.023e18
#a = 1
#b = (1/n) * coeff**1.5 * expTerm
#c = -(1/n) * coeff**1.5 * expTerm
#print (-b + np.sqrt(b**2 - 4*a*c)) / 2*a

#
#optDepth = calcOpticalDepthOpt2(partcoords * (udist*0.01), rho * umass/((0.01*udist)**3), 
#                                0.01 * max(partcoords[:,0])  * (udist*0.01), 10)
#print optDepth
#createOptDepthData4c(partcoords * (udist*0.01), partvels * (udist * 0.01)/utime, u * (udist * 0.01)**2/utime**2, num_parts, 
#                                     1e-15 * 0.001 * umass / 2e30, (7 * 6.96e10), 2.5, 45, 0)

#calcRotatingOptDepths([filename, diskname, IPname], "20d_0.4e", [0,181,15], [0,181,15])
#calcOptDepthsOverTime('/data/rob1g10/SPH/DATA/Animations/JournalClub/processedDISKs/', "20d_0.4e", 90, 90)

#opticalDepthPython((partcoords - Becoords) * (udist*0.01), partvels * (udist * 0.01)/utime, u * (udist * 0.01)**2/utime**2, num_parts, 
#                                 1e-15 * 0.001 * umass / 2e30, (7 * 6.96e10), 2.5, 90, 90)

#plotRotatingOptDepthFigs([filename, diskname, IPname], "../Plots to Show/All Round Image/20d_0.4e_Spiral/", "20d_0.4e", [0,181,15], [0,181,15])
#oD = readOptDepthData4py()
#oD = readOptDepthData4py()
#for i in range(len(R)):
#    if(partcoords[i, 1] < 0.1 and partcoords[i, 1] > -0.1 and partcoords[i, 0] > 0.):
#        plt.scatter(R[i] * udist / (7 * 6.96e10), oD[i], s=0.25)
##plt.ylim(-3e-18, 3e-17)
#plt.xlim(0, 6)
#plt.ylabel("Optical depth")
#plt.xlabel("Radius in stellar radii")

#for theta in range(0, 91, 15):
##    calcOptDepthsOverTime("justBe_halfInj2", '/data/rob1g10/SPH/DATA/DISKS/justBe_halfInj2/processedDISKs/',
##                      434, 9, 0, theta)
#    calcOptDepthsOverTime("80d_0.0e", '/data/rob1g10/SPH/DATA/DISKS/80d_0.0e/processedDISKs/',
#                      61, 9, 0, theta)
#    print theta
#QDRANGE = [1,2]
#calcOptDepthsOverTime("justBe_halfInj2", '/data/rob1g10/SPH/DATA/DISKS/justBe_halfInj2/processedDISKs/',
#              325, 100, 0, 0, QDRANGE)
#calcOptDepthsOverTime("80d_0.0e", '/data/rob1g10/SPH/DATA/DISKS/80d_0.0e/processedDISKs/',
#              61, 100, 0, 0, QDRANGE)
#calcOptDepthsOverTime("80d_0.2e", '/data/rob1g10/SPH/DATA/DISKS/80d_0.2e/processedDISKs/',
#              61, 100, 0, 0, QDRANGE)
#calcOptDepthsOverTime("80d_0.4e", '/data/rob1g10/SPH/DATA/DISKS/80d_0.4e/processedDISKs/',
#              61, 100, 0, 0, QDRANGE)
#theta = phi = 0
#QDRANGE = [1,100]
#calcOptDepthsOverTime("60d_0.0e", '/data/rob1g10/SPH/DATA/DISKS/60d_0.0e/processedDISKs/',
#                  86, 9, 0, theta, QDRANGE)
#calcOptDepthsOverTime("60d_0.2e", '/data/rob1g10/SPH/DATA/DISKS/60d_0.2e/processedDISKs/',
#                  86, 9, 0, theta, QDRANGE)
#calcOptDepthsOverTime("60d_0.4e", '/data/rob1g10/SPH/DATA/DISKS/60d_0.4e/processedDISKs/',
#                  85, 9, 0, theta, QDRANGE)
#calcOptDepthsOverTime("40d_0.0e", '/data/rob1g10/SPH/DATA/DISKS/40d_0.0e/processedDISKs/',
#                  130, 9, 0, theta, QDRANGE)
#calcOptDepthsOverTime("40d_0.2e", '/data/rob1g10/SPH/DATA/DISKS/40d_0.2e/processedDISKs/',
#                  128, 9, 0, theta, QDRANGE)
#calcOptDepthsOverTime("40d_0.4e", '/data/rob1g10/SPH/DATA/DISKS/40d_0.4e/processedDISKs/',
#                  124, 9, 0, theta, QDRANGE)
#calcOptDepthsOverTime("20d_0.0e", '/data/rob1g10/SPH/DATA/DISKS/20d_0.0e/processedDISKs/',
#                  286, 9, 0, theta, QDRANGE)
#calcOptDepthsOverTime("20d_0.2e", '/data/rob1g10/SPH/DATA/DISKS/20d_0.2e/processedDISKs/',
#                  270, 9, 0, theta, QDRANGE)
#calcOptDepthsOverTime("20d_0.4e", '/data/rob1g10/SPH/DATA/DISKS/20d_0.4e/processedDISKs/',
#                  227, 9, 0, theta, QDRANGE)

#QDRANGE = [1,2]
#for theta in range(0, 91, 90):
#    calcOptDepthsOverTime("80d_0.4e_eqTest", '/data/rob1g10/SPH/DATA/DISKS/80d_0.4e_eqTest/processedDISKs/',
#                          356, 1, 0, theta, QDRANGE)

#QDRANGE = [1,2]
#for theta in range(0, 91, 100):
#    calcOptDepthsOverTime("SHEAR", "20d_0.0e_eqTest", '/data/rob1g10/SPH/DATA/DISKS/20d_0.0e_eqTest/processedDISKs/',
#                          387, 1, 0, theta, QDRANGE)

#QDRANGE = [1,8]
#for theta in range(0, 81, 10):
#    calcOptDepthsOverTime("injBurst_40d_0.4e", '/data/rob1g10/SPH/DATA/DISKS/injBurst_40d_0.4e/processedDISKs/',
#                          1, 24, 0, theta, QDRANGE)

#for i in range(10):         
#    QDRANGE = [1,8]
#    for theta in range(0, 91, 10):
#        calcOptDepthsOverTime("1e-11", "1e-11", "isoBe_fromBeginningTenth", '/data/rob1g10/SPH/DATA/DISKS/isoBe_fromBeginningTenth/processedDISKs/',
#                              411+5+i, 1, 0, theta, QDRANGE)
#QDRANGE = [1,2]
#for theta in range(0, 91, 10):
#        calcOptDepthsOverTime("1e-11UP", "1e-11UP", "isoBe_fromBeginningTenth", '/data/rob1g10/SPH/DATA/DISKS/isoBe_fromBeginningTenth/processedDISKs/',
#                              411, 1, 0, theta, QDRANGE)
#QDRANGE = [1,2]
#for theta in range(0, 91, 10):
#        calcOptDepthsOverTime("1e-11DOWN", "1e-11DOWN", "isoBe_fromBeginningTenth", '/data/rob1g10/SPH/DATA/DISKS/isoBe_fromBeginningTenth/processedDISKs/',
#                              411, 1, 0, theta, QDRANGE)
#QDRANGE = [1,35]
#for theta in range(0, 91, 10):
#    calcOptDepthsOverTime("SCIAMA", "SCIAMA", "20d_0.2e_eqTest", '/data/rob1g10/SPH/DATA/DISKS/20d_0.2e_eqTest/processedDISKs/',
#                          387, 6, 0, theta, QDRANGE)
#

#QDRANGE = [1,34]
#for theta in range(0, 91, 10):
#    calcOptDepthsOverTime("SCIAMA", "SCIAMA", "20d_0.2e_eqTest", '/data/rob1g10/SPH/DATA/DISKS/20d_0.2e_eqTest/processedDISKs/',
#                          387, 1, 0, theta, QDRANGE)
#                  
#QDRANGE = [1,25]
#for theta in range(0, 91, 10):
#    calcOptDepthsOverTime("SCIAMA", "SCIAMA", "20d_0.0e_eqTest", '/data/rob1g10/SPH/DATA/DISKS/20d_0.0e_eqTest/processedDISKs/',
#                          387, 1, 0, theta, QDRANGE)
#QDRANGE = [1,7]
#for theta in range(0, 91, 10):
#    calcOptDepthsOverTime("SCIAMA", "SCIAMA", "20d_0.0e_eqTest", '/data/rob1g10/SPH/DATA/DISKS/20d_0.0e_eqTest/processedDISKs/',
#                          388, 1, 0, theta, QDRANGE)
#   
#QDRANGE = [1,39]
#for theta in range(0, 91, 10):
#    calcOptDepthsOverTime("SCIAMA2", "SCIAMA2", "BH20d_0.2e", '/data/rob1g10/SPH/DATA/DISKS/BH20d_0.2e/processedDISKs/',
#                          300, 1, 0, theta, QDRANGE)  
#QDRANGE = [1,5]
#for theta in range(0, 91, 10):
#    calcOptDepthsOverTime("SCIAMA2", "SCIAMA4", "80d_0.4e_eqTest", '/data/rob1g10/SPH/DATA/DISKS/80d_0.4e_eqTest/processedDISKs/',
#                          341, 8, 0, theta, QDRANGE)
#                         
#QDRANGE = [1,5]
#for theta in range(0, 91, 10):
#    calcOptDepthsOverTime("SCIAMA2", "SCIAMA4", "80d_0.2e_eqTest", '/data/rob1g10/SPH/DATA/DISKS/80d_0.2e_eqTest/processedDISKs/',
#                          261, 8, 0, theta, QDRANGE)
#
#QDRANGE = [1,5]
#for theta in range(0, 91, 10):
#    calcOptDepthsOverTime("SCIAMA2", "SCIAMA4", "80d_0.0e_eqTest", '/data/rob1g10/SPH/DATA/DISKS/80d_0.0e_eqTest/processedDISKs/',
#                          261, 8, 0, theta, QDRANGE)

#QDRANGE = [1,15]
#for theta in range(0, 91, 10):
#    calcOptDepthsOverTime("5e-12", "5e-12", "40d_0.4e_eqTest", '/data/rob1g10/SPH/DATA/DISKS/40d_0.4e_eqTest/processedDISKs/',
#                          380, 4, 0, theta, QDRANGE)
#QDRANGE = [1,11]                          
#for theta in range(0, 91, 10):
#    calcOptDepthsOverTime("5e-12", "5e-12", "40d_0.2e_eqTest", '/data/rob1g10/SPH/DATA/DISKS/40d_0.2e_eqTest/processedDISKs/',
#                          380, 1, 0, theta, QDRANGE)
#
#QDRANGE = [1,11]
#for theta in range(0, 91, 10):
#    calcOptDepthsOverTime("5e-12", "5e-12", "40d_0.2e_eqTest", '/data/rob1g10/SPH/DATA/DISKS/40d_0.2e_eqTest/processedDISKs/',
#                          381, 1, 0, theta, QDRANGE)                         
#QDRANGE = [1,11]
#for theta in range(0, 91, 10):
#    calcOptDepthsOverTime("5e-12", "5e-12", "40d_0.2e_eqTest", '/data/rob1g10/SPH/DATA/DISKS/40d_0.2e_eqTest/processedDISKs/',
#                          382, 1, 0, theta, QDRANGE)
#QDRANGE = [1,4]
#for theta in range(0, 91, 10):
#    calcOptDepthsOverTime("5e-12", "5e-12", "40d_0.2e_eqTest", '/data/rob1g10/SPH/DATA/DISKS/40d_0.2e_eqTest/processedDISKs/',
#                          383, 1, 0, theta, QDRANGE)
#QDRANGE = [1,10]
#for theta in range(90, 91, 10):
#    calcOptDepthsOverTime("5e-12", "5e-12", "40d_0.0e_eqTest", '/data/rob1g10/SPH/DATA/DISKS/40d_0.0e_eqTest/processedDISKs/',
#                          380, 1, 0, theta, QDRANGE)
#QDRANGE = [1,10]
#for theta in range(90, 91, 10):
#    calcOptDepthsOverTime("5e-12", "5e-12", "40d_0.0e_eqTest", '/data/rob1g10/SPH/DATA/DISKS/40d_0.0e_eqTest/processedDISKs/',
#                          381, 1, 0, theta, QDRANGE)
#QDRANGE = [1,10]
#for theta in range(90, 91, 10):
#    calcOptDepthsOverTime("5e-12", "5e-12", "40d_0.0e_eqTest", '/data/rob1g10/SPH/DATA/DISKS/40d_0.0e_eqTest/processedDISKs/',
#                          382, 1, 0, theta, QDRANGE)
#QDRANGE = [1,7]
#for theta in range(90, 91, 10):
#    calcOptDepthsOverTime("5e-12", "5e-12", "40d_0.0e_eqTest", '/data/rob1g10/SPH/DATA/DISKS/40d_0.0e_eqTest/processedDISKs/',
#                          383, 1, 0, theta, QDRANGE)
#                          
#QDRANGE = [5,7]
#for theta in range(0, 91, 10):
#    calcOptDepthsOverTime("SCIAMA", "SCIAMA", "60d_0.4e_eqTest", '/data/rob1g10/SPH/DATA/DISKS/60d_0.4e_eqTest/processedDISKs/',
#                          581, 5, 0, theta, QDRANGE)
#                         
#QDRANGE = [1,5]
#for theta in range(0, 91, 10):
#    calcOptDepthsOverTime("SCIAMA", "SCIAMA", "60d_0.2e_eqTest", '/data/rob1g10/SPH/DATA/DISKS/60d_0.2e_eqTest/processedDISKs/',
#                          581, 6, 0, theta, QDRANGE)
#
#QDRANGE = [1,5]
#for theta in range(0, 91, 10):
#    calcOptDepthsOverTime("SCIAMA", "SCIAMA", "60d_0.0e_eqTest", '/data/rob1g10/SPH/DATA/DISKS/60d_0.0e_eqTest/processedDISKs/',
#                          581, 8, 0, theta, QDRANGE)
#
#QDRANGE = [1,66]
#for theta in range(0, 91, 10):
#    calcOptDepthsOverTime("BASICFAST", "TESTNEW", "20d_0.4e_eqTest", '/data/rob1g10/SPH/DATA/DISKS/20d_0.4e_eqTest/processedDISKs/',
#                          387, 3, 0, theta, QDRANGE)                        
#QDRANGE = [1,15]
#for theta in range(0, 91, 10):
#    calcOptDepthsOverTime("BASICFAST", "TESTNEW", "40d_0.4e_eqTest", '/data/rob1g10/SPH/DATA/DISKS/40d_0.4e_eqTest/processedDISKs/',
#                          380, 4, 0, theta, QDRANGE)
                          
#QDRANGE = [1,2]
#for theta in range(0, 91, 10):
#    calcOptDepthsOverTime("SCIAMA", "SCIAMA", "80d_0.4e_eqTest", '/data/rob1g10/SPH/DATA/DISKS/80d_0.4e_eqTest/processedDISKs/',
#                          387, 1, 0, theta, QDRANGE)
#
#QDRANGE = [1,2]
#for theta in range(0, 91, 10):
#    calcOptDepthsOverTime("SCIAMA", "SCIAMA", "80d_0.2e_eqTest", '/data/rob1g10/SPH/DATA/DISKS/80d_0.2e_eqTest/processedDISKs/',
#                          387, 1, 0, theta, QDRANGE)
##                          
#QDRANGE = [1,15]
#for theta in range(0, 91, 10):
#    calcOptDepthsOverTime("SCIAMA", "SCIAMA", "40d_0.4e_eqTest", '/data/rob1g10/SPH/DATA/DISKS/40d_0.4e_eqTest/processedDISKs/',
#                          380, 8, 0, theta, QDRANGE)
#                         
#QDRANGE = [1,2]
#for theta in range(0, 91, 100):
#    calcOptDepthsOverTime("BASICFAST", "TESTNEW", "40d_0.2e_eqTest", '/data/rob1g10/SPH/DATA/DISKS/40d_0.2e_eqTest/processedDISKs/',
#                          380, 1, 0, theta, QDRANGE)
##                         
#QDRANGE = [1,10]
#for theta in range(0, 91, 10):
#    calcOptDepthsOverTime("SCIAMA", "SCIAMA", "40d_0.0e_eqTest", '/data/rob1g10/SPH/DATA/DISKS/40d_0.0e_eqTest/processedDISKs/',
#                          380, 12, 0, theta, QDRANGE)
                         
#
#QDRANGE = [1,2]
#for theta in range(0, 91, 10):
#    calcOptDepthsOverTime("BASIC", "isoBe_fromBeginningTenth", '/data/rob1g10/SPH/DATA/DISKS/isoBe_fromBeginningTenth/processedDISKs/',
#                          411, 1, 0, theta, QDRANGE, dropbox='y')
#QDRANGE = [1,2]
#for theta in range(0, 91, 10):
#    calcOptDepthsOverTime("BASIC", "isoBe_fromBeginningTenth", 'DATA/',
#                          411, 1, 0, theta, QDRANGE, dropbox='y')
#QDRANGE = [1,2]
#for theta in range(0, 91, 10):
#    calcOptDepthsOverTime("BASICRTEMP2E5", "BASIC", "isoBe_fromBeginningTenth", '/data/rob1g10/SPH/DATA/DISKS/isoBe_fromBeginningTenth/processedDISKs/',
#                          411, 1, 0, theta, QDRANGE, dropbox='y')
#QDRANGE = [1,2]
#for theta in range(0, 91, 10):
#    calcOptDepthsOverTime("BASICNOBE", "isoBe_fromBeginningTenth", '/data/rob1g10/SPH/DATA/DISKS/isoBe_fromBeginningTenth/processedDISKs/',
#                          411, 1, 0, theta, QDRANGE, dropbox='y')
#QDRANGE = [1,2]
#for theta in range(0, 91, 10):
#    calcOptDepthsOverTime("TESTNEW", "isoBe_fromBeginningTenth", '/data/rob1g10/SPH/DATA/DISKS/isoBe_fromBeginningTenth/processedDISKs/',
#                          411, 1, 0, theta, QDRANGE, dropbox='y')
#QDRANGE = [1,2]
#for theta in range(0, 91, 10):
#    calcOptDepthsOverTime("TESTEXPH", "isoBe_fromBeginningTenth", '/data/rob1g10/SPH/DATA/DISKS/isoBe_fromBeginningTenth/processedDISKs/',
#                          411, 1, 0, theta, QDRANGE, dropbox='y')
#QDRANGE = [1,2]
#for theta in range(0, 91, 10):
#    calcOptDepthsOverTime("TESTEXP", "isoBe_fromBeginningTenth", '/data/rob1g10/SPH/DATA/DISKS/isoBe_fromBeginningTenth/processedDISKs/',
#                          411, 1, 0, theta, QDRANGE, dropbox='y')
#QDRANGE = [1,2]
#for theta in range(0, 91, 10):
#    calcOptDepthsOverTime("TESTRAD", "TESTNEW", "isoBe_fromBeginningTenth", '/data/rob1g10/SPH/DATA/DISKS/isoBe_fromBeginningTenth/processedDISKs/',
#                          411, 1, 0, theta, QDRANGE, dropbox='y')
#QDRANGE = [1,2]
#for theta in range(0, 91, 10):
#    calcOptDepthsOverTime("TESTBOTH", "isoBe_fromBeginningTenth", '/data/rob1g10/SPH/DATA/DISKS/isoBe_fromBeginningTenth/processedDISKs/',
#                          411, 1, 0, theta, QDRANGE, dropbox='y')
#QDRANGE = [1,2]
#for theta in range(0, 91, 10):
#    calcOptDepthsOverTime("TESTION", "isoBe_fromBeginningTenth", '/data/rob1g10/SPH/DATA/DISKS/isoBe_fromBeginningTenth/processedDISKs/',
#                          411, 1, 0, theta, QDRANGE, dropbox='y')
#QDRANGE = [1,2]
#for theta in range(0, 91, 10):
#    calcOptDepthsOverTime("BASICION6", "BASIC", "isoBe_fromBeginningTenth", '/data/rob1g10/SPH/DATA/DISKS/isoBe_fromBeginningTenth/processedDISKs/',
#                          411, 1, 0, theta, QDRANGE, dropbox='y')
#QDRANGE = [1,2]
#for theta in range(30, 91, 10):
#    calcOptDepthsOverTime("BASICION3", "BASIC", "isoBe_fromBeginningTenth", '/data/rob1g10/SPH/DATA/DISKS/isoBe_fromBeginningTenth/processedDISKs/',
#                          411, 1, 0, theta, QDRANGE, dropbox='y')
#QDRANGE = [1,2]
#for theta in range(70, 91, 10):
#    calcOptDepthsOverTime("BASICION1", "BASIC", "isoBe_fromBeginningTenth", '/data/rob1g10/SPH/DATA/DISKS/isoBe_fromBeginningTenth/processedDISKs/',
#                          411, 1, 0, theta, QDRANGE, dropbox='y')
#QDRANGE = [1,2]
#for theta in range(0, 91, 10):
#    calcOptDepthsOverTime("BASICION05", "BASIC", "isoBe_fromBeginningTenth", '/data/rob1g10/SPH/DATA/DISKS/isoBe_fromBeginningTenth/processedDISKs/',
#                          411, 1, 0, theta, QDRANGE, dropbox='y')
#QDRANGE = [1,2]
#for theta in range(0, 91, 10):
#    calcOptDepthsOverTime("BASICION02", "BASIC", "isoBe_fromBeginningTenth", '/data/rob1g10/SPH/DATA/DISKS/isoBe_fromBeginningTenth/processedDISKs/',
#                          411, 1, 0, theta, QDRANGE, dropbox='y')
#QDRANGE = [1,2]
#for theta in range(0, 91, 10):
#    calcOptDepthsOverTime("BASICION01", "BASIC", "isoBe_fromBeginningTenth", '/data/rob1g10/SPH/DATA/DISKS/isoBe_fromBeginningTenth/processedDISKs/',
#                          411, 1, 0, theta, QDRANGE, dropbox='y')
#QDRANGE = [1,2]
#for theta in range(0, 91, 10):
#    calcOptDepthsOverTime("BASICION005", "BASIC", "isoBe_fromBeginningTenth", '/data/rob1g10/SPH/DATA/DISKS/isoBe_fromBeginningTenth/processedDISKs/',
#                          411, 1, 0, theta, QDRANGE, dropbox='y')
#QDRANGE = [1,2]
#for theta in range(0, 91, 10):
#    calcOptDepthsOverTime("BASICION001", "BASIC", "isoBe_fromBeginningTenth", '/data/rob1g10/SPH/DATA/DISKS/isoBe_fromBeginningTenth/processedDISKs/',
#                          411, 1, 0, theta, QDRANGE, dropbox='y')
#QDRANGE = [1,2]
#for theta in range(0, 91, 10):
#    calcOptDepthsOverTime("BASICION0005", "BASIC", "isoBe_fromBeginningTenth", '/data/rob1g10/SPH/DATA/DISKS/isoBe_fromBeginningTenth/processedDISKs/',
#                          411, 1, 0, theta, QDRANGE, dropbox='y')
#QDRANGE = [1,2]
#for theta in range(0, 91, 10):
#    calcOptDepthsOverTime("BASICION0001", "BASIC", "isoBe_fromBeginningTenth", '/data/rob1g10/SPH/DATA/DISKS/isoBe_fromBeginningTenth/processedDISKs/',
#                         411, 1, 0, theta, QDRANGE, dropbox='y')
#QDRANGE = [1,2]
#for theta in range(0, 91, 10):
#    calcOptDepthsOverTime("TESTION10", "TESTNEW", "isoBe_fromBeginningTenth", '/data/rob1g10/SPH/DATA/DISKS/isoBe_fromBeginningTenth/processedDISKs/',
#                          411, 1, 0, theta, QDRANGE, dropbox='y')
#QDRANGE = [1,2]
#for theta in range(0, 91, 10):
#    calcOptDepthsOverTime("TESTION6", "TESTNEW", "isoBe_fromBeginningTenth", '/data/rob1g10/SPH/DATA/DISKS/isoBe_fromBeginningTenth/processedDISKs/',
#                          411, 1, 0, theta, QDRANGE, dropbox='y')
#QDRANGE = [1,2]
#for theta in range(0, 91, 10):
#    calcOptDepthsOverTime("TESTION3", "TESTNEW", "isoBe_fromBeginningTenth", '/data/rob1g10/SPH/DATA/DISKS/isoBe_fromBeginningTenth/processedDISKs/',
#                          411, 1, 0, theta, QDRANGE, dropbox='y')
#QDRANGE = [1,2]
#for theta in range(0, 91, 10):
#    calcOptDepthsOverTime("TESTION1", "TESTNEW", "isoBe_fromBeginningTenth", '/data/rob1g10/SPH/DATA/DISKS/isoBe_fromBeginningTenth/processedDISKs/',
#                          411, 1, 0, theta, QDRANGE, dropbox='y')
#QDRANGE = [1,2]
#for theta in range(0, 91, 10):
#    calcOptDepthsOverTime("TESTION05", "TESTNEW", "isoBe_fromBeginningTenth", '/data/rob1g10/SPH/DATA/DISKS/isoBe_fromBeginningTenth/processedDISKs/',
#                          411, 1, 0, theta, QDRANGE, dropbox='y')                          
#QDRANGE = [1,2]
#for theta in range(0, 91, 10):
#    calcOptDepthsOverTime("TESTION02", "TESTNEW", "isoBe_fromBeginningTenth", '/data/rob1g10/SPH/DATA/DISKS/isoBe_fromBeginningTenth/processedDISKs/',
#                          411, 1, 0, theta, QDRANGE, dropbox='y')
#QDRANGE = [1,2]
#for theta in range(0, 91, 10):
#    calcOptDepthsOverTime("TESTION01", "TESTNEW", "isoBe_fromBeginningTenth", '/data/rob1g10/SPH/DATA/DISKS/isoBe_fromBeginningTenth/processedDISKs/',
#                          411, 1, 0, theta, QDRANGE, dropbox='y')
#QDRANGE = [1,35]
#for theta in range(0, 91, 10):
#    calcOptDepthsOverTime("BASIC", "20d_0.2e_eqTest", '/data/rob1g10/SPH/DATA/DISKS/20d_0.2e_eqTest/processedDISKs/',
#                          387, 3, 0, theta, QDRANGE)
#QDRANGE = [1,25]
#for theta in range(0, 91, 10):
#    calcOptDepthsOverTime("BASIC", "20d_0.0e_eqTest", '/data/rob1g10/SPH/DATA/DISKS/20d_0.0e_eqTest/processedDISKs/',
#                          387, 5, 0, theta, QDRANGE)
#QDRANGE = [1,8]
#for theta in range(50, 91, 10):
#    calcOptDepthsOverTime("BASIC", "isoBe_fromBeginningTenth", '/data/rob1g10/SPH/DATA/DISKS/isoBe_fromBeginningTenth/processedDISKs/',
#                          411, 13, 0, theta, QDRANGE)                         

                          
#QDRANGE = [1,14]
#for theta in range(0, 91, 10):
#    calcOptDepthsOverTime("40d_0.4e_eqTest", '/data/rob1g10/SPH/DATA/DISKS/40d_0.4e_eqTest/processedDISKs/',
#                          380, 8, 0, theta, QDRANGE)

#QDRANGE = [1,2]
#for theta in [0, 90]:
#    calcOptDepthsOverTime("40d_0.4e_eqTest", '/data/rob1g10/SPH/DATA/DISKS/40d_0.4e_eqTest/processedDISKs/',
#                          3, 1, 0, theta, QDRANGE)


#T1 = np.array([])
##I1 = np.array([])
##for DISKnum in range(22, 57):
#for theta in range(0, 91, 100):
#        #filename = '/data/rob1g10/SPH/DATA/DISKS/20d_0.4e_eqTest/processedDISKs/QGDISK387{}'.format(DISKnum)
#        #diskname = filename[:-26] + "disk001"
#        #IPname = filename[:-11] + "IPDISK" + filename[-5:-2]
#        #Becoords, NScoords, partcoords, partvels, rho, u, num_parts, time, umass, udist, utime = readFILES(filename, diskname, IPname)
##    
##        baseName = "20d_0.4e_eqTest"    
##        CcodeName = "BASIC"
##        obsAngles = [theta, 0]
#    Tsum, partcoords1 = plotRotatingOptDepthFigsDropbox("../../Desktop/SciamaDATA/NEW/BASICOLD_theta{}.txt".format(theta), #"../../Desktop/SciamaDATA/NEW/Min0.9_R8-15_theta{}.txt".format(theta),
#                            #"/data/rob1g10/SPH/DATA/DISKS/{}/optDepthData{}{}_theta{}_phi{}_time{:0.1f}.txt".format(baseName, CcodeName, baseName, obsAngles[0], obsAngles[1], time * utime / (24 * 3600)),
#                            [filename, diskname, IPname], 
#                            "../Plots to Show/IsoBe data tests/optDepthZ/",#"../Plots to Show/All Round Image/isoBe_temp/", 
#                            "BASICOLDTOP", "OptDepth", [0,1,90], [0,11,15], [theta, 0], 15)
###    OD2, partcoords2 = plotRotatingOptDepthFigs([filename, diskname, IPname], 
###                             "../Plots to Show/All Round Image/20d_0.4e_eqTest/", 
###                             "20d_0.4e_eqTest", [0,91,90], [0,11,15], [80, 0])
##    T1 = np.append(T1, Tsum)
###    I1 = np.append(I1, np.sum(1-np.exp(-OD1)))
###    T1 = np.exp(-OD1)
###    T2 = np.exp(-OD2)
###    print T1[-1] #np.sum(T2)
###I1 = np.array([42, 43, 42, 42.5, 44, 46, 48, 47.5, 39, 24])
#I1 = np.array([ 91.92051094,  90.46173838,  88.29814729,  88.47741181,  89.37006442, 89.96901982,  89.9781312,   89.2090977,   87.90792666,  80.54662987])
##plt.plot(range(0, 91, 10), T1 / 1000, label="Simulation")
##print T1 / 1000
#plt.plot(range(0, 91, 10), I1, label="Simulation")
#plt.legend(loc="best")
#plt.xlabel("inclination angle")
#plt.ylabel("-Equivalent Width")
#plt.xlabel("optical depth")
#plt.ylabel("Z")

#for theta in [0, 30, 60, 90]:
# plotRotatingOptDepthHistsDropbox("../../Desktop/SciamaDATA/NEW/BASICOLD_theta",
#                             [filename, diskname, IPname],  [0,91,30], ["blue", "green", "red", "black"], show="y")
#plt.show()

#for obsAngle in [0, 90]:
#    for theta in [0, 90, 180]:
#        for phi in [0, 90]:
#            plotRotatingOptDepthFigs([filename, diskname, IPname], 
#                                     "../Plots to Show/", 
#                                     "isoBe_close2equilibrium", 
#                                     [theta,theta+1,15], [phi,phi+1,15], [obsAngle, 0])
#print "obsAngle =", theta

"""----------------------------------------------------------------------------

            CREATING DATA FOR SINGLE FILES

----------------------------------------------------------------------------"""

theta = 0
phi = 0

#filename = '/data/rob1g10/SPH/DATA/DISKS/20d_0.0e/processedDISKs/QGDISK10703'
#diskname = '/data/rob1g10/SPH/DATA/DISKS/20d_0.0e/disk001'
#IPname = '/data/rob1g10/SPH/DATA/DISKS/20d_0.0e/processedDISKs/IPDISK107'
#
#baseName = "20d_0.0e"
#Becoords, NScoords, partcoords, partvels, rho, u, num_parts, time, umass, udist, utime = readFILES(filename, diskname, IPname)
#createOptDepthData4c(partcoords * (udist*0.01), partvels * (udist * 0.01)/utime, u * (udist * 0.01)**2/utime**2, num_parts, 
#                     1e-15 * 0.001 * umass / 2e30, (7 * 6.96e10), 2.5, theta, phi)
#os.system("./opticalDepth")
#name = "{}_theta{}_phi{}_time{:0.1f}".format(baseName, theta, phi, time * utime / (24 * 3600))
#os.system("mv optDepthData4py.txt /data/rob1g10/SPH/DATA/DISKS/{}/optDepthData{}.txt".format(baseName, name))
#
#
#filename = '/data/rob1g10/SPH/DATA/DISKS/40d_0.0e/processedDISKs/QGDISK04901'
#diskname = '/data/rob1g10/SPH/DATA/DISKS/40d_0.0e/disk001'
#IPname = '/data/rob1g10/SPH/DATA/DISKS/40d_0.0e/processedDISKs/IPDISK049'
#
#baseName = "40d_0.0e"
#Becoords, NScoords, partcoords, partvels, rho, u, num_parts, time, umass, udist, utime = readFILES(filename, diskname, IPname)
#createOptDepthData4c(partcoords * (udist*0.01), partvels * (udist * 0.01)/utime, u * (udist * 0.01)**2/utime**2, num_parts, 
#                     1e-15 * 0.001 * umass / 2e30, (7 * 6.96e10), 2.5, theta, phi)
#os.system("./opticalDepth")
#name = "{}_theta{}_phi{}_time{:0.1f}".format(baseName, theta, phi, time * utime / (24 * 3600))
#os.system("mv optDepthData4py.txt /data/rob1g10/SPH/DATA/DISKS/{}/optDepthData{}.txt".format(baseName, name))

#filename = '/data/rob1g10/SPH/DATA/DISKS/60d_0.0e/processedDISKs/QGDISK03302'
#diskname = '/data/rob1g10/SPH/DATA/DISKS/60d_0.0e/disk001'
#IPname = '/data/rob1g10/SPH/DATA/DISKS/60d_0.0e/processedDISKs/IPDISK033'
#
#baseName = "60d_0.0e"
#Becoords, NScoords, partcoords, partvels, rho, u, num_parts, time, umass, udist, utime = readFILES(filename, diskname, IPname)
#createOptDepthData4c(partcoords * (udist*0.01), partvels * (udist * 0.01)/utime, u * (udist * 0.01)**2/utime**2, num_parts, 
#                     1e-15 * 0.001 * umass / 2e30, (7 * 6.96e10), 2.5, theta, phi)
#os.system("./opticalDepthTime")
#name = "{}_theta{}_phi{}_time{:0.1f}".format(baseName, theta, phi, time * utime / (24 * 3600))
#os.system("mv optDepthData4pyTime.txt /data/rob1g10/SPH/DATA/DISKS/{}/optDepthData{}.txt".format(baseName, name))


#filename = '/data/rob1g10/SPH/DATA/DISKS/80d_0.0e/processedDISKs/QGDISK02309'
#diskname = '/data/rob1g10/SPH/DATA/DISKS/80d_0.0e/disk001'
#IPname = '/data/rob1g10/SPH/DATA/DISKS/80d_0.0e/processedDISKs/IPDISK023'
#
#baseName = "80d_0.0e"
#Becoords, NScoords, partcoords, partvels, rho, u, num_parts, time, umass, udist, utime = readFILES(filename, diskname, IPname)
#createOptDepthData4c(partcoords * (udist*0.01), partvels * (udist * 0.01)/utime, u * (udist * 0.01)**2/utime**2, num_parts, 
#                     1e-15 * 0.001 * umass / 2e30, (7 * 6.96e10), 2.5, theta, phi)
#os.system("./opticalDepthTime")
#name = "{}_theta{}_phi{}_time{:0.1f}".format(baseName, theta, phi, time * utime / (24 * 3600))
#os.system("mv optDepthData4pyTime.txt /data/rob1g10/SPH/DATA/DISKS/{}/optDepthData{}.txt".format(baseName, name))

"""----------------------------------------------------------------------------

            OPTICAL DEPTH DATA CREATION FOR EQUIVALENT WIDTH

----------------------------------------------------------------------------"""
#phi = 90
#for theta in range(0, 91, 15):
#    orbits_0e = [20,40,60,80]
#    DISKs_0e = ["28604", "13001", "08606", "06101"]
#    orbits_2e = [20,40,60,80]
#    DISKs_2e = ["27004", "12804", "08606", "06101"]
#    orbits_4e = [20,40,60,80]
#    DISKs_4e = ["22706", "12404", "08504", "06005"]
#    
#    allOrbits = orbits_0e + orbits_2e + orbits_4e
#    allDISKs = DISKs_0e + DISKs_2e + DISKs_4e
#    
#    
#    for i in range(len(allOrbits)):
#        eccentricity = 2 * (i / 4)
#        filename = '/data/rob1g10/SPH/DATA/DISKS/{}/processedDISKs/QGDISK{}'.format(str(allOrbits[i])+"d_0."+str(eccentricity)+"e", allDISKs[i])
#        diskname = filename[:-26] + "disk001"
#        IPname = filename[:-11] + "IPDISK" + filename[-5:-2]
#        Becoords, NScoords, partcoords, partvels, rho, u, num_parts, time, umass, udist, utime = readFILES(filename, diskname, IPname)
#        createOptDepthData4c(partcoords * (udist*0.01), partvels * (udist * 0.01)/utime, u * (udist * 0.01)**2/utime**2, num_parts, 
#                                     1e-15 * 0.001 * umass / 2e30, (7 * 6.96e10), 2.5, theta, phi)
#        os.system("./opticalDepthFast")
#        baseName = filename[29:37]
#        name = "{}_theta{}_phi{}_time{:0.1f}".format(baseName, theta, phi, time * utime / (24 * 3600))
#        os.system("mv optDepthData4pyFast.txt /data/rob1g10/SPH/DATA/DISKS/{}/optDepthData{}.txt".format(baseName, name))

"""----------------------------------------------------------------------------

            CALCULATING EQUIVALENT WIDTH

----------------------------------------------------------------------------"""

#filename = '/data/rob1g10/SPH/DATA/DISKS/20d_0.0e/processedDISKs/QGDISK10703'
#diskname = '/data/rob1g10/SPH/DATA/DISKS/20d_0.0e/disk001'
#IPname = '/data/rob1g10/SPH/DATA/DISKS/20d_0.0e/processedDISKs/IPDISK107'
#
#filename = '/data/rob1g10/SPH/DATA/DISKS/40d_0.0e/processedDISKs/QGDISK04901'
#diskname = '/data/rob1g10/SPH/DATA/DISKS/40d_0.0e/disk001'
#IPname = '/data/rob1g10/SPH/DATA/DISKS/40d_0.0e/processedDISKs/IPDISK049'
#
#filename = '/data/rob1g10/SPH/DATA/DISKS/60d_0.0e/processedDISKs/QGDISK03302'
#diskname = '/data/rob1g10/SPH/DATA/DISKS/60d_0.0e/disk001'
#IPname = '/data/rob1g10/SPH/DATA/DISKS/60d_0.0e/processedDISKs/IPDISK033'
#
#filename = '/data/rob1g10/SPH/DATA/DISKS/80d_0.0e/processedDISKs/QGDISK02309'
#diskname = '/data/rob1g10/SPH/DATA/DISKS/80d_0.0e/disk001'
#IPname = '/data/rob1g10/SPH/DATA/DISKS/80d_0.0e/processedDISKs/IPDISK023'
#
#baseName = filename[29:37]
#theta = phi = 90
#Becoords, NScoords, partcoords, partvels, rho, u, num_parts, time, umass, udist, utime = readFILES(filename, diskname, IPname)
#oD = readOptDepthData4py("/data/rob1g10/SPH/DATA/DISKS/{}/optDepthData{}_theta{}_phi{}_time{:0.1f}.txt".format(baseName, baseName, theta, phi, time * utime / (24 * 3600)))
#T = np.exp(-oD) / 1e-3
#print(calcEW(partvels, T, np.ones(len(partvels))*0.01, 600, 600, udist, utime, intSteps=200))


"""----------------------------------------------------------------------------

            MAC

----------------------------------------------------------------------------"""

#phi = 0
#for theta in range(0, 9):
#    filename = '/Users/rorybrown/Dropbox/PhD/DISKs/QGDISK28201'
#    diskname = '/Users/rorybrown/Dropbox/PhD/DISKs/disk001'
#    IPname = '/Users/rorybrown/Dropbox/PhD/DISKs/IPDISK282'
#    Becoords, NScoords, partcoords, partvels, rho, u, num_parts, time, umass, udist, utime = readFILES(filename, diskname, IPname)
#    name = "DISKs/optDepthData20d_0.4e_theta{}_phi{}_time2051.2_zoomY.png".format(theta, phi)
#    oD = readOptDepthData4py("/Users/rorybrown/Dropbox/PhD/DISKs/optDepthData20d_0.4e_theta{}_phi{}_time2051.2.txt".format(theta, phi))
#    T = np.exp(-oD) / 1e-3
#    print(theta, phi)
#    print(len(T))
#    plotAndSaveSpectralLine(partvels, T, [phi, theta], 200, udist, 
#            utime, "uniform", name, save='y', xlims="n")