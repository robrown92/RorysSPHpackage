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
    
    f = open("optDepthParams4c{}.txt".format(CcodeName), "w")
    f.write("{}\n".format(numparts))
    f.write("{:.10}\n".format(partmass))
    f.write("{:.10}\n".format(BeRadius))
    f.write("{:.10}\n".format(BeMass))
    f.write("{:.10}\n".format(powerLawIndex))
    f.write("{:.10}\n".format(baseDensity))
    f.write("{:.10}\n".format(theta))
    f.close()
    return


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

def calcOptDepthsOverTime(CcodeName, dataCreationName, baseName, inputdirectory, startDISKs, 
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
        createOptDepthData4c(dataCreationName, #CcodeName,
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
    
    return
    

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


def opticalDepthPython(partcoords, partvels, u, rho, numparts, partmass, BeRadius, 
                         powerLawIndex, theta, phi, intSteps=10):
                           
    opticalDepth = np.zeros(numparts)
    
    CUTOFF = 15
    
    m_e = 9.1e-28 #"mass of electron"
    e = 4.8e-10 #"electron charge"
    lambda0 = 656.28e-7 #"rest wavelength"
    c = 3e10 #"speed of light"
    f = 0.6407 #"absorption oscillator strength of the transition. for n=3 to n=2 Balmer line (emitting at 656nm), this is 3.03x10^-19"
    n_rs = 6.022e18 #"population density of level 2 electrons at stellar radius"
    kb = 1.38e-23 #"Boltzmann constant"
    idealGasConstant = 8.314 #"ideal gas constant"
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
    H = (R / BeRadius) ** 1.5
    
    """
    calculate thermal velocity for all particles
    temperature -> U = 0.5 * (3 * kb * T) _-> 2 * U = 3 * kb * T
    thermal velocity = sqrt (3 * kb * T / m) = sqrt (2 * U / m)
    """
    v_th = np.sqrt(2 * u / partmass)
    
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
            if(np.sqrt(xcoords[i]**2 + ycoords[i]**2) < BeRadius and zcoords[i] < 0
                or intSum > CUTOFF):
                intSum = CUTOFF
#                behindStarCount += 1
                break
            a = z_D[i] + (j*stepSize) # integration lower limit
            b = a + stepSize # integration upper limit
            
            "find nearest particle to end a and get values"
            partA = np.argmin( np.sqrt( (xcoords - xcoords[i])**2 + (ycoords - ycoords[i])**2 + (zcoords - a)**2))
#            print i, partA
#            if(i%10000 == 0):
#                print "{}%".format(i/ 1000)
            "calc modified population density"
            n = (1-ionFrac[partA]) * 6.022e22 * rho[partA]
            "add this to the sum"
            intSum += (b-a) * coefficient * n * np.exp(-0.5 * (v[partA] - v_D[partA]) ** 2  / (v_th[partA] ** 2)) / v_th[partA]  
    
        opticalDepth[i] = intSum
#    print behindStarCount
    f = open("optDepthsFromPy.txt", "w")
    for i in range(numparts):
        f.write("{}\n".format(opticalDepth[i]))
    f.close()
    
    return opticalDepth

#opticalDepthPython(partcoords * udist, partvels * udist/utime, u * udist**2/utime**2, 
#                   rho * umass / udist**3, numparts, 1e-15 * umass / 2e30, (7 * 6.96e10), 
#                         2.5, theta, phi)

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

#Becoords, NScoords, partcoords, partvels, rho, u, num_parts, time, umass, udist, utime = readFILES(filename, diskname, IPname)
#print max(rho) * umass / udist**3
#
#sortOrder = np.argsort(partcoords[:,2])[::-1]
#print sortOrder
#
#partcoords = partcoords[sortOrder]
#partvels = partvels[sortOrder]
#rho = rho[sortOrder]

" calcs by hand -------------------------------------------------------------"

"constants"
#theta = np.deg2rad(0)
#T = 1.3e4
#m_e = 9.1e-28 # mass of electron
#e = 4.8e-10 # electron charge
#lambda0 = 656.28e-7 # rest wavelength
#c = 3e10 # speed of light
#nu0 = c / lambda0 # rest frequency
#f = 0.6407 # absorption oscillator strength of the transition. for n=3 to n=2 Balmer line (emitting at 656nm), this is 0.6407
#kb = 1.38e-16 # Boltzmann constant 
#h = 6.626e-27 # plank constant 
#amu = 1.66e-24 # atomic mass unit
#G = 6.67e-8 # gravitational constant
#M_stellar = 18 * 2e33

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
#ionCoeff = (2 * np.pi * m_e * kb * T) / (h*h)
#ionExp = np.exp( (-13.6 * 1.6e-12) / (kb * T))
#ionNumDens = 6.022e23 * 1e-12#rho
#ionA = 1
#ionB = (1 / ionNumDens) * ionCoeff**1.5 * ionExp
#ionC = -ionB
#ionFrac = (-ionB + np.sqrt(ionB*ionB - 4*ionA*ionC)) / (2*ionA)
#print 1 - ionFrac
#n = ionFrac * 6.022 * pow(10, 23) * rho * umass / udist**3

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

#QDRANGE = [1,66]
#for theta in range(0, 91, 10):
#    calcOptDepthsOverTime("BASIC", "20d_0.4e_eqTest", '/data/rob1g10/SPH/DATA/DISKS/20d_0.4e_eqTest/processedDISKs/',
#                          387, 2, 0, theta, QDRANGE)

#QDRANGE = [1,2]
#for theta in range(0, 91, 10):
#    calcOptDepthsOverTime("BASIC", "BASIC", "isoBe_fromBeginningTenth", '/data/rob1g10/SPH/DATA/DISKS/isoBe_fromBeginningTenth/processedDISKs/',
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
#    calcOptDepthsOverTime("TESTRAD", "isoBe_fromBeginningTenth", '/data/rob1g10/SPH/DATA/DISKS/isoBe_fromBeginningTenth/processedDISKs/',
#                         411, 1, 0, theta, QDRANGE, dropbox='y')
#QDRANGE = [1,2]
#for theta in range(0, 91, 10):
#    calcOptDepthsOverTime("TESTALL", "TESTNEW", "isoBe_fromBeginningTenth", '/data/rob1g10/SPH/DATA/DISKS/isoBe_fromBeginningTenth/processedDISKs/',
#                          411, 1, 0, theta, QDRANGE, dropbox='y')
#QDRANGE = [1,2]
#for theta in range(0, 91, 10):
#    calcOptDepthsOverTime("TESTBOTH", "TESTNEW", "isoBe_fromBeginningTenth", '/data/rob1g10/SPH/DATA/DISKS/isoBe_fromBeginningTenth/processedDISKs/',
#                          411, 1, 0, theta, QDRANGE, dropbox='y')
#QDRANGE = [1,2]
#for theta in range(0, 91, 10):
#    calcOptDepthsOverTime("TESTION", "isoBe_fromBeginningTenth", '/data/rob1g10/SPH/DATA/DISKS/isoBe_fromBeginningTenth/processedDISKs/',
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

#                          
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
#for theta in range(0, 91, 10):
#    Tsum, partcoords1 = plotRotatingOptDepthFigs("BASIC", [filename, diskname, IPname], 
#                             "../Plots to Show/All Round Image/isoBe_fromBeginningTenth/", 
#                             "isoBe_fromBeginningTenth", [0,91,90], [0,11,15], [theta, 0], 25)
##    OD2, partcoords2 = plotRotatingOptDepthFigs([filename, diskname, IPname], 
##                             "../Plots to Show/All Round Image/20d_0.4e_eqTest/", 
##                             "20d_0.4e_eqTest", [0,91,90], [0,11,15], [80, 0])
#    T1 = np.append(T1, Tsum)
##    I1 = np.append(I1, np.sum(1-np.exp(-OD1)))
##    T1 = np.exp(-OD1)
##    T2 = np.exp(-OD2)
#    print T1[-1] #np.sum(T2)
#plt.plot(range(0, 91, 10), T1)
##plt.plot(range(0, 91, 10), I1)
#plt.xlabel("inclination angle")
#plt.ylabel("-EW")

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