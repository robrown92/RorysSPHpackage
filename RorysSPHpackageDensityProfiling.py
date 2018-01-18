# -*- coding: utf-8 -*-
"""
Created on Wed Feb 10 15:35:37 2016

@author: rob1g10
"""

import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import pylab as plt
import os
from scipy.optimize import leastsq, curve_fit
from RorysSPHpackage import readQDISK, readUnits, readFILES

"""----------------------------------------------------------------------------

                    DENSITY PROFILE PLOTTING

----------------------------------------------------------------------------"""


def plotDensityProfileOLD(partcoords, rhos, udist, umass, ylims=[-1e-6, 1e-6], tol=1e-3):
    """
    Produces a plot of the density profile for the Be star disk.

    GEOMETRY
    --------

    z axis is height of disk, y axis is width and x axis is depth.

    INPUT VARIABLES
    ---------------

    partcoords ((npart, 3) float array) - x,y,z coordinates of the disk
                                         particles.

    rhos (npart float array) - density values of the disk particles.

    udist, umass (floats) - units used in simulation

    ylim (2 element float) - defines limits of y axis.

    tol (float) - tolerance either side of x=0 that the function will use particles for
          plotting.

    OUTPUT VARIABLES
    ----------------

    A plot of the disk profile.

    """

    " Initialise arrays "

    ycoordsPLOT = np.array([])
    rhosPLOT = np.array([])

    " Get only particles in one slice of the disk i.e. x=0 "

    for i in range(len(partcoords[:,0])):
        if(np.abs(partcoords[i,0]) <= tol and partcoords[i,1] > 0):
           #and np.abs(partcoords[i,2]) <= tol):
            ycoordsPLOT = np.append(ycoordsPLOT, partcoords[i,1])
            rhosPLOT = np.append(rhosPLOT, rhos[i])
    
    """ Conversion of units """
    rhosPLOT *= umass/(udist**3)
    ycoordsPLOT *= udist / (7 * 695700e3 * 100)  
    
    """ Density Scatter """
    
    p = np.polyfit(np.log10(ycoordsPLOT), np.log10(rhosPLOT), 1)
#    print p
    
#    plt.plot(np.log10(ycoordsPLOT), p[0] * np.log10(ycoordsPLOT) + p[1])
    plt.scatter(ycoordsPLOT, np.log10(rhosPLOT), s=0.75, color=ylims[0], label=ylims[1] +"gradient = {:.3}".format(p[0]))
#    plt.xlim(0, 1.5)
#    plt.ylim(ylims[0], ylims[1])
    plt.xlabel("Distance from centre of Be star in Be radii")
    plt.ylabel("log of g/cm^3")
    plt.legend()
    #plt.title("Plot showing density profile and power law fit (gradient of {:.3})".format(p[0]))
#    plt.title("Density Profile for {}".format(ylims[1]))
    #plt.show()
    return
    
def plotDensityProfile(filenames, ylims=[-1e-6, 1e-6], tol=1e-3, colours = ['b', 'g', 'r', 'turquoise', 'black']):
    """
    Produces a plot of the density profile for the Be star disk.

    GEOMETRY
    --------

    z axis is height of disk, y axis is width and x axis is depth.

    INPUT VARIABLES
    ---------------

    partcoords ((npart, 3) float array) - x,y,z coordinates of the disk
                                         particles.

    rhos (npart float array) - density values of the disk particles.

    udist, umass (floats) - units used in simulation

    ylim (2 element float) - defines limits of y axis.

    tol (float) - tolerance either side of x=0 that the function will use particles for
          plotting.

    OUTPUT VARIABLES
    ----------------

    A plot of the disk profile.

    """

    " Initialise arrays "
    count = 0
    Ms = [-9, -10]
    for filename in filenames:
#        diskname = 'INPUT/disk001'
#        IPname = 'INPUT/IPDISK411'
#        diskname = '../../Desktop/40d_0.0e_eqTest/disk001'
#        IPname = '../../Desktop/40d_0.0e_eqTest/IPDISK' + filename[-5:-2]
        diskname = filename[:-26] + "disk001"
        IPname = filename[:-11] + "IPDISK" + filename[-5:-2]
        Becoords, NScoords, partcoords, partvels, rho, u, num_parts, time, umass, udist, utime = readFILES(filename, diskname, IPname)
        print max(rho * umass/(udist**3) )
        R = np.sqrt(partcoords[:,0]**2 + partcoords[:,1]**2)
        Rs = np.array([])
        rhos = np.array([])
        for i in range(len(R)):
            if(R[i] > 1.0 * 0.0844):
                Rs = np.append(Rs, R[i])
                rhos = np.append(rhos, rho[i])
        
        """ Conversion of units """
        rhos *= umass/(udist**3)
        Rs *= udist / (7 * 695700e3 * 100)
        
        Rs *= 7
        
        """ Density Scatter """
        
        #p = np.polyfit(np.log10(Rs), np.log10(rhos), 1)
    #    print p
        
        plt.scatter(Rs, np.log10(rhos), s=0.75, color=colours[count], label="{}".format(filename[-33:-36])) #label=ylims[1] +"gradient = {:.3}".format(p[0]))
        #plt.plot(np.log10(Rs), p[0] * np.log10(Rs) + p[1], color='r')
    #    plt.xlim(0, 1.5)
    #    plt.ylim(ylims[0], ylims[1])
        
    #    plt.title("Plot showing density profile and power law fit (gradient of {:.3})".format(p[0]))
    #    plt.title("Plot showing density profile for {}".format(ylims[1]))
        #plt.show()
#        plt.text(6, -9 - (count*0.5), "{}, rho0={:0.3}".format(Ms[count], max(rho * umass/(udist**3) )), color=colours[count])
        count += 1 
    plt.xlabel(r"radius (solar radii)")
    plt.ylabel(r"log density (g cm$^{-3}$)")
    plt.legend(loc="best")
    plt.xlim(1,10)
    plt.show()
    return

def plotEquatorialDensityProfile(filenames, ylims=[-1e-6, 1e-6], tol=1e-4):
    """
    Produces a plot of the density profile for the Be star disk.

    GEOMETRY
    --------

    z axis is height of disk, y axis is width and x axis is depth.

    INPUT VARIABLES
    ---------------

    partcoords ((npart, 3) float array) - x,y,z coordinates of the disk
                                         particles.

    rhos (npart float array) - density values of the disk particles.

    udist, umass (floats) - units used in simulation

    ylim (2 element float) - defines limits of y axis.

    tol (float) - tolerance either side of x=0 that the function will use particles for
          plotting.

    OUTPUT VARIABLES
    ----------------

    A plot of the disk profile.

    """

    " Initialise arrays "
    colour = 0
    start = 0
    for filename in filenames[start::4]:
#        diskname = 'INPUT/disk001'
#        IPname = 'INPUT/IPDISK411'
        diskname = filename[:-26] + "disk001"
        IPname = filename[:-11] + "IPDISK" + filename[-5:-2]
        Becoords, NScoords, partcoords, partvels, rho, u, num_parts, time, umass, udist, utime = readFILES(filename, diskname, IPname)

        R = np.sqrt(partcoords[:,0]**2 + partcoords[:,1]**2) * udist / (7 * 695700e3 * 100)  
        Rs = np.array([])
        rhos = np.array([])
        for i in range(len(R)):
            if(R[i] > 1. and np.abs(partcoords[i,2]) < tol):
                Rs = np.append(Rs, R[i])
                rhos = np.append(rhos, rho[i])
        
        """ Conversion of units """
        rhos *= umass/(udist**3)
        #Rs *= udist / (7 * 695700e3 * 100)  
        
        """ Density Scatter """
        
        #p = np.polyfit(np.log10(Rs), np.log10(rhos), 1)
    #    print p
        
        sortOrder = np.argsort(Rs)         
        
        plt.plot(Rs[sortOrder], np.log10(rhos)[sortOrder],  label = [r'1.4M$_{\odot}$, e=0.0',
                                                                   r'1.4M$_{\odot}$, e=0.2',
                                                                   r'1.4M$_{\odot}$, e=0.4',
                                                                   r'1.4M$_{\odot}$, e=0.6',
                      
                                                                   r'6M$_{\odot}$, e=0.0',
                                                                   r'6M$_{\odot}$, e=0.2',
                                                                   r'6M$_{\odot}$, e=0.4',
                                                                   r'6M$_{\odot}$, e=0.6',
                     
                                                                   r'10M$_{\odot}$, e=0.0',
                                                                   r'10M$_{\odot}$, e=0.2',
                                                                   r'10M$_{\odot}$, e=0.4',
                                                                   r'10M$_{\odot}$, e=0.6',][start::4][colour], 
                   color=['b', 'g', 'r', 'y', 'black', 'brown', 'turquoise', 'pink', 'purple'][colour]) #label=ylims[1] +"gradient = {:.3}".format(p[0]))
        colour += 1
        #plt.plot(np.log10(Rs), p[0] * np.log10(Rs) + p[1], color='r')
#        plt.xlim(0, 1.5)
    #    plt.ylim(ylims[0], ylims[1])
        plt.xlabel(r"R (R$_\odot$)")
        plt.ylabel(r"log density (g cm$^{-3}$)")
        plt.legend(loc='lower left')
    #    plt.title("Plot showing density profile and power law fit (gradient of {:.3})".format(p[0]))
    #    plt.title("Plot showing density profile for {}".format(ylims[1]))
        #plt.show()
    return Rs, rhos

def plotDensityProfileRiviniusandCarciofi(partcoords, rhos, udist, umass, ylims=[-1e-6, 1e-6], tol=5e-5):
    """
    Produces a plot of the density profile for the Be star disk to check
    Rivinius and Carciofi paper power law:
    
    rho = rho0 * (r / R*)^-n
    
    n = 2.5 - 4 is expected
    
    rho0 of around 10^-10 - 10^-13 is expected.

    GEOMETRY
    --------

    z axis is height of disk, y axis is width and x axis is depth.

    INPUT VARIABLES
    ---------------

    partcoords ((npart, 3) float array) - x,y,z coordinates of the disk
                                         particles.

    rhos (npart float array) - density values of the disk particles.

    udist, umass (floats) - units used in simulation

    ylim (2 element float) - defines limits of y axis.

    tol (float) - tolerance either side of x=0 that the function will use particles for
          plotting.

    OUTPUT VARIABLES
    ----------------

    A plot of the disk profile.

    """

    " Initialise arrays "

    radcoordsPLOT = np.array([])
    rhosPLOT = np.array([])

    " Get only particles in one slice of the disk i.e. x=0 "

    for i in range(len(partcoords[:,0])):
        if(np.abs(partcoords[i,2]) <= tol):
           #and np.abs(partcoords[i,2]) <= tol):
            radcoordsPLOT = np.append(radcoordsPLOT, np.sqrt(partcoords[i,1]**2 + partcoords[i,0]**2))
            rhosPLOT = np.append(rhosPLOT, rhos[i])
    
    """ Conversion of units """
    rhosPLOT *= umass/(udist**3)
    radcoordsPLOT *= udist / (7 * 695700e3 * 100)  
    
    """ Density Scatter """
    
    p = np.polyfit(np.log10(radcoordsPLOT), np.log10(rhosPLOT), 1)
    print p
    
    plt.plot(np.log10(radcoordsPLOT), p[0] * np.log10(radcoordsPLOT) + p[1])
    plt.scatter(np.log10(radcoordsPLOT), np.log10(rhosPLOT), s=0.75, color=ylims[0], label=ylims[1]+"gradient = {:.3}".format(p[0]))
#    plt.xlim(0, 1.5)
#    plt.ylim(ylims[0], ylims[1])
    plt.xlabel("log of distance from centre of Be star in Be Radii")
    plt.ylabel("log of kg/m^3")
    plt.legend()
    #plt.title("Plot showing density profile and power law fit (gradient of {:.3})".format(p[0]))
#    plt.title("Plot showing density profile for 24 & 83 day period")
    #plt.show()
    return
    

def plotDensityProfileTouhamiandGies(partcoords, rhos, u, udist, umass, utime, ylims=[-1e-6, 1e-6], tol=5e-6):
    """
    Produces a plot of the density profile for the Be star disk to check
    Touhami and Gies paper (2011) power law:
    
    rho = rho0 * R^-n exp(-0.5*(Z / H(R))^2)
    
    where H(R) = (cs / Vk) * R^ 3/2
    
    cs is speed of sound, Vk is keplerian velocity, rho0 is base gas density,
    R is radial coordinate and Z is vertical cylindrical coordinate.
    
    n = 2.5 - 4 is expected
    
    rho0 of around 10^-10 - 10^-13 is expected.

    GEOMETRY
    --------

    z axis is height of disk, y axis is width and x axis is depth.

    INPUT VARIABLES
    ---------------

    partcoords ((npart, 3) float array) - x,y,z coordinates of the disk
                                         particles.

    rhos (npart float array) - density values of the disk particles.
    
    u (npart float array) - internal energy of particles.

    udist, umass (floats) - units used in simulation

    ylim (2 element float) - defines limits of y axis.

    tol (float) - tolerance either side of x=0 that the function will use particles for
          plotting.

    OUTPUT VARIABLES
    ----------------

    A plot of the disk profile.

    """
    
    rhosPLOT = rhos * umass/(udist**3)
    radcoordsPLOT = np.sqrt(partcoords[:,1]**2 + partcoords[:,0]**2) * udist / (7 * 695700e3)
    zcoordsPLOT = partcoords[:,2] * udist / (7 * 695700e3)
    uPLOT = u * umass * udist**2 / utime**2
    
    fit = np.array([])
    breaks = np.array([0])
    radcoordsFIT = np.array([])
    rhosFIT = np.array([])
    uFIT = np.array([])
    zcoordsFIT = np.array([])
    for j in range(10):
        for i in range(len(partcoords[:,0])):
            if(np.abs(partcoords[i,2]) <= i*0.1*max(zcoordsPLOT)+tol and np.abs(partcoords[i,2]) >= i*0.1*max(zcoordsPLOT)-tol):
               #and np.abs(partcoords[i,2]) <= tol):
                radcoordsFIT = np.append(radcoordsFIT, np.sqrt(partcoords[i,1]**2 + partcoords[i,0]**2))
                rhosFIT = np.append(rhosFIT, rhosPLOT[i]) 
                uFIT = np.append(uPLOT, uPLOT[i])
                zcoordsFIT = np.append(zcoordsFIT, zcoordsPLOT)
        breaks = np.append(breaks, len(uFIT))
    cs = np.sqrt(2 * uFIT / 3) 
    vk = 702189.2                   # sqrt(G M / R)
    H = (cs * vk) * radcoordsPLOT ** 1.5
    fit = np.append(fit, max(rhosPLOT) * radcoordsFIT**-3.5 * np.exp(-0.5*(zcoordsPLOT / H)**2))

    
    """ Density Scatter """ 
    
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    ax.set_xlabel('Radial coordinate in Be star Radii')
    ax.set_ylabel('Log of density in kg/m^3')
    ax.set_zlabel('Z coordinate in Be star Radii')
    
    ax.scatter(radcoordsPLOT, np.log10(rhosPLOT), zcoordsPLOT, s=0.75, color=ylims[0], label=ylims[1])
    for i in range(10):
        ax.plot(radcoordsFIT[breaks[i]:breaks[i+1]], np.log10(fit[breaks[i]:breaks[i+1]]), zcoordsFIT[breaks[i]:breaks[i+1]])
    ax.view_init(elev=0, azim=90)
#    plt.xlim(0, 1.5)
#    plt.ylim(ylims[0], ylims[1])
#    plt.xlabel("Distance from centre of Be star in Be Radii")
#    plt.ylabel("log of kg/m^3")
#    plt.legend()
#    plt.title("Plot showing density profile and power law fit (gradient of {:.3})".format(p[0]))
#    plt.title("Plot showing density profile for 24 & 83 day period")
#    plt.show()
    return 
 
   
def TouhamiandGiesfunc(n, rho0, Z, R, u):
    cs = np.sqrt(2 * u / 3) 
    vk = 702189.2                   
    H = (cs * vk) * R ** 1.5
    return  np.exp(-0.5 * (Z / H) ** 2) * rho0 * R ** -n
    
def plotZdependence(partcoords, rho, udist, umass, Rrange, ylims=[-1e-6, 1e-6], tol=1e-3):
    """
    Produces a plot of the density profile for the Be star disk in Z and rho 
    space.

    GEOMETRY
    --------

    z axis is height of disk, y axis is width and x axis is depth.

    INPUT VARIABLES
    ---------------

    partcoords ((npart, 3) float array) - x,y,z coordinates of the disk
                                         particles.

    rhos (npart float array) - density values of the disk particles.

    udist, umass (floats) - units used in simulation
    
    Rrange (array of two floats) - minimum and maximum for 
                                   range of radius from be star to check disk.

    ylim (2 element float) - defines limits of y axis.

    tol (float) - tolerance either side of x=0 that the function will use particles for
          plotting.

    OUTPUT VARIABLES
    ----------------

    A plot of the disk profile.

    """

    " Initialise arrays "

    Rs = np.sqrt(partcoords[:,0]**2 + partcoords[:,1]**2)
    
    """ Conversion of units """
    rho *= umass/(udist**3)
    Rs *= udist / (7 * 695700e3 * 100)
    Z = partcoords[:,2] * udist / (7 * 695700e3 * 100)
    
    """ Density Scatter """
    Zs = np.array([])
    rhos = np.array([])
    for i in range(len(Rs)):
        if(Rs[i] < Rrange[1] and Rs[i] > Rrange[0]):
            Zs = np.append(Zs, Z[i])
            rhos = np.append(rhos, rho[i])
    
    plt.scatter(Zs, rhos, s=0.75, color=ylims[0]) #label=ylims[1] +"gradient = {:.3}".format(p[0]))
    #plt.plot(np.log10(Rs), p[0] * np.log10(Rs) + p[1], color='r')
#    plt.xlim(0, 1.5)
    plt.ylim(min(rhos), max(rhos))
    plt.xlabel("Height of disk particle in Be radii")
    plt.ylabel("log of g/cm^3")
    plt.legend()
    #plt.title("Plot showing density profile and power law fit (gradient of {:.3})".format(p[0]))
    plt.title("Dependence of rho on Z for {}".format(ylims[1]))
    #plt.show()
    return
    
    
def plotPowerLawVariationWithR(partcoords, rhos, u, Rlims, Rsteps, udist, utime, fitType):
    """
    Produces a plot of the variation of the power law index, n, with radius
    included in the calculation.

    GEOMETRY
    --------

    z axis is height of disk, y axis is width and x axis is depth.

    INPUT VARIABLES
    ---------------

    partcoords ((npart, 3) float array) - x,y,z coordinates of the disk
                                         particles.

    rhos (npart float array) - density values of the disk particles.
    
    u (npart float array) - internal energy values of the disk particles.

    udist, utime (floats) - units used in simulation
    
    Rlims (2 float array) - max and min R to perform calculation between 
                            in stellar radii.
    
    Rsteps (int) - number of steps from 1 stellar radii to Rmax.

    fitType (string) - determines what type of fitting algorithm to use 
                       (leastsq or polyfit)

    OUTPUT VARIABLES
    ----------------

    A plot of the disk profile.

    """
    
    if(fitType == "polyfit"):    
        
        Rsteparray = np.linspace(Rlims[0] * (7 * 695700e3 * 100) / udist, Rlims[1] * (7 * 695700e3 * 100) / udist, Rsteps)
        n = np.array([])
        
        Rs = np.sqrt(partcoords[:,0]**2 + partcoords[:,1]**2)   
        
        for i in range(len(Rsteparray)):
            Rtemp = np.array([])
            Ztemp = np.array([])
            utemp = np.array([])
            rhotemp = np.array([])
            for j in range(len(Rs)):
                if(Rs[j] <= Rsteparray[i]):
                    Rtemp = np.append(Rtemp, Rs[j])
                    rhotemp = np.append(rhotemp, rhos[j])
                    
            p = np.polyfit(np.log10(Rtemp * udist), np.log10(rhotemp * umass/(udist**3)), 1)
            
            n = np.append(n, p[0])
        plt.plot(Rsteparray * udist / (7 * 695700e3 * 100), n)
    #    plt.plot(range(12, 26), [4 for i in range(12, 26)])
    #    plt.plot(range(12, 26), [3.5 for i in range(12, 26)])
    #    plt.plot(range(12, 26), [2 for i in range(12, 26)])    
        plt.xlabel("Maximum radius particles are included for fitting")
        plt.ylabel("Power law index n given by fit")

    if(fitType == "leastsq"):    
    
        def TouhamiandGiesfuncLOCAL(n):
            cs = np.sqrt(2 * utemp / 3) 
            vk = 7.021892e7 / (udist / utime)                
            H = (cs * vk) * Rtemp ** 1.5
            return  np.exp(-0.5 * (Ztemp / H) ** 2) * rho0 * Rtemp ** n
        
        Rsteparray = np.linspace(Rlims[0] * (7 * 695700e3 * 100) / udist, Rlims[1] * (7 * 695700e3 * 100) / udist, Rsteps)
        n = np.array([])
        
        Rs = np.sqrt(partcoords[:,0]**2 + partcoords[:,1]**2)
    
        rho0 = max(rhos)    
        
        for i in range(len(Rsteparray)):
            Rtemp = np.array([])
            Ztemp = np.array([])
            utemp = np.array([])
            rhotemp = np.array([])
            for j in range(len(Rs)):
                if(Rs[j] <= Rsteparray[i]):
                    Rtemp = np.append(Rtemp, Rs[j])
                    Ztemp = np.append(Ztemp, partcoords[j,2])
                    utemp = np.append(utemp, u[j])
                    rhotemp = np.append(rhotemp, rhos[j])
            
            n = np.append(n, leastsq(TouhamiandGiesfuncLOCAL, 3.5)[0])
     
        plt.plot(Rsteparray * udist / (7 * 695700e3 * 100), n)
    #    plt.plot(range(12, 26), [4 for i in range(12, 26)])
    #    plt.plot(range(12, 26), [3.5 for i in range(12, 26)])
    #    plt.plot(range(12, 26), [2 for i in range(12, 26)])    
        plt.xlabel("Maximum radius particles are included for fitting")
        plt.ylabel("Power law index n given by fit")
        
def plotPowerLawWindowWithR(filename, Rlims, Rsteps, tol=10.0):
    """
    Produces a plot of the variation of the power law index, n, with radius
    included in the calculation.

    GEOMETRY
    --------

    z axis is height of disk, y axis is width and x axis is depth.

    INPUT VARIABLES
    ---------------

    
    Rlims (2 float array) - max and min R to perform calculation between 
                            in stellar radii.
    
    Rsteps (int) - number of steps from 1 stellar radii to Rmax.



    OUTPUT VARIABLES
    ----------------

    A plot of the disk profile.

    """
    
    diskname = filename[:-26] + "disk001"
    IPname = filename[:-11] + "IPDISK" + filename[-5:-2]
    Becoords, NScoords, partcoords, partvels, rhos, u, num_parts, time, umass, udist, utime = readFILES(filename, diskname, IPname)
  
        
    Rsteparray = np.linspace(Rlims[0] * (7 * 695700e3 * 100) / udist, Rlims[1] * (7 * 695700e3 * 100) / udist, Rsteps)
    n = np.array([])
    
    Rs = np.sqrt(partcoords[:,0]**2 + partcoords[:,1]**2)   
    
    ndiff = np.array([])
    for i in range(1, len(Rsteparray)):
        Rtemp = np.array([])
        Ztemp = np.array([])
        utemp = np.array([])
        rhotemp = np.array([])
        for j in range(len(Rs)):
            if(Rs[j] <= Rsteparray[i]):
                Rtemp = np.append(Rtemp, Rs[j])
                rhotemp = np.append(rhotemp, rhos[j])
        if(len(Rtemp) < 2):
            break
        p = np.polyfit(np.log10(Rtemp * udist), np.log10(rhotemp * umass/(udist**3)), 1)
        
        n = np.append(n, p[0])
#        print n[-1]
#        print Rsteparray[i] * udist / (7 * 695700e3 * 100)
        if(len(n) > 2):
            if(np.abs(n[-1]) / np.abs(n[-2]) > tol):
#                print Rtemp[-1] * udist / (7 * 695700e3 * 100)
                print "YES"
#            deltaY = np.log10(rhotemp * umass/(udist**3))[-1] - np.log10(rhotemp * umass/(udist**3))[-2]
#            deltaX = np.log10(Rtemp * udist)[-1] - np.log10(Rtemp * udist)[-2]
#            n = np.append(n, deltaY / deltaX)
    NDIFF = 0
    for i in range(2, len(n)):
        ndiff = np.append(ndiff, n[i] - n[i-2])
    for i in range(2, len(ndiff)):
        if(ndiff[i-2] > ndiff[i-1] and ndiff[i-1] < ndiff[i] and ndiff[i] > -0.5 and ndiff[i] < NDIFF):
#            print ndiff[i]
            NDIFF = ndiff[i]
            NDIFFR = Rsteparray[i] * udist / (7 * 695700e3 * 100)
#            print ndiff[i]
    print filename[18:-26], NDIFFR

#    print Rsteparray[np.argmin(ndiff)] * udist / (7 * 695700e3 * 100)
#    print Rsteparray[np.argmax(ndiff)] * udist / (7 * 695700e3 * 100)
    plt.plot(Rsteparray[:len(ndiff)] * udist / (7 * 695700e3 * 100), ndiff)
#    plt.plot(range(12, 26), [4 for i in range(12, 26)])
#    plt.plot(range(12, 26), [3.5 for i in range(12, 26)])
#    plt.plot(range(12, 26), [2 for i in range(12, 26)])    
    plt.xlabel("Radius in stellar radii")
    plt.ylabel("Derivative of power law index n")


def plotDiskFlare(partcoords, udist, Rsteps):
    """
    Produces a plot of the disk flare in the z direction with radius
    included in the calculation.

    GEOMETRY
    --------

    z axis is height of disk, y axis is width and x axis is depth.

    INPUT VARIABLES
    ---------------

    partcoords ((npart, 3) float array) - x,y,z coordinates of the disk
                                         particles.

    udist (float) - units used in simulation
    
    Rsteps (int) - number of steps from 1 stellar radii to Rmax.

    tol (float) - tolerance either side of x=0 that the function will use particles for
          plotting.

    OUTPUT VARIABLES
    ----------------

    A plot of the disk flaring.

    """
    
    R = np.sqrt(partcoords[:,0]**2 + partcoords[:,1]**2)
    tol = 0.5 * max(R) / Rsteps
    Rsteparray = np.linspace((7 * 695700e3 * 100) / udist, max(R), Rsteps)
    Rs = np.array([])
    Zs = np.array([])
    for j in range(len(Rsteparray)):
        Rtemp = np.array([])
        Ztemp = np.array([])
        for i in range(len(partcoords)):
#            print np.abs(Rsteparray[j] - R[i]) < tol
            if(np.abs(Rsteparray[j] - R[i]) < tol):
                Rtemp = np.append(Rtemp, R[i])
                Ztemp = np.append(Ztemp, partcoords[i, 2])
        Rs = np.append(Rs, np.sum(Rtemp) / len(Rtemp))
        Zs = np.append(Zs, np.sum(Ztemp) / len(Ztemp))
        
#    plt.scatter(Rs * udist / (7 * 695700e3 * 100), Zs * udist / (7 * 695700e3 * 100))   
    plt.scatter(R * udist / (7 * 695700e3 * 100), partcoords[:,2] * udist / (7 * 695700e3 * 100), s=0.5) 
    plt.xlabel("Radius from Be star in stellar radii")
    plt.ylabel("Maximum height of particle at given radius") 
 
   
def plotFlareLineWithRadius(partcoords, udist,  Rsteps):
    """
    Produces a plot of the density profile for the Be star disk in Z and rho 
    space.

    GEOMETRY
    --------

    z axis is height of disk, y axis is width and x axis is depth.

    INPUT VARIABLES
    ---------------

    partcoords ((npart, 3) float array) - x,y,z coordinates of the disk
                                         particles.

    rhos (npart float array) - density values of the disk particles.

    udist, umass (floats) - units used in simulation
    
    Rsteps (int) - number of steps from minimum to maximum radius of be star.

    ylim (2 element float) - defines limits of y axis.

    tol (float) - tolerance either side of x=0 that the function will use particles for
          plotting.

    OUTPUT VARIABLES
    ----------------

    A plot of the disk profile.

    """

    " Initialise arrays "

    Rs = np.sqrt(partcoords[:,0]**2 + partcoords[:,1]**2)
    
    """ Conversion of units """
    Rs *= udist / (7 * 695700e3 * 100)
    Z = partcoords[:,2] * udist / (7 * 695700e3 * 100)
    
    R = np.array([])
    EW = np.array([])
    EWneg = np.array([])
    Zmaxs = np.array([])
    for j in np.delete(np.linspace(min(Rs), max(Rs), Rsteps), -1):
        j2 = j + ( max(Rs) - min(Rs) ) / Rsteps
        R = np.append(R, j)
        Zs = np.array([])
        for i in range(len(Rs)):
            if(Rs[i] >= j and Rs[i] < j2):
                Zs = np.append(Zs, Z[i])
        EW = np.append(EW, max(Zs))
        EWneg = np.append(EWneg, min(Zs))
        print max(Zs)
        Zmaxs = np.append(Zmaxs, max(Zs))

    plt.plot(R, EW, color='b')
    plt.plot(R, EWneg, color='b')
    plt.xlabel("Radius of disk particle in Be radii")
    plt.ylabel("Height of disk particle in Be radii")
    #plt.show()
    return Zmaxs

def plotDensityProfileFromPartNumbers(filenames, Rlims, Rsteps, pltargs=[0], colour='b'):
    """
    Produces a plot of the density profile for the Be star disk in Z and rho 
    space.

    GEOMETRY
    --------

    z axis is height of disk, y axis is width and x axis is depth.

    INPUT VARIABLES
    ---------------

    partcoords ((npart, 3) float array) - x,y,z coordinates of the disk
                                         particles.

    rhos (npart float array) - density values of the disk particles.

    udist, umass (floats) - units used in simulation
    
    Rsteps (int) - number of steps from minimum to maximum radius of be star.

    ylim (2 element float) - defines limits of y axis.

    tol (float) - tolerance either side of x=0 that the function will use particles for
          plotting.

    OUTPUT VARIABLES
    ----------------

    A plot of the disk profile.

    """
    
    for filename in filenames:
        diskname = filename[:-26] + "disk001"
        IPname = filename[:-11] + "IPDISK" + filename[-5:-2]
#        diskname = '../../Desktop/40d_0.0e_eqTest/disk001'
#        IPname = '../../Desktop/40d_0.0e_eqTest/IPDISK' + filename[-5:-2]
        Becoords, NScoords, partcoords, partvels, rhos, u, num_parts, time, umass, udist, utime = readFILES(filename, diskname, IPname)

        " Initialise arrays "
    
        Rs = np.sqrt(partcoords[:,0]**2 + partcoords[:,1]**2)
        
        """ Conversion of units """
        Rs *= udist / (7 * 695700e3 * 100)
        Z = partcoords[:,2] * udist / (7 * 695700e3 * 100)
        
        R = np.array([])
        numpartsDens = np.array([])
        for j in np.delete(np.linspace(Rlims[0], Rlims[1], Rsteps), -1):
            j2 = j + ( max(Rs) - min(Rs) ) / Rsteps
            R = np.append(R, j + 0.5 * (j2 - j))
            num = 0
            Zs = np.array([])
            for i in range(len(Rs)):
                if(Rs[i] >= j and Rs[i] < j2):
                    Zs = np.append(Zs, Z[i])
                    num += 1
            if(num == 0):
                numpartsDens = np.append(numpartsDens, 0)
            else:
                volume = 2 * np.pi * (j2**2 - j**2) * max(Zs)
                numpartsDens = np.append(numpartsDens, num / volume)
        
        numpartsDens *= max(rhos) * umass/(udist**3) / max(numpartsDens)
        if(pltargs[0] == 0):
            plt.plot(R, np.log10(numpartsDens), label='Density per volume')#, color=colour)
        else:
            plt.plot(np.log10(R), np.log10(numpartsDens), label=pltargs[0]+' Density per volume')#, color=colour)
        plt.xlabel("Radius of disk particle in Be radii")
        plt.ylabel("Density in g per cube cm")
        #plt.xlim(0, 25)
    #    plt.title("Density and Number counted density")
    #    plt.legend()
    return

def plotDensityProfileFromPartNumbersWithFit(partcoords, rhos, udist, umass, Rlims, Rsteps, pltargs=[0], colour='b'):
    """
    Produces a plot of the density profile for the Be star disk in Z and rho 
    space.

    GEOMETRY
    --------

    z axis is height of disk, y axis is width and x axis is depth.

    INPUT VARIABLES
    ---------------

    partcoords ((npart, 3) float array) - x,y,z coordinates of the disk
                                         particles.

    rhos (npart float array) - density values of the disk particles.

    udist, umass (floats) - units used in simulation
    
    Rsteps (int) - number of steps from minimum to maximum radius of be star.

    ylim (2 element float) - defines limits of y axis.

    tol (float) - tolerance either side of x=0 that the function will use particles for
          plotting.

    OUTPUT VARIABLES
    ----------------

    A plot of the disk profile.

    """

    " Initialise arrays "

    Rs = np.sqrt(partcoords[:,0]**2 + partcoords[:,1]**2)
    
    """ Conversion of units """
    Rs *= udist / (7 * 695700e3 * 100)
    Z = partcoords[:,2] * udist / (7 * 695700e3 * 100)
    
    R = np.array([])
    numpartsDens = np.array([])
    for j in np.delete(np.linspace(Rlims[0], Rlims[1], Rsteps), -1):
        j2 = j + ( max(Rs) - min(Rs) ) / Rsteps
        R = np.append(R, j + 0.5 * (j2 - j))
        num = 0
        Zs = np.array([])
        for i in range(len(Rs)):
            if(Rs[i] >= j and Rs[i] < j2):
                Zs = np.append(Zs, Z[i])
                num += 1
        if(num == 0):
            numpartsDens = np.append(numpartsDens, 0)
        else:
            volume = 2 * np.pi * (j2**2 - j**2) * max(Zs)
            numpartsDens = np.append(numpartsDens, num / volume)
    
    numpartsDens *= max(rhos) * umass/(udist**3) / max(numpartsDens)
    
    "polyfit"
#    " Create arrays for fit "
#    R = np.sqrt(partcoords[:,0]**2 + partcoords[:,1]**2)
#    delArr = np.array([])
#    for i in range(len(T)):
#        if(R[i] > Rlims[1]):
#            delArr = np.append(delArr, i)
#            
#    xcoords = np.delete(xcoords, delArr)
#    ycoords = np.delete(ycoords, delArr)
#    zcoords = np.delete(zcoords, delArr)
#    T = np.delete(T, delArr)
    p = np.polyfit(np.log10(R), np.log10(numpartsDens), 1)
    print "polyfit:", p[0], p[1]
#    print max(R), max(np.log10(R))
    
    if(pltargs[0] == 0):
        plt.plot(np.log10(R), np.log10(numpartsDens), label='Density per volume', color=colour)
        plt.plot(np.log10(R), np.log10(R)*(-2.98) - 10.79)
    else:
        plt.plot(np.log10(R), np.log10(numpartsDens), label=pltargs[0]+' Density per volume', color=colour)
    plt.xlabel("log Radius of disk particle in Be radii")
    plt.ylabel("Density in g per cube cm")
    plt.xlim(0, 1.5)
    plt.ylim(-16.5, -11)
#    plt.title("Density and Number counted density")
#    plt.legend()
    return

def plotColumnDensityProfileFromPartNumbers(partcoords, rhos, udist, umass, Rlims, Rsteps):
    """
    Produces a plot of the density profile for the Be star disk in Z and rho 
    space.

    GEOMETRY
    --------

    z axis is height of disk, y axis is width and x axis is depth.

    INPUT VARIABLES
    ---------------

    partcoords ((npart, 3) float array) - x,y,z coordinates of the disk
                                         particles.

    rhos (npart float array) - density values of the disk particles.

    udist, umass (floats) - units used in simulation
    
    Rsteps (int) - number of steps from minimum to maximum radius of be star.

    ylim (2 element float) - defines limits of y axis.

    tol (float) - tolerance either side of x=0 that the function will use particles for
          plotting.

    OUTPUT VARIABLES
    ----------------

    A plot of the disk profile.

    """

    " Initialise arrays "

    Rs = np.sqrt(partcoords[:,0]**2 + partcoords[:,1]**2)
    
    """ Conversion of units """
    Rs *= udist / (7 * 695700e3 * 100)
    Z = partcoords[:,2] * udist / (7 * 695700e3 * 100)
    
    R = np.array([])
    numpartsDens = np.array([])
    for j in np.delete(np.linspace(Rlims[0], Rlims[1], Rsteps), -1):
        j2 = j + ( max(Rs) - min(Rs) ) / Rsteps
        R = np.append(R, j + 0.5 * (j2 - j))
        num = 0
        for i in range(len(Rs)):
            if(Rs[i] >= j and Rs[i] < j2):
                num += 1

        area = np.pi * (j2**2 - j**2) 
        numpartsDens = np.append(numpartsDens, num / area)
    
    numpartsDens *= max(rhos) * umass/(udist**3) / max(numpartsDens)
    p = np.polyfit(np.log10(R), np.log10(numpartsDens), 1)
    print p
    plt.plot(np.log10(R), np.log10(numpartsDens), color='g', label='Column density')
    plt.xlabel("Radius of disk particle in Be radii")
    plt.ylabel("Density in g per cube cm")
    plt.title("Density and Number counted density")
    plt.legend()
#    plt.xlim(1, 10)
    return 

"""----------------------------------------------------------------------------

                    2D HIST

----------------------------------------------------------------------------"""  

def plotDensity2dHist(filename, phi, theta, show='y', save='n', outFilename="figure", cbar=[0,60],
                 xylims=[-25,25,-25,25]):
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
 
    diskname = filename[:-26] + "disk001"
    IPname = filename[:-11] + "IPDISK" + filename[-5:-2]
    Becoords, NScoords, partcoords, partvels, rho, u, num_parts, time, umass, udist, utime = readFILES(filename, diskname, IPname)
    Becoords *= udist / (7 * 6.96e10)
    NScoords *= udist / (7 * 6.96e10)
    partcoords *= udist / (7 * 6.96e10)

    BINS = 200
    h, xedges, yedges = np.histogram2d(partcoords[:,0], partcoords[:,1], bins=BINS)
#    print xedges - xedges2
#    print yedges - yedges2
    
    plt.pcolor(xedges, yedges, h, cmap='RdBu')
    plt.xlim(xylims[0], xylims[1])
    plt.ylim(xylims[2], xylims[3])
#    plt.clim(cbar[0], cbar[1])
#    plt.ylabel("EW")
#    plt.xlabel("Time / days")
    plt.colorbar()
    
    if(show == 'y'):
        plt.show()
    
    if(save == 'y'):
        plt.savefig(outFilename)#, figsize=(18, 16), dpi=100)

    if(show == 'n'):
        plt.close()
        
"""----------------------------------------------------------------------------

                    RING DENSITY

----------------------------------------------------------------------------"""  


def plotRingDensity(filename, colour='b'):
    """
    Produces a plot of the density profile for the Be star disk in Z and rho 
    space.

    GEOMETRY
    --------

    z axis is height of disk, y axis is width and x axis is depth.

    INPUT VARIABLES
    ---------------

    partcoords ((npart, 3) float array) - x,y,z coordinates of the disk
                                         particles.

    rhos (npart float array) - density values of the disk particles.

    udist, umass (floats) - units used in simulation
    
    Rsteps (int) - number of steps from minimum to maximum radius of be star.

    ylim (2 element float) - defines limits of y axis.

    tol (float) - tolerance either side of x=0 that the function will use particles for
          plotting.

    OUTPUT VARIABLES
    ----------------

    A plot of the disk profile.

    """
    
    diskname = filename[:-26] + "disk001"
    IPname = filename[:-11] + "IPDISK" + filename[-5:-2]
    Becoords, NScoords, partcoords, partvels, rho, u, num_parts, time, umass, udist, utime = readFILES(filename, diskname, IPname)
    Becoords *= udist / (7 * 6.96e10)
    NScoords *= udist / (7 * 6.96e10)
    partcoords *= udist / (7 * 6.96e10)

    " Initialise arrays "
    R = np.sqrt(partcoords[:,0]**2 + partcoords[:,1]**2)
    Rs = np.array([])
    phi = np.array([])
    x = np.array([])
    y = np.array([])
    rhos = np.array([])
    
    Rmax = 1.05
    Rmin = 1.0
    zFac = 0.03
    for i in range(len(R)):
        if(R[i] > Rmin and R[i] < Rmax 
           and partcoords[i,2] > -zFac*(Rmax-Rmin) 
           and partcoords[i,2] < zFac*(Rmax-Rmin)):
            Rs = np.append(Rs, R[i])
            rhos = np.append(rhos, rho[i])
            phi = np.append(phi, np.arctan2(partcoords[i,0], partcoords[i,1]))
            x = np.append(x, partcoords[i,0])
            y = np.append(y, partcoords[i,1])
    
    """ Conversion of units """
    Z = partcoords[:,2]
    rhos *= umass/(udist**3)
    
    sortOrder = np.argsort(phi)
    phi = phi[sortOrder]
    rhos = rhos[sortOrder]
    print phi
    
    """ SINE CURVE FIT """
    
#    N = 1000 # number of data points
#    t = np.linspace(0, 4*np.pi, N)
#    data = 3.0*np.sin(t+0.001) + 0.5 + np.random.randn(N) # create artificial data with noise
    N = len(rhos)
    data = rhos
    t = phi
    
    guess_mean = np.mean(data)
    guess_std = 3*np.std(data)/(2**0.5)
    guess_phase = 0
    
    # we'll use this to plot our first estimate. This might already be good enough for you
    data_first_guess = guess_std*np.sin(t+guess_phase) + guess_mean
    
    # Define the function to optimize, in this case, we want to minimize the difference
    # between the actual data and our "guessed" parameters
    optimize_func = lambda x: x[0]*np.sin(t+x[1]) + x[2] - data
    est_std, est_phase, est_mean = leastsq(optimize_func, [guess_std, guess_phase, guess_mean])[0]
    
    # recreate the fitted curve using the optimized parameters
    data_fit = est_std*np.sin(t+est_phase) + est_mean
    
    plt.plot(data, '.')
    plt.plot(data_fit, label='after fitting')
    plt.plot(data_first_guess, label='first guess')
    plt.legend()
    plt.show()
#    print FFIT

#    plt.scatter(phi, rhos)
#    plt.plot(phi, FFIT)    
#    plt.plot(x,y)
#    plt.xlabel("Radius of disk particle in Be radii")
#    plt.ylabel("Density in g per cube cm")
#    plt.ylim(1e-10, 1e-14)
#    plt.title("Density and Number counted density")
#    plt.legend()
    return


"""----------------------------------------------------------------------------

                    FOURIER SERIES

----------------------------------------------------------------------------""" 

time = np.linspace(0, 8*np.pi, 100)
y = np.array([])
for i in range(len(time)):
    if(i % 5 == 0):
        y = np.append(y, 0)
    else:
        y = np.append(y, 1)

def cn(y, n, time):
    period = 2*np.pi
    c = y*np.exp(-1j*2*n*np.pi*time/period)
    return c.sum()/c.size

def f(x, y, Nh, time):
    period = 2*np.pi
    f = np.array([2*cn(y, i, time)*np.exp(1j*2*i*np.pi*x/period) for i in range(1,Nh+1)])
    return f.sum()

y2 = np.array([f(t,y,2,time).real for t in time])

#plt.plot(time, y)
#plt.plot(time, y2)

    
"""----------------------------------------------------------------------------

                    RUNNING STUFF

----------------------------------------------------------------------------"""  

#filename = '/data/rob1g10/SPH/DATA/DISKS/injection_off_a1.5e11/processedDISKs/QGDISK14701'
#filename = '/data/rob1g10/SPH/DATA/DISKS/20d_0.4e/processedDISKs/QGDISK28001'
#filename = '/data/rob1g10/SPH/DATA/DISKS/20d_0.4e/processedDISKs/QGDISK28201'
#filename = '/data/rob1g10/SPH/DATA/DISKS/justBe_halfInj2/processedDISKs/QGDISK50301'
#filename = '/data/rob1g10/SPH/DATA/DISKS/injection_off_justBe/processedDISKs/QGDISK17801'
#filename = '/data/rob1g10/SPH/DATA/DISKS/40d_0.4e/processedDISKs/QGDISK13401'
#filename = '/data/rob1g10/SPH/DATA/DISKS/80d_0.0e/processedDISKs/QGDISK06101'
#filename = '/data/rob1g10/SPH/injFiddle/sph3D/processedDISKs/QGDISK32701'
#filename = '/data/rob1g10/SPH/DATA/DISKS/OkazakiRecreate/processedDISKs/QGDISK08301'
#filename = '/data/rob1g10/SPH/DATA/DISKS/20d_0.4e/processedDISKs/QGDISK28001'
#filename = '/data/rob1g10/SPH/DATA/DISKS/sph3D_f90_radialTemp_q0.5/processedDISKs/QGDISK08101'
#filename = '/data/rob1g10/SPH/DATA/DISKS/isoBe_injFixed/processedDISKs/QGDISK04901'
#filename = '/data/rob1g10/SPH/DATA/DISKS/isoBe_fromBeginningTenth/processedDISKs/QGDISK41101'
#filename = '/data/rob1g10/SPH/DATA/DISKS/20d_0.0e_eqTest/processedDISKs/QGDISK38701'

# filenames = ['../../Documents/PhD/Okazaki SPH/RadialTempSPH/images/QGDISK01102',
#             '../../Documents/PhD/Okazaki SPH/RadialTempSPH/images/QGDISK01104',
#              '../../Documents/PhD/Okazaki SPH/RadialTempSPH/images/QGDISK01106'
#              ]
# i=0
# names = ['02', '04', '06']
# colours = ['b','g', 'r']
# for filename in filenames:
#     """ Read File """
    
#     try:
#         Becoords, NScoords, partcoords, partvels, rho, u, num_parts, time, umass, udist, utime = readFILES(filename, diskname)
#     except:
#         Becoords, NScoords, partcoords, partvels, rho, u, num_parts, umass, udist, utime = readFILES(filename, diskname)
#         time = 1
# #    name = filename[29:37]
#     plotDensityProfileFromPartNumbers(partcoords, rho, udist, umass, [1, 2.5], 200, colour=colours[i])
#     i+=1

""" Read File """
#diskname = filename[:-26] + "disk001"
#IPname = filename[:-11] + "IPDISK" + filename[-5:-2]
#filename = 'DATA/QGDISK41101'
#diskname = 'INPUT/disk001'
#IPname = 'INPUT/IPDISK411'
#Becoords, NScoords, partcoords, partvels, rho, u, num_parts, time, umass, udist, utime = readFILES(filename, diskname, IPname)
#print num_parts

"""Touhami and Gies"""

#R = np.sqrt(partcoords[:,0]**2 + partcoords[:,1]**2)
#R *= udist / (7 * 6.96e10)

Rtemp = np.array([])
Ztemp = np.array([])
utemp = np.array([])
rhotemp = np.array([])
#for j in range(len(R)):
#    if(R[j] <= 1.5):
#        Rtemp = np.append(Rtemp, R[j])
#        rhotemp = np.append(rhotemp, rho[i])
        
#p = np.polyfit(np.log10(Rtemp), np.log10(rhotemp), 1)
#pALL = np.polyfit(np.log10(R), np.log10(rhos), 1)
#plt.scatter(Rtemp, rhotemp)
#plt.show()
#print np.polyfit(np.log10(R), np.log10(rho), 1)[0]
#Z = partcoords[:,2]
#rho0 = max(rho)
#
#Rtemp = np.array([])
#Ztemp = np.array([])
#utemp = np.array([])
#for i in range(len(R)):
#    if(np.abs(R[i]) <= 1.17):
#        Rtemp = np.append(Rtemp, R[i])
#        Ztemp = np.append(Ztemp, Z[i])
#        utemp = np.append(utemp, u[i])
#        
#def TouhamiandGiesfunc2(n):
#    cs = np.sqrt(2 * utemp / 3) 
#    vk = 7.021892e7 / (udist / utime)                
#    H = (cs * vk) * Rtemp ** 1.5
#    return  np.exp(-0.5 * (Ztemp / H) ** 2) * rho0 * Rtemp ** -n
#    
#n = leastsq(TouhamiandGiesfunc2, 3.5)
#print n
    
#fig = plt.figure()
#ax = fig.add_subplot(111, projection='3d')
#
#ax.set_xlabel('Radial coordinate in Be star Radii')
#ax.set_ylabel('Log of density in kg/m^3')
#ax.set_zlabel('Z coordinate in Be star Radii')
#
#ax.scatter(R, np.log10(rho), Z, s=0.75, color='b', label='data')
#ax.scatter(R, np.exp(-0.5 * (Ztemp / H) ** 2) * rho0 * Rtemp ** -, Z)
#ax.view_init(elev=0, azim=90)
#udist = 5.7700e12
#utime = 2.8344e05
#umass = 3.583e34   
    
"""Rivinius and Carciofi"""

#plotDensityProfileRiviniusandCarciofi(partcoords, rho, udist, umass, ylims=['b', 'justBe'])

"""Normal"""

#rho0 = max(rho)
#
#xtemp = np.array([])
#ytemp = np.array([])
#ztemp = np.array([])
#rhotemp = np.array([])
#utemp = np.array([])
#for i in range(len(partcoords)):
##if(rho[i] > rho0 * 1.22e-4):
#    xtemp = np.append(xtemp, partcoords[i,0])
#    ytemp = np.append(ytemp, partcoords[i,1])
#    ztemp = np.append(ztemp, partcoords[i,2])
#    rhotemp =  np.append(rhotemp, rho[i])
#    utemp = np.append(utemp, u[i])
#
#partcoordstemp = np.zeros((len(xtemp), 3))
#for i in range(len(xtemp)):
#    partcoordstemp[i] = [xtemp[i], ytemp[i], ztemp[i]]
#  
#plotDensityProfile(['/data/rob1g10/SPH/DATA/DISKS/M1.4_20d_0.0e/processedDISKs/QGDISK19315',
#                    '/data/rob1g10/SPH/DATA/DISKS/M6_17.4d_0.0e/processedDISKs/QGDISK39533',
#                    '/data/rob1g10/SPH/DATA/DISKS/M10_15.8d_0.0e/processedDISKs/QGDISK32745',],
#                    ylims=['b', 'just Be '])
                
#plotDensityProfile(['/data/rob1g10/SPH/DATA/DISKS/isoBe_M15_Mdot9/processedDISKs/QGDISK08701',
#                    '/data/rob1g10/SPH/DATA/DISKS/isoBe_M15_Mdot10/processedDISKs/QGDISK08101'],
#                    ylims=['b', 'just Be '])

#plotDensityProfile(['/data/rob1g10/SPH/DATA/DISKS/M1.4_43d_0.4e/processedDISKs/QGDISK51901',
#                    '/data/rob1g10/SPH/DATA/DISKS/M6_37.5d_0.4e/processedDISKs/QGDISK16709',
#                    '/data/rob1g10/SPH/DATA/DISKS/M10_34.1d_0.4e/processedDISKs/QGDISK15001'],
#                    ylims=['b', 'just Be '])
                    
#plotDensityProfile(['/data/rob1g10/SPH/DATA/DISKS/M14_20d_0.0e/processedDISKs/QGDISK35901',
#                   '/data/rob1g10/SPH/DATA/DISKS/M10_20d_0.0e/processedDISKs/QGDISK30001',
#                   '/data/rob1g10/SPH/DATA/DISKS/M8_20d_0.0e/processedDISKs/QGDISK36001',
#                   '/data/rob1g10/SPH/DATA/DISKS/M6_20d_0.0e/processedDISKs/QGDISK36001',
#                   '/data/rob1g10/SPH/DATA/DISKS/M1.4_20d_0.0e/processedDISKs/QGDISK36001'],
#                    ylims=['b', 'just Be '])

#plotDensityProfile(['/data/rob1g10/SPH/DATA/DISKS/20d_0.2e_eqTest/processedDISKs/QGDISK13801',
#                    '/data/rob1g10/SPH/DATA/DISKS/BH20d_0.292e/processedDISKs/QGDISK13801',
#                    '/data/rob1g10/SPH/DATA/DISKS/BH16.7d_0.2e/processedDISKs/QGDISK36101',
#                    '/data/rob1g10/SPH/DATA/DISKS/BH14.3d_0.2e/processedDISKs/QGDISK36101',
#                    '/data/rob1g10/SPH/DATA/DISKS/BH24d_0.2e/processedDISKs/QGDISK43101',
#                    '/data/rob1g10/SPH/DATA/DISKS/BH28d_0.2e/processedDISKs/QGDISK43101'], ylims=['b', 'just Be '])

#plotDensityProfile(['/data/rob1g10/SPH/DATA/DISKS/isoBe_fromBeginningTenth/processedDISKs/QGDISK41101'], ylims=['b', 'just Be '])
#plotDensityProfile(['DATA/QGDISK41101'], ylims=['b', 'just Be '])
#plotDensityProfile(['../../Desktop/40d_0.0e_eqTest/QGDISK38108'], ylims=['b', 'just Be '])
plotEquatorialDensityProfile(['/data/rob1g10/SPH/DATA/DISKS/M1.4_20d_0.0e/processedDISKs/QGDISK19315',
                         '/data/rob1g10/SPH/DATA/DISKS/M1.4_28d_0.2e/processedDISKs/QGDISK30918',
                         '/data/rob1g10/SPH/DATA/DISKS/M1.4_43d_0.4e/processedDISKs/QGDISK53903',
                         '/data/rob1g10/SPH/DATA/DISKS/M1.4_79.1d_0.6e/processedDISKs/QGDISK38204',
                         
                         '/data/rob1g10/SPH/DATA/DISKS/M6_17.4d_0.0e/processedDISKs/QGDISK39533',
                         '/data/rob1g10/SPH/DATA/DISKS/M6_24.3d_0.2e/processedDISKs/QGDISK36042',
                         '/data/rob1g10/SPH/DATA/DISKS/M6_37.5d_0.4e/processedDISKs/QGDISK15836',
                         '/data/rob1g10/SPH/DATA/DISKS/M6_68.8d_0.6e/processedDISKs/QGDISK27501',
                         
                         '/data/rob1g10/SPH/DATA/DISKS/M10_15.8d_0.0e/processedDISKs/QGDISK32745',
                         '/data/rob1g10/SPH/DATA/DISKS/M10_22.1d_0.2e/processedDISKs/QGDISK37237',
                         '/data/rob1g10/SPH/DATA/DISKS/M10_34.1d_0.4e/processedDISKs/QGDISK14118',
                         '/data/rob1g10/SPH/DATA/DISKS/M10_62.6d_0.6e/processedDISKs/QGDISK13413',
                         ], 
                         ylims=['b', 'just Be '])
#Rs, rhos = plotEquatorialDensityProfile(['DATA/QGDISK41101'], ylims=['b', 'just Be '], tol=1e-3)
#x = np.log10(Rs)
#y = np.log10(rhos)
#p1 = np.polyfit(x, y, 1)
#p2 = np.polyfit(x, y, 2)
#p3 = np.polyfit(x, y, 3)
#p4 = np.polyfit(x, y, 4)
#pfit = [-0.09690832,  -0.45785672,  -0.53858783,  -2.14967644, -10.75693424]
xfit = np.linspace(1, 25, 100)
yfit = np.log10(np.linspace(1, 25, 100))
#plt.plot(xfit, 10**(pfit[0]*yfit**4 + pfit[1]*yfit**3 + pfit[2]*yfit**2 + pfit[3]*yfit + pfit[4]), color="black")
#plt.scatter(10**x, 10**(p1[0]*x + p1[1]), color="black", s=0.5)
#plt.scatter(10**x, 10**(p2[0]*x**2 + p2[1]*x + p2[2]), color="green", s=0.5)
#plt.scatter(10**x, 10**(p3[0]*x**3 + p3[1]*x**2 + p3[2]*x + p3[3]), color="red", s=0.5)
#plt.scatter(x, (p4[0]*x**4 + p4[1]*x**3 + p4[2]*x**2 + p4[3]*x + p4[4]), color="yellow", s=0.5)
#plt.scatter(x, 2e-11*x**-2.8, color="yellow", s=0.5)
#plt.ylim(-2e-12, 2e-11)
#print p1
#print p4
#print Rs
   
#R = np.sqrt(partcoordstemp[:,0]**2 + partcoordstemp[:,1]**2)
#Z = partcoordstemp[:,2]
#
#Rtemp = np.array([])
#Ztemp = np.array([])
#utemp = np.array([])
#for i in range(len(R)):
#    if(np.abs(R[i]) <= 1.17):
#        Rtemp = np.append(Rtemp, R[i])
#        Ztemp = np.append(Ztemp, Z[i])
#        utemp = np.append(utemp, u[i])
#
#n = leastsq(TouhamiandGiesfunc2, 3.5)
#print n

#plt.scatter(np.log10(TouhamiandGiesfunc2(-3.5)), np.log10(rhotemp))
#plt.xlim(0, 6e-8)
#plt.ylim(0, 7e-7)

""" Column Density """

#plotColumnDensityProfileFromPartNumbers(partcoords, rho, udist, umass, [1, 10], 200)

""" Counting number of parts for comparison to SPH density """
#plotDensityProfileFromPartNumbers(['../../Desktop/40d_0.0e_eqTest/QGDISK38001'], [1, 25], 20)
#plotDensityProfileFromPartNumbers(['/data/rob1g10/SPH/DATA/DISKS/isoBe_fromBeginningTenth/processedDISKs/QGDISK50101',
#                                   '/data/rob1g10/SPH/DATA/DISKS/isoBe_fromBeginningTenth/processedDISKs/QGDISK39101',
#                                   '/data/rob1g10/SPH/DATA/DISKS/isoBe_fromBeginningTenth/processedDISKs/QGDISK39501'],
#                                   [1, 25], 200)
#plotDensityProfileFromPartNumbers(['/data/rob1g10/SPH/DATA/DISKS/20d_0.4e_eqTest/processedDISKs/QGDISK33020',
#                                   '/data/rob1g10/SPH/DATA/DISKS/20d_0.4e_eqTest/processedDISKs/QGDISK33021',
#                                   '/data/rob1g10/SPH/DATA/DISKS/20d_0.4e_eqTest/processedDISKs/QGDISK33022',
#                                   '/data/rob1g10/SPH/DATA/DISKS/20d_0.4e_eqTest/processedDISKs/QGDISK33023'],
#                                   [1, 15], 200)
#plotDensityProfileFromPartNumbers(['/data/rob1g10/SPH/DATA/DISKS/20d_0.2e_eqTest/processedDISKs/QGDISK13801',
#                    '/data/rob1g10/SPH/DATA/DISKS/BH20d_0.292e/processedDISKs/QGDISK13801',
#                    '/data/rob1g10/SPH/DATA/DISKS/BH16.7d_0.2e/processedDISKs/QGDISK36101',
#                    '/data/rob1g10/SPH/DATA/DISKS/BH14.3d_0.2e/processedDISKs/QGDISK36101'],
#                                   [1, 15], 200)
#plotDensityProfileFromPartNumbers(['/data/rob1g10/SPH/DATA/DISKS/20d_0.0e_eqTest/processedDISKs/QGDISK38701',
#                                   '/data/rob1g10/SPH/DATA/DISKS/20d_0.2e_eqTest/processedDISKs/QGDISK38701',
#                                   '/data/rob1g10/SPH/DATA/DISKS/20d_0.4e_eqTest/processedDISKs/QGDISK38701'],
#                                   [1, 25], 115)

#plotDensityProfileFromPartNumbers([#'/data/rob1g10/SPH/DATA/DISKS/M1.4_20d_0.0e/processedDISKs/QGDISK36403',
#            #'/data/rob1g10/SPH/DATA/DISKS/M6_17.4d_0.0e/processedDISKs/QGDISK12513',
#            #'/data/rob1g10/SPH/DATA/DISKS/M10_15.8d_0.0e/processedDISKs/QGDISK32014'],
#            
##            '/data/rob1g10/SPH/DATA/DISKS/M1.4_28d_0.2e/processedDISKs/QGDISK30019',
##            '/data/rob1g10/SPH/DATA/DISKS/M6_24.3d_0.2e/processedDISKs/QGDISK35017',
##            '/data/rob1g10/SPH/DATA/DISKS/M10_22.1d_0.2e/processedDISKs/QGDISK36017'],
#            '/data/rob1g10/SPH/DATA/DISKS/M1.4_20d_0.0e/processedDISKs/QGDISK44801',
#            '/data/rob1g10/SPH/DATA/DISKS/M6_17.4d_0.0e/processedDISKs/QGDISK40401',
#            '/data/rob1g10/SPH/DATA/DISKS/M10_15.8d_0.0e/processedDISKs/QGDISK33601'],
#            
##            '/data/rob1g10/SPH/DATA/DISKS/M1.4_43d_0.4e/processedDISKs/QGDISK44801',
##            '/data/rob1g10/SPH/DATA/DISKS/M6_37.5d_0.4e/processedDISKs/QGDISK06144',
##            '/data/rob1g10/SPH/DATA/DISKS/M10_34.1d_0.4e/processedDISKs/QGDISK14134'],
#            [1, 25], 115)
#plotDensityProfileFromPartNumbersWithFit(partcoords, rho, udist, umass, [1, 10], 200)
#print max(R)
#plotDensityProfileFromPartNumbers(['/data/rob1g10/SPH/DATA/DISKS/M14_20d_0.0e/processedDISKs/QGDISK35901',
#                                   '/data/rob1g10/SPH/DATA/DISKS/M10_20d_0.0e/processedDISKs/QGDISK30001',
#                                   '/data/rob1g10/SPH/DATA/DISKS/M8_20d_0.0e/processedDISKs/QGDISK36001',
#                                   '/data/rob1g10/SPH/DATA/DISKS/M6_20d_0.0e/processedDISKs/QGDISK36001',
#                                   '/data/rob1g10/SPH/DATA/DISKS/M1.4_20d_0.0e/processedDISKs/QGDISK36001'],
#                                   [1, 25], 200)

#plotDensityProfileFromPartNumbers(['/data/rob1g10/SPH/DATA/DISKS/isoBe_M15_Mdot9/processedDISKs/QGDISK08701',
#                                   '/data/rob1g10/SPH/DATA/DISKS/isoBe_M15_Mdot10/processedDISKs/QGDISK07901',
#                                   '/data/rob1g10/SPH/DATA/DISKS/isoBe_M15_Mdot11/processedDISKs/QGDISK08301'],
#                                   [1, 25], 200)

"""Z dependence"""

#plotZdependence(partcoords, rho, udist, umass, [1, 5], ylims=['b', 'isolated Be'])
#plotZdependence(partcoords, rho, udist, umass, [5, 5.1], ylims=['b', 'isolated Be'])
#plotZdependence(partcoords, rho, udist, umass, [10, 10.1], ylims=['b', 'isolated Be'])
#plotZdependence(partcoords, rho, udist, umass, [15, 15.1], ylims=['b', 'isolated Be'])
#plotZdependence(partcoords, rho, udist, umass, [20, 20.5], ylims=['b', 'isolated Be'])

""" Power law index dependence on R """

#plotPowerLawVariationWithR(partcoords, rho, u, [2, 25], 50, udist, utime, fitType="polyfit")
#plotPowerLawWindowWithR('/data/rob1g10/SPH/DATA/DISKS/20d_0.0e_eqTest/processedDISKs/QGDISK39024', [2, 25], 115)
#plotPowerLawWindowWithR('/data/rob1g10/SPH/DATA/DISKS/20d_0.2e_eqTest/processedDISKs/QGDISK39034', [2, 25], 115)
#plotPowerLawWindowWithR('/data/rob1g10/SPH/DATA/DISKS/20d_0.4e_eqTest/processedDISKs/QGDISK39030', [2, 25], 115)
#plotPowerLawWindowWithR('/data/rob1g10/SPH/DATA/DISKS/40d_0.0e_eqTest/processedDISKs/QGDISK38005', [2, 25], 115)
#plotPowerLawWindowWithR('/data/rob1g10/SPH/DATA/DISKS/40d_0.2e_eqTest/processedDISKs/QGDISK38004', [2, 25], 115)
#plotPowerLawWindowWithR('/data/rob1g10/SPH/DATA/DISKS/40d_0.4e_eqTest/processedDISKs/QGDISK38907', [2, 25], 115)
#plotPowerLawWindowWithR('/data/rob1g10/SPH/DATA/DISKS/60d_0.0e_eqTest/processedDISKs/QGDISK58302', [2, 25], 115)
#plotPowerLawWindowWithR('/data/rob1g10/SPH/DATA/DISKS/60d_0.2e_eqTest/processedDISKs/QGDISK58004', [2, 25], 115)
#plotPowerLawWindowWithR('/data/rob1g10/SPH/DATA/DISKS/60d_0.4e_eqTest/processedDISKs/QGDISK58201', [2, 25], 115)
#plotPowerLawWindowWithR('/data/rob1g10/SPH/DATA/DISKS/80d_0.0e_eqTest/processedDISKs/QGDISK26502', [2, 25], 115)
#plotPowerLawWindowWithR('/data/rob1g10/SPH/DATA/DISKS/80d_0.2e_eqTest/processedDISKs/QGDISK26704', [2, 25], 115)
#plotPowerLawWindowWithR('/data/rob1g10/SPH/DATA/DISKS/80d_0.4e_eqTest/processedDISKs/QGDISK36001', [2, 25], 115)

#findDiskSize(['/data/rob1g10/SPH/DATA/DISKS/20d_0.0e_eqTest/processedDISKs/QGDISK39024',
#              '/data/rob1g10/SPH/DATA/DISKS/20d_0.2e_eqTest/processedDISKs/QGDISK39034',
#              '/data/rob1g10/SPH/DATA/DISKS/20d_0.4e_eqTest/processedDISKs/QGDISK39030',
#              '/data/rob1g10/SPH/DATA/DISKS/40d_0.0e_eqTest/processedDISKs/QGDISK38005',
#              '/data/rob1g10/SPH/DATA/DISKS/40d_0.2e_eqTest/processedDISKs/QGDISK38004',
#              '/data/rob1g10/SPH/DATA/DISKS/40d_0.4e_eqTest/processedDISKs/QGDISK38907',
#              '/data/rob1g10/SPH/DATA/DISKS/60d_0.0e_eqTest/processedDISKs/QGDISK58302',
#              '/data/rob1g10/SPH/DATA/DISKS/60d_0.2e_eqTest/processedDISKs/QGDISK58004',
#              '/data/rob1g10/SPH/DATA/DISKS/60d_0.4e_eqTest/processedDISKs/QGDISK58201',
#              '/data/rob1g10/SPH/DATA/DISKS/80d_0.0e_eqTest/processedDISKs/QGDISK26502',
#              '/data/rob1g10/SPH/DATA/DISKS/80d_0.2e_eqTest/processedDISKs/QGDISK26704',
#              '/data/rob1g10/SPH/DATA/DISKS/80d_0.4e_eqTest/processedDISKs/QGDISK36001',
#              '/data/rob1g10/SPH/DATA/DISKS/100d_0.0e_eqTest/processedDISKs/QGDISK38102',
#              '/data/rob1g10/SPH/DATA/DISKS/100d_0.2e_eqTest/processedDISKs/QGDISK38702',
#              '/data/rob1g10/SPH/DATA/DISKS/100d_0.4e_eqTest/processedDISKs/QGDISK38707'])
#plotPowerLawWindowWithR('/data/rob1g10/SPH/DATA/DISKS/40d_0.2e_eqTest/processedDISKs/QGDISK38701', [2, 15], 50)
#plotPowerLawWindowWithR('/data/rob1g10/SPH/DATA/DISKS/40d_0.4e_eqTest/processedDISKs/QGDISK38701', [2, 15], 50)

""" Disk Flare """

#plotDiskFlare(partcoords, udist, 30)
#Zmaxs = plotFlareLineWithRadius(partcoords, udist, 50)
#x = np.linspace(1, 25, len(Zmaxs))
#p1 = np.polyfit(x, Zmaxs, 1)
#p2 = np.polyfit(x, Zmaxs, 2)
#p3 = np.polyfit(x, Zmaxs, 3)
#p4 = np.polyfit(x, Zmaxs, 4)
#plt.plot(x, p1[0]*x + p1[1])
#plt.plot(x, p2[0]*x**2 + p2[1]*x + p2[2])
#plt.plot(x, p3[0]*x**3 + p3[1]*x**2 + p3[2]*x + p3[3])
#plt.plot(x, p4[0]*x**4 + p4[1]*x**3 + p4[2]*x**2 + p4[3]*x + p4[4], color="red")
#print p4

""" 2D Hist """

#plotDensity2dHist('/data/rob1g10/SPH/DATA/DISKS/isoBe_injFixed/processedDISKs/QGDISK04901',
#             0, 0, show='y', save='n', outFilename='bleh', cbar=[0,8],
#             xylims=[-25,25,-25,25])

""" Ring Density """

#plotRingDensity('/data/rob1g10/SPH/DATA/DISKS/sph3D_f90_radialTemp_q0.5/processedDISKs/QGDISK06001')

""" OTHER """
#plt.scatter(R, 1e-5 * np.sqrt(2*u / 3) * (0.01*udist / utime), s=0.5)
#plt.scatter(R, partve, s=0.5)
#escapeVel = 1e-5 * np.sqrt(2 * 6.67e-11 * 18 * 2e30 / (np.linspace(1,1.2,20) * 7 * 695700e3))
#kepVel = escapeVel / np.sqrt(2)
#plt.plot(np.linspace(1,7,len(kepVel)), escapeVel)
#plt.plot(np.linspace(1,7,len(kepVel)), kepVel)
#plt.xlim(1,10)
plt.show()
