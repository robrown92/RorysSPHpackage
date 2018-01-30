import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import pylab as plt
import os
from scipy.optimize import leastsq, curve_fit
from RorysSPHpackageStripdisks import createAccretionArray
from RorysSPHpackage import readQDISK, readUnits, readFILES, readQADISK, findPeriastron, calcQDISKsForOrbit, findBinaryPhase, findDiskSize

def rotation_matrix(axis, theta):
    """
    Return the rotation matrix associated with counterclockwise rotation about
    the given axis by theta radians.
    """
    axis = np.asarray(axis)
    axis = axis/np.sqrt(np.dot(axis, axis))
    a = np.cos(theta/2.0)
    b, c, d = -axis*np.sin(theta/2.0)
    aa, bb, cc, dd = a*a, b*b, c*c, d*d
    bc, ad, ac, ab, bd, cd = b*c, a*d, a*c, a*b, b*d, c*d
    return np.array([[aa+bb-cc-dd, 2*(bc+ad), 2*(bd-ac)],
                     [2*(bc-ad), aa+cc-bb-dd, 2*(cd+ab)],
                     [2*(bd+ac), 2*(cd-ab), aa+dd-bb-cc]])


def plot_System_NSacc_Obscuration_DiscArea(baseName, pngdirectory, DISKstart, DISKend, DISKsteps,
                                           orientation=[0,0], obscTol=1.0, areaSteps=10):
                                               
    maxRhos = np.array([])
    minRhos = np.array([])
    times = np.array([])
    
    NSacc = np.array([])
    rhoObsc = np.array([])
    discAreaArray = np.array([])
    
    theta = np.deg2rad(-orientation[0])
    phi = np.deg2rad(-orientation[1])
    
    "Assemble full arrays of quantities to plot of each time step"
    for DISKnum in range(DISKstart, DISKend):
        for DISKstep in range(1, DISKsteps):
            filename = '/data/rob1g10/SPH/DATA/DISKS/{}/processedDISKs/QGDISK{:03}{:02}'.format(baseName, DISKnum, DISKstep)
            diskname = filename[:-26] + "disk001"
            IPname = filename[:-11] + "IPDISK" + filename[-5:-2]
            Becoords, NScoords, partcoords, partvels, rho, u, num_parts, time, umass, udist, utime = readFILES(filename, diskname, IPname)
            Becoords *= udist / (7 * 6.96e10)
            NScoords *= udist / (7 * 6.96e10)
            partcoords *= udist / (7 * 6.96e10)
            rho *= umass / (udist ** 3)
        
            xcoords = np.array(partcoords[:,0])
            ycoords = np.array(partcoords[:,1])
            zcoords = np.array(partcoords[:,2])
            NScoords = np.array(NScoords)
            
            times = np.append(times, time * utime / (24 * 3600))
            maxRhos = np.append(maxRhos, np.max(rho))
            minRhos = np.append(minRhos, np.min(rho))
            
            "Particle Rotation------------------------------------------------"    
            
            rotSolution = np.dot(rotation_matrix([0, 0, 1], phi), [xcoords, ycoords, zcoords])
            
            xrot1, yrot1, zrot1 = rotSolution[0], rotSolution[1], rotSolution[2]
            
            rotSolution = np.dot(rotation_matrix([1, 0, 0], theta), [xrot1, yrot1, zrot1])
            
            xSol, ySol, zSol = rotSolution[0], rotSolution[1], rotSolution[2]
            
            "Now NS-----------------------------------------------------------"
            
            rotSolution = np.dot(rotation_matrix([0, 0, 1], phi), [NScoords[0], NScoords[1], NScoords[2]])
            
            NSxrot1, NSyrot1, NSzrot1 = rotSolution[0], rotSolution[1], rotSolution[2]
            
            rotSolution = np.dot(rotation_matrix([1, 0, 0], theta), [NSxrot1, NSyrot1, NSzrot1])
            
            NSxSol, NSySol, NSzSol = rotSolution[0], rotSolution[1], rotSolution[2]
            
            "Obscuration------------------------------------------------------"
            
            rhoSum = np.sum( rho[ (xSol <= NSxSol + obscTol)
                              * (xSol >= NSxSol - obscTol) 
                              * (ySol <= NSySol + obscTol) 
                              * (ySol >= NSySol - obscTol) 
                              * (zSol >= NSzSol)
                            ])            
            
            if(rhoSum < 0):
                rhoSum = 0
            
            rhoObsc = np.append(rhoObsc, rhoSum) 
            
            "Disc Area--------------------------------------------------------"
            "Step through bands of y coordinates and find max and min x coords to find area"
            
            discArea = 0
            stepSize = (np.max(ySol) - np.min(ySol)) / areaSteps
            for ybands in np.linspace(np.min(ySol), np.max(ySol), areaSteps):
                xband = xSol[ (ySol >= ybands) * (ySol < ybands + stepSize) ]
                if(len(xband) > 0):
                    discArea += stepSize * (np.max(xband) - np.min(xband))
            
            discAreaArray = np.append(discAreaArray, discArea)
            
            
    times -= times[0]
    
    "Alternate colors"
    rhomax = np.max(maxRhos)
    rhomin = np.min(minRhos)
    cm = plt.cm.get_cmap('rainbow')
#    cm = plt.cm.get_cmap('hot')
#    colours = np.zeros(len(rho))
#    for i in range(len(rho)):
#        rgrad = round(1 - rho[i]/rhomax, 1)
#        if(rgrad > 1):
#            rgrad = 1
#        colours[i] = rgrad
        
    count = -1
    for DISKnum in range(DISKstart+1, DISKend-1):
        for DISKstep in range(1,DISKsteps):
            count += 1
            if(count+10 >= len(times)):
                break
            print "Making image", count+1, "of", len(times)
            filename = '/data/rob1g10/SPH/DATA/DISKS/{}/processedDISKs/QGDISK{:03}{:02}'.format(baseName, DISKnum, DISKstep)
            diskname = filename[:-26] + "disk001"
            IPname = filename[:-11] + "IPDISK" + filename[-5:-2]
            Becoords, NScoords, partcoords, partvels, rho, u, num_parts, time, umass, udist, utime = readFILES(filename, diskname, IPname)
            Becoords *= udist / (7 * 6.96e10)
            NScoords *= udist / (7 * 6.96e10)
            partcoords *= udist / (7 * 6.96e10)
            rho *= umass / (udist ** 3)
            

            """3D Plot of System"""
        
            fig = plt.figure(figsize=(15, 15))
            ax3D = fig.add_subplot(2, 2, 1, projection='3d')
#            print len(rho), len(xcoords)
            sc = ax3D.scatter(partcoords[:,0], partcoords[:,1], partcoords[:,2], s=0.5, 
                              c=rho, vmin=rhomin, vmax=rhomax, cmap=cm, edgecolor='')
        
            ax3D.scatter(NScoords[0], NScoords[1], NScoords[2], color='b')
        
            ax3D.set_xlim(-11, 11)
            ax3D.set_ylim(-11, 11)
            ax3D.set_zlim(-11, 11)
        
            ax3D.set_xlabel('x')
            ax3D.set_ylabel('y')
            ax3D.set_zlabel('z')
        
            ax3D.view_init(elev=orientation[0]+90, azim=orientation[1]+90)
            
            " NS Accretion "  
            
        #    os.system("ls /data/rob1g10/SPH/DATA/DISKS/{}/disk* > lsdisks.txt".format(filename))
        #    
        #    f = open("lsdisks.txt", 'r')
        #    lines = f.readlines()
        #	
        #    numdisks = int(lines[-1][-3:])
        #    
        #    f.close()
        #
        #    disks = np.array([])        
        #    
        #    for ND in range(1, numdisks+1):
        #        disks = np.append(disks, '/data/rob1g10/SPH/DATA/DISKS/{}/disk{:03}'.format(baseNames[BN], ND))
            
        
        #    accTimes, accArray = createAccretionArray(disks,"NScapture", None, 'day', tol=0.01)
        #    accArray *= 1e-14 * umass / (24*3600)
        #    accTimes = accTimes.astype(float)
            
            axAccr = fig.add_subplot(2, 2, 2)
        
        #    axAccr.plot
        
            axDisc.set_xlabel('Time (days)')
            axDisc.set_ylabel('NS accretion')
        
            " Obscuration "
            
            axObsc = fig.add_subplot(2, 2, 3)
            
            axObsc.scatter(times[count+5], rhoObsc[count+5], color='r')
            axObsc.plot(times[count:count+10], rhoObsc[count:count+10])
            axObsc.set_ylim(np.min(rhoObsc), np.max(rhoObsc))
            
            axDisc.set_xlabel('Time (days)')
            axDisc.set_ylabel('Obscuration')
        
            " Disc Size "
            
            axDisc = fig.add_subplot(2, 2, 4)
                
            axDisc.scatter(times[count+5], discAreaArray[count+5], color='r')
            axDisc.plot(times[count:count+10], discAreaArray[count:count+10])
            axDisc.set_ylim(np.min(discAreaArray), np.max(discAreaArray))
            
            axDisc.set_xlabel('Time (days)')
            axDisc.set_ylabel('Disc Area')
            
            " Save file "
            
            outFilename = "{}AniPart{:04}".format(pngdirectory, count+1)
            plt.savefig(outFilename, figsize=(18, 16), dpi=100)
            
            os.system("convert -delay 20 -loop 0 {}*.png {}.gif".format(pngdirectory, baseName))
    
    return


def plotVisibleDiscArea(baseName, DISKstart, DISKend, DISKsteps,
                        orientation=[0,0], areaSteps=10):

    times = np.array([])
    
    discAreaArray = np.array([])
    
    theta = np.deg2rad(-orientation[0])
    phi = np.deg2rad(-orientation[1])
    
    "Assemble full arrays of quantities to plot of each time step"
    for DISKnum in range(DISKstart, DISKend):
        for DISKstep in range(1, DISKsteps):
            filename = '/data/rob1g10/SPH/DATA/DISKS/{}/processedDISKs/QGDISK{:03}{:02}'.format(baseName, DISKnum, DISKstep)
            diskname = filename[:-26] + "disk001"
            IPname = filename[:-11] + "IPDISK" + filename[-5:-2]
            Becoords, NScoords, partcoords, partvels, rho, u, num_parts, time, umass, udist, utime = readFILES(filename, diskname, IPname)
            Becoords *= udist / (7 * 6.96e10)
            NScoords *= udist / (7 * 6.96e10)
            partcoords *= udist / (7 * 6.96e10)
            rho *= umass / (udist ** 3)
        
            xcoords = np.array(partcoords[:,0])
            ycoords = np.array(partcoords[:,1])
            zcoords = np.array(partcoords[:,2])
            NScoords = np.array(NScoords)
            
            times = np.append(times, time * utime / (24 * 3600))
            
            "Particle Rotation------------------------------------------------"    
            
            rotSolution = np.dot(rotation_matrix([0, 0, 1], phi), [xcoords, ycoords, zcoords])
            
            xrot1, yrot1, zrot1 = rotSolution[0], rotSolution[1], rotSolution[2]
            
            rotSolution = np.dot(rotation_matrix([1, 0, 0], theta), [xrot1, yrot1, zrot1])
            
            xSol, ySol, zSol = rotSolution[0], rotSolution[1], rotSolution[2]   
            
            "Disc Area--------------------------------------------------------"
            "Step through bands of y coordinates and find max and min x coords to find area"
            
            discArea = 0
            stepSize = (np.max(ySol) - np.min(ySol)) / areaSteps
            for ybands in np.linspace(np.min(ySol), np.max(ySol), areaSteps):
                xband = xSol[ (ySol >= ybands) * (ySol < ybands + stepSize) ]
                if(len(xband) > 0):
                    discArea += stepSize * (np.max(xband) - np.min(xband))
            
            discAreaArray = np.append(discAreaArray, discArea)
            
    plt.plot(times, discAreaArray)
    
    return


def plotNSObscuration(baseName, DISKstart, DISKend, DISKsteps,
                        orientation=[0,0], obscTol=1.0):

    times = np.array([])
    
    rhoObsc = np.array([])
    
    theta = np.deg2rad(-orientation[0])
    phi = np.deg2rad(-orientation[1])
    
    "Assemble full arrays of quantities to plot of each time step"
    for DISKnum in range(DISKstart, DISKend):
        for DISKstep in range(1, DISKsteps):
            filename = '/data/rob1g10/SPH/DATA/DISKS/{}/processedDISKs/QGDISK{:03}{:02}'.format(baseName, DISKnum, DISKstep)
            diskname = filename[:-26] + "disk001"
            IPname = filename[:-11] + "IPDISK" + filename[-5:-2]
            Becoords, NScoords, partcoords, partvels, rho, u, num_parts, time, umass, udist, utime = readFILES(filename, diskname, IPname)
            Becoords *= udist / (7 * 6.96e10)
            NScoords *= udist / (7 * 6.96e10)
            partcoords *= udist / (7 * 6.96e10)
            rho *= umass / (udist ** 3)
        
            xcoords = np.array(partcoords[:,0])
            ycoords = np.array(partcoords[:,1])
            zcoords = np.array(partcoords[:,2])
            NScoords = np.array(NScoords)
            
            times = np.append(times, time * utime / (24 * 3600))
            
            "Particle Rotation------------------------------------------------"    
            
            rotSolution = np.dot(rotation_matrix([0, 0, 1], phi), [xcoords, ycoords, zcoords])
            
            xrot1, yrot1, zrot1 = rotSolution[0], rotSolution[1], rotSolution[2]
            
            rotSolution = np.dot(rotation_matrix([1, 0, 0], theta), [xrot1, yrot1, zrot1])
            
            xSol, ySol, zSol = rotSolution[0], rotSolution[1], rotSolution[2]   
            
            "Now NS-----------------------------------------------------------"
            
            rotSolution = np.dot(rotation_matrix([0, 0, 1], phi), [NScoords[0], NScoords[1], NScoords[2]])
            
            NSxrot1, NSyrot1, NSzrot1 = rotSolution[0], rotSolution[1], rotSolution[2]
            
            rotSolution = np.dot(rotation_matrix([1, 0, 0], theta), [NSxrot1, NSyrot1, NSzrot1])
            
            NSxSol, NSySol, NSzSol = rotSolution[0], rotSolution[1], rotSolution[2]
            
            "Obscuration------------------------------------------------------"
            
            rhoSum = np.sum( rho[ (xSol <= NSxSol + obscTol)
                              * (xSol >= NSxSol - obscTol) 
                              * (ySol <= NSySol + obscTol) 
                              * (ySol >= NSySol - obscTol) 
                              * (zSol >= NSzSol)
                            ])            
            
            if(rhoSum < 0):
                rhoSum = 0
            
            rhoObsc = np.append(rhoObsc, rhoSum) 
            
    plt.plot(times, rhoObsc)
    
    return


plot_System_NSacc_Obscuration_DiscArea('SXP5.05_M-11_alpha0.3', 
                                       '/data/rob1g10/SPH/DATA/Animations/SXP5.05_allPlot/',
                                       555, 592, 13, orientation=[-30,110])
                                       
#plotVisibleDiscArea('SXP5.05_M-11_alpha0.3', 555, 592, 13, orientation=[-60,110])

#plotNSObscuration('SXP5.05_M-11_alpha0.3', 555, 592, 13, orientation=[-30,110], obscTol=1.0)
