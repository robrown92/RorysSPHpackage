import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import pylab as plt
import os
from scipy.optimize import leastsq
from RorysSPHpackageStripdisks import stripdiskForDISKs, createAccretionArrayForDISKs
#from RorysSPHpackageOpticalDepth import readOptDepthData4py 

"""
 Columns in ASCII file produced in imagep3DSP that have been discovered,
 there are 15 in total with a number of rows equal to the number of active
 particles.

     columns 1,2 & 3 - x, y & z coords or vels?
     column 10 - h (as described by the code creating it)

"""

"""----------------------------------------------------------------------------

                imagep3DSP stuff

----------------------------------------------------------------------------"""


def create3DSPinput(location, start, end):
    """
    Produces an input file for imagep3DSP.

    INPUT VARIABLES
    ---------------

    start (int) - DISK file to start from (1 - 999)
    
    end (int) - DISK file to end at (1 - 999)

    OUTPUT VARIABLES
    ----------------

    An input file for imagep3DSP.

    """
    
    numOfDISKs = abs(end - start) + 1
    
    f = open(location, 'w')
    
    f.write('{}\n'.format(numOfDISKs))
    
    for i in range(start, end+1):
        f.write('DISK{0:03}\n'.format(i))
        
    f.write('y\n')
    f.write('y\n')
    f.write('n\n')
    f.write('1\n')
    f.write('n\n')
    f.write('0\n')

    f.close()
    
    return

"""----------------------------------------------------------------------------

                PLOTTING FUNCTIONS

----------------------------------------------------------------------------"""

def plotAndSaveFig(filename, time, outFilename, plotType, minRho='n', maxRho='n',
                   orientation=[0,0], show='n', save='n', NS='y', xylims='None', T="None"):
    """
    Produces a scatter plot in 3D or 2D of the particles and stars.

    GEOMETRY
    --------

    z axis is height of disk, y axis is width and x axis is depth.

    INPUT VARIABLES
    ---------------
    
    filename (list of 3 strings) - name of the filename with the particle 
                                   data on it and name of disk file for 
                                   units and then the IPdisk name.
    
    time (float) - time for the disk dump being plotted.

    outFilename (string) - name of the file that the plot will be saved as.

    plotType (string 3d/2d) - chooses between 3D or 2D plot.

    orientation (2 element float array) - gives the polar and azimuthal angles
                                          that the plot will adhere to.
                                          (only in 3d currently)

    show (y / n) - chooses whether to display the plot produced or not.

    OUTPUT VARIABLES
    ----------------

    A plot of the particles.

    """
    
    Becoords, NScoords, partcoords, partvels, rho, u, num_parts = readQDISK(filename[0])
    umass, udist, utime = readUnits(filename[1])
    Becoords *= udist / (7 * 6.96e10)
    NScoords *= udist / (7 * 6.96e10)
    partcoords *= udist / (7 * 6.96e10)
    rho *= umass / (udist ** 3)
    if(T == "None"):
        T = np.zeros_like(rho)

    xcoords = partcoords[:,0]
    ycoords = partcoords[:,1]
    zcoords = partcoords[:,2]
    R = np.sqrt(xcoords**2 + ycoords**2)
    
#    rhomax = max(rho)
#    colours = np.zeros((len(rho), 3))
#    for i in range(len(rho)):
#        rgrad = round(1 - rho[i]/rhomax, 1)
#        if(rgrad > 1):
#            rgrad = 1
#        colours[i] = [rgrad, 0, 0]

    if(maxRho !='n'):
        rhomax = maxRho
    else:
        rhomax = max(rho)
    if(maxRho !='n'):
        rhomin = minRho
    elif(min(rho) < 0):
        rhomin = 0
    else:
        rhomin = min(rho)
        
    "Alternate colors"
    
    cm = plt.cm.get_cmap('rainbow')
    colours = np.zeros(len(rho))
    for i in range(len(rho)):
        rgrad = round(1 - rho[i]/rhomax, 1)
        if(rgrad > 1):
            rgrad = 1
        colours[i] = rgrad

    if(plotType == '2d'):

        """ Density Scatter """
        partx = np.zeros(len(partcoords[:,0]))
        party = np.zeros(len(partcoords[:,0]))
        partz = np.zeros(len(partcoords[:,0]))
        NSx = NScoords[0]
        NSy = NScoords[1]
        NSz = NScoords[2]
        
        phi = np.deg2rad(orientation[0])
        theta = np.deg2rad(orientation[1])
        
#        " Rotate in z axis "
#        for i in range(len(partcoords)):
#            partx[i] = xcoords[i] * np.cos(phi) - ycoords[i] * np.sin(phi)
#            party[i] = xcoords[i] * np.sin(phi) + ycoords[i] * np.cos(phi)
#        
#        " Rotate in y axis "
#        for i in range(len(partcoords)):
#            party[i] = partcoords[i,1] * np.cos(theta) - partcoords[i,2] * np.sin(theta)
#            partz[i] = partcoords[i,1] * np.sin(theta) + partcoords[i,2] * np.cos(theta)

        
        NSx = NScoords[0] * np.cos(phi) - NScoords[1] * np.sin(phi)
        NSy = NScoords[0] * np.sin(phi) + NScoords[1] * np.cos(phi) 
        NSy = NScoords[1] * np.cos(theta) - NScoords[2] * np.sin(theta)
        NSz = NScoords[1] * np.sin(theta) + NScoords[2] * np.cos(theta)
        plt.figure(figsize=(18, 16))

#        plt.scatter(xcoords * np.cos(phi) - ycoords * np.sin(phi), 
#                    (xcoords * np.sin(phi) + ycoords * np.cos(phi))*np.cos(theta) - zcoords * np.sin(theta), 
#                    color=colours, s = 0.25)
        sc = plt.scatter(xcoords * np.cos(phi) - ycoords * np.sin(phi), 
#                         partcoords[:,1] * np.sin(theta) + partcoords[:,2] * np.cos(theta) ,
                    (xcoords * np.sin(phi) + ycoords * np.cos(phi))*np.cos(theta) - zcoords * np.sin(theta), 
                    c=rho, s = 1, vmin=rhomin, vmax=rhomax, cmap=cm, edgecolor='')
        cb = plt.colorbar(sc)
        cb.ax.set_ylabel("Density in g per cm cubed")
#        plt.scatter(Becoords[0], Becoords[1], color='b')
        if(NS == 'y'):
            plt.scatter(NScoords[0] * np.cos(phi) - NScoords[1] * np.sin(phi), 
                        (NScoords[0] * np.sin(phi) + NScoords[1] * np.cos(phi))*np.cos(theta) - NScoords[2] * np.sin(theta),
                        color='b')
        else:
            plt.scatter(NScoords[0] * np.cos(phi) - NScoords[1] * np.sin(phi), 
                        (NScoords[0] * np.sin(phi) + NScoords[1] * np.cos(phi))*np.cos(theta) - NScoords[2] * np.sin(theta),
                        color='r', s=0.25)
        #plt.text(-14.0, 14.0, 't = {} days'.format(int(time) - 1856), size=20)
        if(xylims != 'None'):
            plt.xlim(xylims[0], xylims[1])
            plt.ylim(xylims[2], xylims[3])
        plt.xlabel(r"R / Be star radius")
        plt.ylabel(r"R / Be star radius")
#        xcirc2 = [16.84 * np.sin(np.radians(i)) for i in range(0, 360)]
#        ycirc2 = [16.84 * np.cos(np.radians(i)) for i in range(0, 360)]
#        plt.plot(xcirc2, ycirc2, label='n=2')
#        xcirc3 = [14.75 * np.sin(np.radians(i)) for i in range(0, 360)]
#        ycirc3 = [14.75 * np.cos(np.radians(i)) for i in range(0, 360)]
#        plt.plot(xcirc3, ycirc3, label='n=3.5')
#        xcirc4 = [14.32 * np.sin(np.radians(i)) for i in range(0, 360)]
#        ycirc4 = [14.32 * np.cos(np.radians(i)) for i in range(0, 360)]
#        plt.plot(xcirc4, ycirc4, label='n=4')
#        plt.legend()
        
        
        " Just Be "
        
#        plt.scatter(NScoords[0], NScoords[1], color='r', s=0.25)

        " For looking at injected particles "

        #plt.xlim(-50, 100)
        #plt.ylim(-100, 100)
#
#        plt.scatter(xcoords[2:98652], ycoords[2:98652], color=colours, s = 0.5)
#        plt.scatter(xcoords[98652:], ycoords[98652:], color='g', s = 2)

        plt.gca().set_aspect('equal', adjustable='box')
    
        if(show == 'y'):
            plt.show()

    elif(plotType == '3d'):
        
        R = np.sqrt(partcoords[:,0]**2 + partcoords[:,1]**2)

        """3D Plot of System"""

        fig = plt.figure(figsize=(15, 15))
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter(xcoords[2:], ycoords[2:], zcoords[2:], color=colours, s=0.5)
        """ Use to see injected particles """
        #ax.scatter(xcoords[2:20536], ycoords[2:20536], zcoords[2:20536], color='r', s=0.5)
        #ax.scatter(xcoords[20536:], ycoords[20536:], zcoords[20536:], color='g', s=0.5)

#        ax.scatter(Becoords[0], Becoords[1], Becoords[2], color='b')
        if(NS == 'y'):
            ax.scatter(NScoords[0], NScoords[1], NScoords[2], color='b')
        else:
            ax.scatter(NScoords[0], NScoords[1], NScoords[2], color='r', s=0.5)

        ax.set_xlim(-int(max(R)), int(max(R)))
        ax.set_ylim(-int(max(R)), int(max(R)))
        ax.set_zlim(-int(max(R)), int(max(R)))

        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')

        ax.view_init(elev=orientation[0], azim=orientation[1])

        if(show == 'y'):
            plt.show()
    
    elif(plotType == 'plusHistDepths'):
        """ Density Scatter """
        
        f, axarr = plt.subplots(2, figsize=(16, 16))

        axarr[0].scatter(xcoords[2:], ycoords[2:], color=colours, s = 0.25)
#        axarr[0].scatter(Becoords[0], Becoords[1], color='b', s=80)
        axarr[0].scatter(NScoords[0], NScoords[1], color='b')
#        axarr[0].text(-25.0, 13.0, 't = {} days'.format(int(time) - 1965), size=20)
        axarr[0].set_xlim(-25, 25)
        axarr[0].set_ylim(-11, 11)
        axarr[0].set_xlabel(r"R / Be star radius")
        axarr[0].set_ylabel(r"R / Be star radius")

        """ Histogram of vels """
        
        depth = [3e4, 1.7e4, 1.14e4, 5.8e3] #20d 0.4e
        
        for i in range(0,4):
            
            Becoords, NScoords, partcoords, partvels, rho, u, num_parts = readQDISK(filename)
        
            Becoords /= udist / (7 * 6.96e10)
            NScoords /= udist / (7 * 6.96e10)
            partcoords /= udist / (7 * 6.96e10)
            
            unitvel = 0.01 * udist / utime
            
            diskMaxZ = 0.5
            diskMaxY = 2.

            partvelsbin = binVels(partcoords, Becoords, partvels,
                                  depth[i], 45, [10, 10], diskMaxY, diskMaxZ)
    #        vels = partvelsbin * unitvel
        
            angstroms = partvelsbin * unitvel / 3e8   
            angstroms += 1
            angstroms *= 6562.8
            counts, bins = np.histogram(angstroms, bins=int(300 * float(len(partvelsbin)) / len(partvels)))
            print int(100 * float(len(partvelsbin)) / len(partvels))
            
    #        axarr[1].hist(partvelsbin, bins=100, normed='True')
            axarr[1].plot(bins[0:-1], counts, label="{}%".format(int(100 * float(len(partvelsbin)) / len(partvels))))
            axarr[1].set_xlabel('Angstroms')
            axarr[1].set_ylabel('Counts')
#            axarr[1].set_xlim(6550, 6575)
    #        axarr[1].set_ylim(0, 0.5)
        axarr[1].legend()

        if(show == 'y'):
            plt.show()
            
    elif(plotType == 'plusHist'):
        """ Density Scatter """
        
        f, axarr = plt.subplots(1, 2, figsize=(22, 8))

        """ Density Scatter """
        partx = np.zeros(len(partcoords[:,0]))
        party = np.zeros(len(partcoords[:,0]))
        partz = np.zeros(len(partcoords[:,0]))
        NSx = NScoords[0]
        NSy = NScoords[1]
        NSz = NScoords[2]
        
        phi = np.deg2rad(orientation[0])
        theta = np.deg2rad(orientation[1])
        
#        " Rotate in z axis "
#        for i in range(len(partcoords)):
#            partx[i] = xcoords[i] * np.cos(phi) - ycoords[i] * np.sin(phi)
#            party[i] = xcoords[i] * np.sin(phi) + ycoords[i] * np.cos(phi)
#        
#        " Rotate in y axis "
#        for i in range(len(partcoords)):
#            party[i] = partcoords[i,1] * np.cos(theta) - partcoords[i,2] * np.sin(theta)
#            partz[i] = partcoords[i,1] * np.sin(theta) + partcoords[i,2] * np.cos(theta)

        
        NSx = NScoords[0] * np.cos(phi) - NScoords[1] * np.sin(phi)
        NSy = NScoords[0] * np.sin(phi) + NScoords[1] * np.cos(phi) 
        NSy = NScoords[1] * np.cos(theta) - NScoords[2] * np.sin(theta)
        NSz = NScoords[1] * np.sin(theta) + NScoords[2] * np.cos(theta)

#        plt.scatter(xcoords * np.cos(phi) - ycoords * np.sin(phi), 
#                    (xcoords * np.sin(phi) + ycoords * np.cos(phi))*np.cos(theta) - zcoords * np.sin(theta), 
#                    color=colours, s = 0.25)
        sc = axarr[0].scatter(xcoords * np.cos(phi) - ycoords * np.sin(phi), 
                    (xcoords * np.sin(phi) + ycoords * np.cos(phi))*np.cos(theta) - zcoords * np.sin(theta), 
                    c=rho, s = 1, vmin=rhomin, vmax=rhomax, cmap=cm, edgecolor='')
        cb = plt.colorbar(sc, ax=axarr[0])
        cb.ax.set_ylabel("Density in g per cm cubed")
        #axarr[0].scatter(NScoords[0] * np.cos(phi) - NScoords[1] * np.sin(phi), 
        #            (NScoords[0] * np.sin(phi) + NScoords[1] * np.cos(phi))*np.cos(theta) - NScoords[2] * np.sin(theta),
        #            color='b')
                    
        if(xylims != 'None'):
            axarr[0].set_xlim(xylims[0], xylims[1])
            axarr[0].set_ylim(xylims[2], xylims[3])
        axarr[0].set_xlabel(r"R / Be star radius")
        axarr[0].set_ylabel(r"R / Be star radius")

        """ Histogram of vels """
        
        phi = 0
        #theta = 45
        
        BINS = 200
        
        """ Arrange velocities about y-axis so we only require z component"""
    
        partvelx = partvels[:,0] * np.cos(phi) - partvels[:, 1] * np.sin(phi)
        partvely = (partvels[:, 0] * np.sin(phi) + partvels[:, 1] * np.cos(phi)) * np.cos(theta) - partvels[:,2] * np.sin(theta)
        partvelz = partvels[:, 1] * np.sin(theta) + partvels[:, 2] * np.cos(theta)
        
        unitvel = 0.01 * udist / utime
        
        angstroms = partvelz * unitvel / 3e8   
        angstroms += 1
        angstroms *= 6562.8
        
        counts, bins = np.histogram(angstroms, bins=BINS)
        binDiff = bins[1] - bins[0]
        
        " Sort angstroms into the bins manually and scale based on emissivity and transmissivity "
        countsPlot = np.zeros_like(counts)
        for i in range(len(angstroms)):
            binNum = (angstroms[i] - bins[0]) / binDiff
            #print(int(binNum))
            if(int(binNum) < BINS-1):
                countsPlot[int(binNum)] += 1
            else:
                countsPlot[BINS-1] += 1
                
        countsPlot = countsPlot / float(np.sum(countsPlot))
        
        axarr[1].plot(bins[0:-1], countsPlot)#, label=name)
    
        
#        axarr[1].plot(bins[0:-1], counts)
#        axarr[1].hist(angstroms, bins=100)
        axarr[1].set_xlabel('Angstroms')
#        axarr[1].set_ylabel('Counts')
        axarr[1].set_ylim(0, 0.015)
        #axarr[1].set_xlim(6556, 6570)
        axarr[1].set_xlim(6545, 6580)
        axarr[1].set_yticks([])
        axarr[1].legend()

        if(show == 'y'):
            plt.show()
            
    elif(plotType == 'justHist'):
        """ Histogram of vels """
        
        depth = [3e4, 1.7e4, 1.14e4, 5.8e3] #20d 0.4e
        
        for i in range(0,4):
            
            Becoords, NScoords, partcoords, partvels, rho, u, num_parts = readQDISK(filename)
        
            Becoords /= udist / (7 * 6.96e10)
            NScoords /= udist / (7 * 6.96e10)
            partcoords /= udist / (7 * 6.96e10)
            
            unitvel = 0.01 * udist / utime
            
            diskMaxZ = 0.5
            diskMaxY = 2.

            partvelsbin = binVels(partcoords, Becoords, partvels,
                                  depth[i], 45, [10, 10], diskMaxY, diskMaxZ)
    #        vels = partvelsbin * unitvel
        
            angstroms = partvelsbin * unitvel / 3e8   
            angstroms += 1
            angstroms *= 6562.8
            counts, bins = np.histogram(angstroms, bins=int(200 * float(len(partvelsbin)) / len(partvels)))
            print int(100 * float(len(partvelsbin)) / len(partvels))
            
    #        axarr[1].hist(partvelsbin, bins=100, normed='True')
            plt.plot(bins[0:-1], counts, label="{}%".format(int(100 * float(len(partvelsbin)) / len(partvels))))
            plt.xlabel('Angstroms')
            plt.ylabel('Counts')
            plt.xlim(6550, 6575)
        plt.legend()

        if(show == 'y'):
            plt.show()
            
    elif(plotType == 'NSaccretion'):
        
        plt.figure(figsize=(18,8))
    
        """ Density Scatter """
        partx = np.zeros(len(partcoords[:,0]))
        party = np.zeros(len(partcoords[:,0]))
        partz = np.zeros(len(partcoords[:,0]))
        NSx = NScoords[0]
        NSy = NScoords[1]
        NSz = NScoords[2]
        
        phi = np.deg2rad(orientation[0])
        theta = np.deg2rad(orientation[1])
        
#        " Rotate in z axis "
#        for i in range(len(partcoords)):
#            partx[i] = xcoords[i] * np.cos(phi) - ycoords[i] * np.sin(phi)
#            party[i] = xcoords[i] * np.sin(phi) + ycoords[i] * np.cos(phi)
#        
#        " Rotate in y axis "
#        for i in range(len(partcoords)):
#            party[i] = partcoords[i,1] * np.cos(theta) - partcoords[i,2] * np.sin(theta)
#            partz[i] = partcoords[i,1] * np.sin(theta) + partcoords[i,2] * np.cos(theta)

        
        NSx = NScoords[0] * np.cos(phi) - NScoords[1] * np.sin(phi)
        NSy = NScoords[0] * np.sin(phi) + NScoords[1] * np.cos(phi) 
        NSy = NScoords[1] * np.cos(theta) - NScoords[2] * np.sin(theta)
        NSz = NScoords[1] * np.sin(theta) + NScoords[2] * np.cos(theta)

#        plt.scatter(xcoords * np.cos(phi) - ycoords * np.sin(phi), 
#                    (xcoords * np.sin(phi) + ycoords * np.cos(phi))*np.cos(theta) - zcoords * np.sin(theta), 
#                    color=colours, s = 0.25)
        plt.subplot(121)
        sc = plt.scatter(xcoords * np.cos(phi) - ycoords * np.sin(phi), 
                    (xcoords * np.sin(phi) + ycoords * np.cos(phi))*np.cos(theta) - zcoords * np.sin(theta), 
                    c=rho, s = 1, vmin=rhomin, vmax=rhomax, cmap=cm, edgecolor='')
        cb = plt.colorbar(sc)
        cb.ax.set_ylabel("Density in g per cm cubed")
        plt.scatter(NScoords[0] * np.cos(phi) - NScoords[1] * np.sin(phi), 
                    (NScoords[0] * np.sin(phi) + NScoords[1] * np.cos(phi))*np.cos(theta) - NScoords[2] * np.sin(theta),
                    color='b')
                    
        if(xylims != 'None'):
            plt.xlim(xylims[0], xylims[1])
            plt.ylim(xylims[2], xylims[3])
        plt.xlabel(r"R / Be star radius")
        plt.ylabel(r"R / Be star radius")
    
        """ Histogram of vels """
            
        Becoords, NScoords, partcoords, partvels, rho, u, num_parts = readQDISK(filename[0])
    
        Becoords /= udist / (7 * 6.96e10)
        NScoords /= udist / (7 * 6.96e10)
        partcoords /= udist / (7 * 6.96e10)
        
#        phi = np.deg2rad(90)
#        theta = np.deg2rad(90)
#        
#        BINS=100
#        
#        partvelx = partvels[:,0] * np.cos(phi) - partvels[:, 1] * np.sin(phi)
#        partvely = (partvels[:, 0] * np.sin(phi) + partvels[:, 1] * np.cos(phi)) * np.cos(theta) - partvels[:,2] * np.sin(theta)
#        partvelz = partvels[:, 1] * np.sin(theta) + partvels[:, 2] * np.cos(theta)
#        
#        unitvel = 0.01 * udist / utime
#        
#        angstroms = partvelz * unitvel / 3e8   
#        angstroms += 1
#        angstroms *= 6562.8
#        
#        counts, bins = np.histogram(angstroms, bins=BINS)
#        binDiff = bins[1] - bins[0]
#        
#        baseName = "20d_0.4e"
#        FileName = "/data/rob1g10/SPH/DATA/DISKS/{}/optDepthData{}_theta{}_phi{}_time{:0.1f}.txt".format(baseName, baseName, 90, 90, time)# * utime / (24 * 3600))
#        
#        oD = np.array([])
#        " Open file "
#        f = open(FileName, 'r')
#    
#        for line in f:
#            line = line.strip()
#            oD = np.append(oD, float(line))
#
#        f.close()
#        
#        T = np.exp(-oD) / 1e-3
#        
#        " Sort angstroms into the bins manually and scale based on emissivity and transmissivity "
#        countsPlot = np.zeros_like(counts)
#        for i in range(len(angstroms)):
#            binNum = (angstroms[i] - bins[0]) / binDiff
#            if(int(binNum) < BINS-1):
#                try:
#                    countsPlot[int(binNum)] += T[i] # emissScaling[i]
#                except:
#                    pass
#            else:
#                countsPlot[BINS-1] += T[i] #* emissScaling[i]
#                
#        countsPlot = countsPlot / float(np.sum(countsPlot))
#        
#        plt.subplot(224)
#        plt.plot(bins[0:-1], countsPlot)#, label=name)
#    
#        plt.xlabel('Angstroms')
##        plt.ylabel('Counts')
#        plt.yticks([])
#        plt.ylim(0, 0.03)
#        plt.xlim(6540, 6580)
        
        """ NS Accretion """
    
        plt.subplot(122)
        # 8 values per QGDISK
        NSarray = createAccretionArrayForDISKs(['/data/rob1g10/SPH/DATA/DISKS/M10_34.1d_0.4e/disk001', 
                                                '/data/rob1g10/SPH/DATA/DISKS/M10_34.1d_0.4e/disk002'],
                                                'NS', 141, 142)
        NSarray *= 1e-15 * 0.01 * umass / 2e30
        
        QGDISKNum = int(filename[0][-5:-2])
        QGstepNum = int(filename[0][-2:])
        currentElement = (QGDISKNum - 265) * 8 + QGstepNum
        
        IPfilenames = np.array([])
        accreteTimes = np.array([])
        inputdirectory = "/data/rob1g10/SPH/DATA/Animations/M10_34.1d_0.4e/processedDISKs/"
    
        os.system("ls {}IP* > lsIP.txt".format(inputdirectory))
        
        """ IPDISKs """    
        
        f = open("lsIP.txt", "r")
        
        lines = f.readlines()
    
        for line in lines:
            f2 = open("{}".format(line[0:-1]))
            temp = f2.readlines()
            f2.close()
            if(len(temp) != 0):
                IPfilenames = np.append(IPfilenames, line[0:-1])
            
        f.close() 
        
        os.system("rm lsIP.txt")
        
        for IP in IPfilenames:
            accreteTimes = np.append(accreteTimes, readIPDISK(IP))
        accreteTimes *= utime/(24*3600)
        
        plt.plot(accreteTimes[:currentElement] - int(min(accreteTimes)), NSarray[:currentElement])  
        plt.ylim(0, 102 * 1e-15 * 0.01 * umass / 2e30)
        plt.xlim(0, int(max(accreteTimes)) - int(min(accreteTimes)))
        plt.xlabel("Days")
        plt.ylabel("NS Accretion in Solar Masses")
#        print accreteTimes[QGstepNum:QGstepNum+8]
#        print NSarray[currentElement-8:currentElement]                                                 
        
    
        if(show == 'y'):
            plt.show()

    elif(plotType == 'OptDepth'):

        """ Density Scatter """
        partx = np.zeros(len(partcoords[:,0]))
        party = np.zeros(len(partcoords[:,0]))
        partz = np.zeros(len(partcoords[:,0]))
        NSx = NScoords[0]
        NSy = NScoords[1]
        NSz = NScoords[2]
        
        phi = np.deg2rad(orientation[0])
        theta = np.deg2rad(orientation[1])
        
#        " Rotate in z axis "
#        for i in range(len(partcoords)):
#            partx[i] = xcoords[i] * np.cos(phi) - ycoords[i] * np.sin(phi)
#            party[i] = xcoords[i] * np.sin(phi) + ycoords[i] * np.cos(phi)
#        
#        " Rotate in y axis "
#        for i in range(len(partcoords)):
#            party[i] = partcoords[i,1] * np.cos(theta) - partcoords[i,2] * np.sin(theta)
#            partz[i] = partcoords[i,1] * np.sin(theta) + partcoords[i,2] * np.cos(theta)

        
        NSx = NScoords[0] * np.cos(phi) - NScoords[1] * np.sin(phi)
        NSy = NScoords[0] * np.sin(phi) + NScoords[1] * np.cos(phi) 
        NSy = NScoords[1] * np.cos(theta) - NScoords[2] * np.sin(theta)
        NSz = NScoords[1] * np.sin(theta) + NScoords[2] * np.cos(theta)
        plt.figure(figsize=(18, 16))

        "rotate arrays"        
        
        xcoords = partcoords[:,0] * np.cos(phi) - partcoords[:,1] * np.sin(phi)
        ycoords = (partcoords[:,0] * np.sin(phi) + partcoords[:,1] * np.cos(phi))*np.cos(theta) - partcoords[:,2] * np.sin(theta)
        zcoords = partcoords[:,1] * np.sin(theta) + partcoords[:,2] * np.cos(theta)
        
#        """ sort arrays """
#        
#        sortOrder = np.argsort(zcoords)[::-1]
#        
#        xcoords = xcoords[sortOrder]
#        ycoords = ycoords[sortOrder]
#        zcoords = zcoords[sortOrder]
        
        "delete low emission particles"
        
#        delArr = np.array([])
#        for i in range(len(T)):
#            if(T[i] <= 1e-3):
#            if(R[i] > 10):
#                delArr = np.append(delArr, i)
                
#        xcoords = np.delete(xcoords, delArr)
#        ycoords = np.delete(ycoords, delArr)
#        zcoords = np.delete(zcoords, delArr)
#        T = np.delete(T, delArr)
        
#        print len(delArr), len(xcoords)
        
        "reverse order of all arrays"
        #xcoords = xcoords[-1::-1]
        #ycoords = ycoords[-1::-1]
        #zcoords = zcoords[-1::-1]
        #T = T[-1::-1]
#        print len(zcoords)
        
#        plt.scatter(xcoords * np.cos(phi) - ycoords * np.sin(phi), 
#                    (xcoords * np.sin(phi) + ycoords * np.cos(phi))*np.cos(theta) - zcoords * np.sin(theta), 
#                    color=colours, s = 0.25)
        sc = plt.scatter(xcoords, ycoords, 
                    c=T, s = 2, vmin=0, vmax=1, cmap=cm, edgecolor='')
        cb = plt.colorbar(sc)
        cb.ax.set_ylabel("Transmission")
#        plt.scatter(Becoords[0], Becoords[1], color='b')
        if(NS == 'y'):
            plt.scatter(NScoords[0] * np.cos(phi) - NScoords[1] * np.sin(phi), 
                        (NScoords[0] * np.sin(phi) + NScoords[1] * np.cos(phi))*np.cos(theta) - NScoords[2] * np.sin(theta),
                        color='b')
        else:
            plt.scatter(NScoords[0] * np.cos(phi) - NScoords[1] * np.sin(phi), 
                        (NScoords[0] * np.sin(phi) + NScoords[1] * np.cos(phi))*np.cos(theta) - NScoords[2] * np.sin(theta),
                        color='r', s=0.25)
        #plt.text(-14.0, 14.0, 't = {} days'.format(int(time) - 1856), size=20)
        if(xylims != 'None'):
            plt.xlim(xylims[0], xylims[1])
            plt.ylim(xylims[2], xylims[3])
        plt.xlabel(r"R / Be star radius")
        plt.ylabel(r"R / Be star radius")

        plt.gca().set_aspect('equal', adjustable='box')
    
        if(show == 'y'):
            plt.show()

    elif(plotType == 'OptDepthZ'):

        """ Density Scatter """
        partx = np.zeros(len(partcoords[:,0]))
        party = np.zeros(len(partcoords[:,0]))
        partz = np.zeros(len(partcoords[:,0]))
        NSx = NScoords[0]
        NSy = NScoords[1]
        NSz = NScoords[2]
        
        phi = np.deg2rad(orientation[0])
        theta = np.deg2rad(orientation[1])
        
#        " Rotate in z axis "
#        for i in range(len(partcoords)):
#            partx[i] = xcoords[i] * np.cos(phi) - ycoords[i] * np.sin(phi)
#            party[i] = xcoords[i] * np.sin(phi) + ycoords[i] * np.cos(phi)
#        
#        " Rotate in y axis "
#        for i in range(len(partcoords)):
#            party[i] = partcoords[i,1] * np.cos(theta) - partcoords[i,2] * np.sin(theta)
#            partz[i] = partcoords[i,1] * np.sin(theta) + partcoords[i,2] * np.cos(theta)

        
        NSx = NScoords[0] * np.cos(phi) - NScoords[1] * np.sin(phi)
        NSy = NScoords[0] * np.sin(phi) + NScoords[1] * np.cos(phi) 
        NSy = NScoords[1] * np.cos(theta) - NScoords[2] * np.sin(theta)
        NSz = NScoords[1] * np.sin(theta) + NScoords[2] * np.cos(theta)
        plt.figure(figsize=(18, 16))

        "rotate arrays"        
        
        xcoords = partcoords[:,0] * np.cos(phi) - partcoords[:,1] * np.sin(phi)
        ycoords = (partcoords[:,0] * np.sin(phi) + partcoords[:,1] * np.cos(phi))*np.cos(theta) - partcoords[:,2] * np.sin(theta)
        zcoords = partcoords[:,1] * np.sin(theta) + partcoords[:,2] * np.cos(theta)
        
#        """ sort arrays """
#        
#        sortOrder = np.argsort(zcoords)[::-1]
#        
#        xcoords = xcoords[sortOrder]
#        ycoords = ycoords[sortOrder]
#        zcoords = zcoords[sortOrder]
        
        "delete low emission particles"
        
#        delArr = np.array([])
#        for i in range(len(T)):
#            if(T[i] <= 1e-3):
#            if(R[i] > 10):
#                delArr = np.append(delArr, i)
                
#        xcoords = np.delete(xcoords, delArr)
#        ycoords = np.delete(ycoords, delArr)
#        zcoords = np.delete(zcoords, delArr)
#        T = np.delete(T, delArr)
        
#        print len(delArr), len(xcoords)
        
        "reverse order of all arrays"
        xcoords = xcoords[-1::-1]
        ycoords = ycoords[-1::-1]
        zcoords = zcoords[-1::-1]
        T = T[-1::-1]
#        print len(zcoords)

        "sort specific range of particles along x=0"

        delArr = np.array([])
        for i in range(len(xcoords)):
            if(xcoords[i] > 5e-2 or xcoords[i] < -5e-2):# or np.sqrt(xcoords[i]**2 + ycoords[i]**2) < 1):
                delArr = np.append(delArr, i)
                
        xcoords = np.delete(xcoords, delArr)
        ycoords = np.delete(ycoords, delArr)
        zcoords = np.delete(zcoords, delArr)
        T = np.delete(T, delArr)        
        
#        plt.scatter(xcoords * np.cos(phi) - ycoords * np.sin(phi), 
#                    (xcoords * np.sin(phi) + ycoords * np.cos(phi))*np.cos(theta) - zcoords * np.sin(theta), 
#                    color=colours, s = 0.25)
        sc = plt.scatter(np.sqrt(xcoords**2 + ycoords**2), zcoords, 
                    c=T, s = 50, vmin=0, vmax=max(T), cmap=cm, edgecolor='')
        cb = plt.colorbar(sc)
        cb.ax.set_ylabel("Transmission")
#        plt.scatter(Becoords[0], Becoords[1], color='b')
        #plt.text(-14.0, 14.0, 't = {} days'.format(int(time) - 1856), size=20)
        
        plt.xlim(-15, 15)
        plt.ylim(-15, 15)
        plt.xlabel(r"R / Be star radius")
        plt.ylabel(r"R / Be star radius")

        #plt.gca().set_aspect('equal', adjustable='box')
    
        if(show == 'y'):
            plt.show()

    elif(plotType == 'OptDepthPlusHist'):

        f, axarr = plt.subplots(1, 2, figsize=(16, 8))

        """ Density Scatter """
        partx = np.zeros(len(partcoords[:,0]))
        party = np.zeros(len(partcoords[:,0]))
        partz = np.zeros(len(partcoords[:,0]))
        NSx = NScoords[0]
        NSy = NScoords[1]
        NSz = NScoords[2]
        
        phi = np.deg2rad(orientation[0])
        theta = np.deg2rad(orientation[1])
        
#        " Rotate in z axis "
#        for i in range(len(partcoords)):
#            partx[i] = xcoords[i] * np.cos(phi) - ycoords[i] * np.sin(phi)
#            party[i] = xcoords[i] * np.sin(phi) + ycoords[i] * np.cos(phi)
#        
#        " Rotate in y axis "
#        for i in range(len(partcoords)):
#            party[i] = partcoords[i,1] * np.cos(theta) - partcoords[i,2] * np.sin(theta)
#            partz[i] = partcoords[i,1] * np.sin(theta) + partcoords[i,2] * np.cos(theta)

        
        NSx = NScoords[0] * np.cos(phi) - NScoords[1] * np.sin(phi)
        NSy = NScoords[0] * np.sin(phi) + NScoords[1] * np.cos(phi) 
        NSy = NScoords[1] * np.cos(theta) - NScoords[2] * np.sin(theta)
        NSz = NScoords[1] * np.sin(theta) + NScoords[2] * np.cos(theta)

        "rotate arrays"        
        
        xcoords = partcoords[:,0] * np.cos(phi) - partcoords[:,1] * np.sin(phi)
        ycoords = (partcoords[:,0] * np.sin(phi) + partcoords[:,1] * np.cos(phi))*np.cos(theta) - partcoords[:,2] * np.sin(theta)
        zcoords = partcoords[:,1] * np.sin(theta) + partcoords[:,2] * np.cos(theta)
        
#        """ sort arrays """
#        
#        sortOrder = np.argsort(zcoords)[::-1]
#        
#        xcoords = xcoords[sortOrder]
#        ycoords = ycoords[sortOrder]
#        zcoords = zcoords[sortOrder]
        
        "delete low emission particles"
        
#        delArr = np.array([])
#        for i in range(len(T)):
#            if(T[i] <= 1e-3):
#            if(R[i] > 10):
#                delArr = np.append(delArr, i)
                
#        xcoords = np.delete(xcoords, delArr)
#        ycoords = np.delete(ycoords, delArr)
#        zcoords = np.delete(zcoords, delArr)
#        T = np.delete(T, delArr)
        
#        print len(delArr), len(xcoords)
        
        "reverse order of all arrays"
        xcoords = xcoords[-1::-1]
        ycoords = ycoords[-1::-1]
        zcoords = zcoords[-1::-1]
        T = T[-1::-1]
#        print len(zcoords)
        
#        plt.scatter(xcoords * np.cos(phi) - ycoords * np.sin(phi), 
#                    (xcoords * np.sin(phi) + ycoords * np.cos(phi))*np.cos(theta) - zcoords * np.sin(theta), 
#                    color=colours, s = 0.25)
        sc = axarr[0].scatter(xcoords, ycoords, 
                    c=rho, s = 2, vmin=0, vmax=1, cmap=cm, edgecolor='')
        cb = plt.colorbar(sc, ax=axarr[0])
        cb.ax.set_ylabel("Density in g per cm cubed")
#        plt.scatter(Becoords[0], Becoords[1], color='b')
        if(NS == 'y'):
            axarr[0].scatter(NScoords[0] * np.cos(phi) - NScoords[1] * np.sin(phi), 
                        (NScoords[0] * np.sin(phi) + NScoords[1] * np.cos(phi))*np.cos(theta) - NScoords[2] * np.sin(theta),
                        color='b')
        else:
            axarr[0].scatter(NScoords[0] * np.cos(phi) - NScoords[1] * np.sin(phi), 
                        (NScoords[0] * np.sin(phi) + NScoords[1] * np.cos(phi))*np.cos(theta) - NScoords[2] * np.sin(theta),
                        color='r', s=0.25)
        #plt.text(-14.0, 14.0, 't = {} days'.format(int(time) - 1856), size=20)
        if(xylims != 'None'):
            axarr[0].set_xlim(xylims[0], xylims[1])
            axarr[0].set_ylim(xylims[2], xylims[3])
        axarr[0].set_xlabel(r"R / Be star radius")
        axarr[0].set_ylabel(r"R / Be star radius")

        """ Histogram of vels """
        
        phi = 0
        theta = 45
        
        BINS = 200
        
        " Transmission "
        trans = T
        #print np.sum(trans)
    
        " Convert angle to radians "
    
        theta = np.deg2rad(45)
        phi = np.deg2rad(0)
    
        """ Arrange velocities about y-axis so we only require z component"""
    
        partvelx = partvels[:,0] * np.cos(phi) - partvels[:, 1] * np.sin(phi)
        partvely = (partvels[:, 0] * np.sin(phi) + partvels[:, 1] * np.cos(phi)) * np.cos(theta) - partvels[:,2] * np.sin(theta)
        partvelz = partvels[:, 1] * np.sin(theta) + partvels[:, 2] * np.cos(theta)
        
        unitvel = 0.01 * udist / utime
        
        angstroms = partvelz * unitvel / 3e8   
        angstroms += 1
        angstroms *= 6562.8
        
        counts, bins = np.histogram(angstroms, bins=BINS)
        binDiff = bins[1] - bins[0]
        
        " Sort angstroms into the bins manually and scale based on emissivity and transmissivity "
        countsPlot = np.zeros_like(counts)
        for i in range(len(angstroms)):
            binNum = (angstroms[i] - bins[0]) / binDiff
            #print(int(binNum))
            if(R[i]<100):
                if(int(binNum) < BINS-1):
                    countsPlot[int(binNum)] += trans[i] #* emissScaling[i]
                else:
                    countsPlot[BINS-1] += trans[i] #* emissScaling[i]
                
        countsPlot = countsPlot# / float(np.sum(countsPlot))
        
        axarr[1].plot(bins[0:-1], countsPlot)#, label=name)
    
        
#        axarr[1].plot(bins[0:-1], counts)
#        axarr[1].hist(angstroms, bins=100)
        axarr[1].set_xlabel('Angstroms')
#        axarr[1].set_ylabel('Counts')
        axarr[1].set_ylim(0, 0.015)
        #axarr[1].set_xlim(6556, 6570)
        axarr[1].set_xlim(6545, 6580)
        axarr[1].set_yticks([])
        axarr[1].legend()

        plt.gca().set_aspect('equal', adjustable='box')
    
        if(show == 'y'):
            plt.show()

    elif(plotType == 'HalphaAngleEffect'):

        phi = np.deg2rad(orientation[0])
        theta = np.deg2rad(orientation[1])

        theta1 = orientation[1]
        phi1 = orientation[0]

        """ Density Scatter """
        
        f, axarr = plt.subplots(1, 2, figsize=(16, 8))

        """ Image """

        axarr[0].plot(np.cos(np.deg2rad(range(0, 360, 10))), np.sin(np.deg2rad(range(0, 360, 10))), color='black')
        axarr[0].plot([1, 25, 25, 1], [0, 0.03*25, -0.03*25, 0], color='black')
        axarr[0].plot([-1, -25, -25, -1], [0, 0.03*25, -0.03*25, 0], color='black')
        axarr[0].plot([0, 0], [0, 25], '-', color='black')
        #arrowy = np.array([3.5, 1.5, 2, 2, 1.5])
        #arrowx = np.array([0, 0, 0.5, 0.5, 0])
        arrowy = np.array([0, 25])
        arrowx = np.array([0, 0])
        axarr[0].plot(-arrowx * np.cos(theta) + arrowy * np.sin(theta),
                      arrowx * np.sin(theta) + arrowy * np.cos(theta), color='red' )

        axarr[0].text(-arrowx[1] * np.cos(theta) + arrowy[1] * np.sin(theta),
                      arrowx[1] * np.sin(theta) + arrowy[1] * np.cos(theta), "{} degrees".format(theta1) )

        axarr[0].set_xlim(-26, 26)
        axarr[0].set_ylim(-26, 26)
        axarr[0].axis("off")
        #axarr[0].set_xlabel("R / stellar radii")
        #axarr[0].set_ylabel("R / stellar radii")


        """ Histogram of vels """
        
        phi = 0
        #theta = 45
        
        BINS = 200
        
        """ Arrange velocities about y-axis so we only require z component"""
    
        partvelx = partvels[:,0] * np.cos(phi) - partvels[:, 1] * np.sin(phi)
        partvely = (partvels[:, 0] * np.sin(phi) + partvels[:, 1] * np.cos(phi)) * np.cos(theta) - partvels[:,2] * np.sin(theta)
        partvelz = partvels[:, 1] * np.sin(theta) + partvels[:, 2] * np.cos(theta)
        
        unitvel = 0.01 * udist / utime
        
        angstroms = partvelz * unitvel / 3e8   
        angstroms += 1
        angstroms *= 6562.8
        
        counts, bins = np.histogram(angstroms, bins=BINS)
        binDiff = bins[1] - bins[0]
        
        " Sort angstroms into the bins manually and scale based on emissivity and transmissivity "
        countsPlot = np.zeros_like(counts)
        for i in range(len(angstroms)):
            binNum = (angstroms[i] - bins[0]) / binDiff
            #print(int(binNum))
            if(int(binNum) < BINS-1):
                countsPlot[int(binNum)] += 1
            else:
                countsPlot[BINS-1] += 1
                
        countsPlot = countsPlot / float(np.sum(countsPlot))
        
        axarr[1].plot(bins[0:-1], countsPlot)#, label=name)
    
        
#        axarr[1].plot(bins[0:-1], counts)
#        axarr[1].hist(angstroms, bins=100)
        axarr[1].set_xlabel('Angstroms')
#        axarr[1].set_ylabel('Counts')
        axarr[1].set_ylim(0, 0.015)
        #axarr[1].set_xlim(6556, 6570)
        axarr[1].set_xlim(6545, 6580)
        axarr[1].set_yticks([])
        axarr[1].legend()

        if(show == 'y'):
            plt.show()
            
    elif(plotType == 'ODAngleEffect'):

        phi = np.deg2rad(orientation[0])
        theta = np.deg2rad(orientation[1])

        theta1 = orientation[1]
        phi1 = orientation[0]

        """ Density Scatter """
        
        f, axarr = plt.subplots(1, 2, figsize=(22, 8))

        """ Image """

        axarr[0].plot(np.cos(np.deg2rad(range(0, 360, 10))), np.sin(np.deg2rad(range(0, 360, 10))), color='black')
        axarr[0].plot([1, 25, 25, 1], [0, 0.03*25, -0.03*25, 0], color='black')
        axarr[0].plot([-1, -25, -25, -1], [0, 0.03*25, -0.03*25, 0], color='black')
        axarr[0].plot([0, 0], [0, 25], '-', color='black')
        #arrowy = np.array([3.5, 1.5, 2, 2, 1.5])
        #arrowx = np.array([0, 0, 0.5, 0.5, 0])
        arrowy = np.array([0, 25])
        arrowx = np.array([0, 0])
        axarr[0].plot(-arrowx * np.cos(theta) + arrowy * np.sin(theta),
                      arrowx * np.sin(theta) + arrowy * np.cos(theta), color='red' )

        axarr[0].text(-arrowx[1] * np.cos(theta) + arrowy[1] * np.sin(theta),
                      arrowx[1] * np.sin(theta) + arrowy[1] * np.cos(theta), "{} degrees".format(theta1) )

        axarr[0].set_xlim(-26, 26)
        axarr[0].set_ylim(-26, 26)
        axarr[0].axis("off")
        #axarr[0].set_xlabel("R / stellar radii")
        #axarr[0].set_ylabel("R / stellar radii")


        """ Density Scatter """
        partx = np.zeros(len(partcoords[:,0]))
        party = np.zeros(len(partcoords[:,0]))
        partz = np.zeros(len(partcoords[:,0]))
        NSx = NScoords[0]
        NSy = NScoords[1]
        NSz = NScoords[2]
        
        phi = np.deg2rad(orientation[0])
        theta = np.deg2rad(orientation[1])
        
#        " Rotate in z axis "
#        for i in range(len(partcoords)):
#            partx[i] = xcoords[i] * np.cos(phi) - ycoords[i] * np.sin(phi)
#            party[i] = xcoords[i] * np.sin(phi) + ycoords[i] * np.cos(phi)
#        
#        " Rotate in y axis "
#        for i in range(len(partcoords)):
#            party[i] = partcoords[i,1] * np.cos(theta) - partcoords[i,2] * np.sin(theta)
#            partz[i] = partcoords[i,1] * np.sin(theta) + partcoords[i,2] * np.cos(theta)

        
        NSx = NScoords[0] * np.cos(phi) - NScoords[1] * np.sin(phi)
        NSy = NScoords[0] * np.sin(phi) + NScoords[1] * np.cos(phi) 
        NSy = NScoords[1] * np.cos(theta) - NScoords[2] * np.sin(theta)
        NSz = NScoords[1] * np.sin(theta) + NScoords[2] * np.cos(theta)

        "rotate arrays"        
        
        xcoords = partcoords[:,0] 
        ycoords = partcoords[:,1] 
        zcoords = partcoords[:,2]
        
#        """ sort arrays """
#        
#        sortOrder = np.argsort(zcoords)[::-1]
#        
#        xcoords = xcoords[sortOrder]
#        ycoords = ycoords[sortOrder]
#        zcoords = zcoords[sortOrder]
        
        "delete low emission particles"
        
#        delArr = np.array([])
#        for i in range(len(T)):
#            if(T[i] <= 1e-3):
#            if(R[i] > 10):
#                delArr = np.append(delArr, i)
                
#        xcoords = np.delete(xcoords, delArr)
#        ycoords = np.delete(ycoords, delArr)
#        zcoords = np.delete(zcoords, delArr)
#        T = np.delete(T, delArr)
        
#        print len(delArr), len(xcoords)
        
        "reverse order of all arrays"
        xcoords = xcoords[-1::-1]
        ycoords = ycoords[-1::-1]
        zcoords = zcoords[-1::-1]
        T = T[-1::-1]
#        print len(zcoords)
        
#        plt.scatter(xcoords * np.cos(phi) - ycoords * np.sin(phi), 
#                    (xcoords * np.sin(phi) + ycoords * np.cos(phi))*np.cos(theta) - zcoords * np.sin(theta), 
#                    color=colours, s = 0.25)
        sc = axarr[1].scatter(xcoords, ycoords, 
                    c=T, s = 2, vmin=0, vmax=1, cmap=cm, edgecolor='')
        cb = plt.colorbar(sc, ax=axarr[1])
        cb.ax.set_ylabel("Transmission")
#        plt.scatter(Becoords[0], Becoords[1], color='b')
        if(NS == 'y'):
            axarr[1].scatter(NScoords[0], NScoords[1],
                        color='black')
        else:
            axarr[1].scatter(NScoords[0], NScoords[1],
                        color='r', s=0.25)
        #plt.text(-14.0, 14.0, 't = {} days'.format(int(time) - 1856), size=20)
        if(xylims != 'None'):
            axarr[1].set_xlim(xylims[0], xylims[1])
            axarr[1].set_ylim(xylims[2], xylims[3])
        axarr[1].set_xlabel(r"R / Be star radius")
        axarr[1].set_ylabel(r"R / Be star radius")

        if(show == 'y'):
            plt.show()
            
#    elif(plotType == 'EW'):
#        
#        plt.figure(figsize=(18,8))
#    
#        """ Density Scatter """
#        partx = np.zeros(len(partcoords[:,0]))
#        party = np.zeros(len(partcoords[:,0]))
#        partz = np.zeros(len(partcoords[:,0]))
#        NSx = NScoords[0]
#        NSy = NScoords[1]
#        NSz = NScoords[2]
#        
#        phi = np.deg2rad(orientation[0])
#        theta = np.deg2rad(orientation[1])
#        
##        " Rotate in z axis "
##        for i in range(len(partcoords)):
##            partx[i] = xcoords[i] * np.cos(phi) - ycoords[i] * np.sin(phi)
##            party[i] = xcoords[i] * np.sin(phi) + ycoords[i] * np.cos(phi)
##        
##        " Rotate in y axis "
##        for i in range(len(partcoords)):
##            party[i] = partcoords[i,1] * np.cos(theta) - partcoords[i,2] * np.sin(theta)
##            partz[i] = partcoords[i,1] * np.sin(theta) + partcoords[i,2] * np.cos(theta)
#
#        
#        NSx = NScoords[0] * np.cos(phi) - NScoords[1] * np.sin(phi)
#        NSy = NScoords[0] * np.sin(phi) + NScoords[1] * np.cos(phi) 
#        NSy = NScoords[1] * np.cos(theta) - NScoords[2] * np.sin(theta)
#        NSz = NScoords[1] * np.sin(theta) + NScoords[2] * np.cos(theta)
#
##        plt.scatter(xcoords * np.cos(phi) - ycoords * np.sin(phi), 
##                    (xcoords * np.sin(phi) + ycoords * np.cos(phi))*np.cos(theta) - zcoords * np.sin(theta), 
##                    color=colours, s = 0.25)
#        plt.subplot(121)
#        sc = plt.scatter(xcoords * np.cos(phi) - ycoords * np.sin(phi), 
#                    (xcoords * np.sin(phi) + ycoords * np.cos(phi))*np.cos(theta) - zcoords * np.sin(theta), 
#                    c=rho, s = 1, vmin=rhomin, vmax=rhomax, cmap=cm, edgecolor='')
#        cb = plt.colorbar(sc)
#        cb.ax.set_ylabel("Density in g per cm cubed")
#        plt.scatter(NScoords[0] * np.cos(phi) - NScoords[1] * np.sin(phi), 
#                    (NScoords[0] * np.sin(phi) + NScoords[1] * np.cos(phi))*np.cos(theta) - NScoords[2] * np.sin(theta),
#                    color='b')
#                    
#        if(xylims != 'None'):
#            plt.xlim(xylims[0], xylims[1])
#            plt.ylim(xylims[2], xylims[3])
#        plt.xlabel(r"R / Be star radius")
#        plt.ylabel(r"R / Be star radius")
#    
#        """ Histogram of vels """
#            
#        Becoords, NScoords, partcoords, partvels, rho, u, num_parts = readQDISK(filename[0])
#    
#        Becoords /= udist / (7 * 6.96e10)
#        NScoords /= udist / (7 * 6.96e10)
#        partcoords /= udist / (7 * 6.96e10)
#        
#        phi = np.deg2rad(90)
#        theta = np.deg2rad(90)
#        
#        BINS=100
#        
#        partvelx = partvels[:,0] * np.cos(phi) - partvels[:, 1] * np.sin(phi)
#        partvely = (partvels[:, 0] * np.sin(phi) + partvels[:, 1] * np.cos(phi)) * np.cos(theta) - partvels[:,2] * np.sin(theta)
#        partvelz = partvels[:, 1] * np.sin(theta) + partvels[:, 2] * np.cos(theta)
#        
#        unitvel = 0.01 * udist / utime
#        
#        angstroms = partvelz * unitvel / 3e8   
#        angstroms += 1
#        angstroms *= 6562.8
#        
#        counts, bins = np.histogram(angstroms, bins=BINS)
#        binDiff = bins[1] - bins[0]
#        
#        baseName = "injBurst_40d_0.4e"
#        FileName = "/data/rob1g10/SPH/DATA/DISKS/{}/optDepthData{}_theta{}_phi{}_time{:0.1f}.txt".format(baseName, baseName, 90, 90, time)# * utime / (24 * 3600))
#        
#        oD = np.array([])
#        " Open file "
#        f = open(FileName, 'r')
#    
#        for line in f:
#            line = line.strip()
#            oD = np.append(oD, float(line))
#
#        f.close()
#        
#        T = np.exp(-oD) / 1e-3
#        
#        " Sort angstroms into the bins manually and scale based on emissivity and transmissivity "
#        countsPlot = np.zeros_like(counts)
#        for i in range(len(angstroms)):
#            binNum = (angstroms[i] - bins[0]) / binDiff
#            if(int(binNum) < BINS-1):
#                try:
#                    countsPlot[int(binNum)] += T[i] # emissScaling[i]
#                except:
#                    pass
#            else:
#                countsPlot[BINS-1] += T[i] #* emissScaling[i]
#                
#        countsPlot = countsPlot / float(np.sum(countsPlot))
#        
#        plt.subplot(224)
#        plt.plot(bins[0:-1], countsPlot)#, label=name)
#    
#        plt.xlabel('Angstroms')
##        plt.ylabel('Counts')
#        plt.yticks([])
#        plt.ylim(0, 0.03)
#        plt.xlim(6540, 6580)
#        
#        """ EW """
#    
#        plt.subplot(222)
#        # 8 values per QGDISK
#        NSarray = createAccretionArrayForDISKs(['/data/rob1g10/SPH/DATA/DISKS/20d_0.4e/disk004', 
#                                                '/data/rob1g10/SPH/DATA/DISKS/20d_0.4e/disk005'],
#                                                'NS', 265, 285)
#        for i in range():                                       
#        OD = readOptDepthData4py("/data/rob1g10/SPH/DATA/DISKS/{}/optDepthData{}_theta{}_phi{}_time{:0.1f}.txt".format(baseName, baseName, 0, 0, time * utime / (24 * 3600)))
#        trans = (np.exp(-OD)) * 1e3
#        
#        QGDISKNum = int(filename[0][-5:-2])
#        QGstepNum = int(filename[0][-2:])
#        currentElement = (QGDISKNum - 265) * 8 + QGstepNum
#        
#        
#        IPfilenames = np.array([])
#        accreteTimes = np.array([])
#        inputdirectory = "/data/rob1g10/SPH/DATA/Animations/MSSL/processedDISKs/"
#    
#        os.system("ls {}IP* > lsIP.txt".format(inputdirectory))
#        
#        """ IPDISKs """    
#        
#        f = open("lsIP.txt", "r")
#        
#        lines = f.readlines()
#    
#        for line in lines:
#            f2 = open("{}".format(line[0:-1]))
#            temp = f2.readlines()
#            f2.close()
#            if(len(temp) != 0):
#                IPfilenames = np.append(IPfilenames, line[0:-1])
#            
#        f.close() 
#        
#        os.system("rm lsIP.txt")
#        
#        for IP in IPfilenames:
#            accreteTimes = np.append(accreteTimes, readIPDISK(IP))
#        accreteTimes *= utime/(24*3600)
#        
#        plt.plot(accreteTimes[:currentElement] - int(min(accreteTimes)), NSarray[:currentElement])  
#        plt.ylim(0, 102 * 1e-15 * 0.01 * umass / 2e30)
#        plt.xlim(0, int(max(accreteTimes)) - int(min(accreteTimes)))
#        plt.xlabel("Days")
#        plt.ylabel("NS Accretion in Solar Masses")
##        print accreteTimes[QGstepNum:QGstepNum+8]
##        print NSarray[currentElement-8:currentElement]                                                 
#        
#    
#        if(show == 'y'):
#            plt.show()

    else:
        print('ERROR: plotType must be "plusHist", "3d" or "2d".')

    """ Save fig """

    if(save == 'y'):
        if(plotType == 'NSaccretion'):
            plt.savefig(outFilename, figsize=(36, 16), dpi=100, bbox_inches='tight')
        else:
            plt.savefig(outFilename, figsize=(18, 16), dpi=100)

    if(show == 'n'):
        plt.close()
    
    return partcoords[:,1] * np.sin(theta) + partcoords[:,2] * np.cos(theta)


def plotRotatingFigs(filenames, pngdirectory, baseName):
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
    
    for theta in range(0, 181, 15):
        for phi in range(0, 181, 15):
            for zoom in [20, 10, 5]:
                Becoords, NScoords, partcoords, partvels, rho, u, num_parts, time, umass, udist, utime = readFILES(filenames[0], filenames[1], filenames[2])
                name = "{}{}_theta{}_phi{}_zoom{}_time{:0.1f}.png".format(pngdirectory, baseName, theta, phi, zoom, time * utime / (24 * 3600))
                plotAndSaveFig([filenames[0], filenames[1]] , time * utime / (24 * 3600), name, '2d', orientation=[0 + phi, 0 + theta], show='n', save='y', xylims=[-zoom, zoom, -zoom, zoom])

    return

def createAnimation(animationName, pngdirectory, inputdirectory, animType):
    """
    Produces an animation using pngs created and placed in a given folder.

    GEOMETRY
    --------

    z axis is height of disk, y axis is width and x axis is depth.

    INPUT VARIABLES
    ---------------

    animationName (string) - what the animation will be named.

    pngdirectory (string) - directory where png files are found.

    directory (string) - directory where the QGDISK and IPDISK files reside.
    
    animtype (string) - type to pass to plotAndSaveFig.
    
    utime, orbitPeriod (floats) - unit of time and orbit period for simulation.

    OUTPUT VARIABLES
    ----------------

    An animation of the respective files.

    """
    
    QGfilenames = np.array([])
    IPfilenames = np.array([])
    
    os.system("ls {}QG* > lsQG.txt".format(inputdirectory))
    os.system("ls {}IP* > lsIP.txt".format(inputdirectory))
    
    umass, udist, utime = readUnits("{}disk001".format(inputdirectory))
    
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
    
    """ IPDISKs """
    
    f = open("lsIP.txt", "r")
    
    lines = f.readlines()
    
    for line in lines:
        IPfilenames = np.append(IPfilenames, line[0:-1])

    f.close()

    time = np.array([])
    
    for i in range(len(IPfilenames)):
        time = np.append(time, readIPDISK(IPfilenames[i]))
    
    """ Read every file and create an image for it """

    for i in range(len(QGfilenames)):
        if(i <= 8):
            name = "{}AniPart000{}".format(pngdirectory, i+1)
        if(i > 8 and i <= 98):
            name = "{}AniPart00{}".format(pngdirectory, i+1)
        if(i > 98 and i < 1000):
            name = "{}AniPart0{}".format(pngdirectory, i+1)
        if(i >= 1000 and i < 10000):
            name = "{}AniPart{}".format(pngdirectory, i+1)
        print "Making image", i+1, "of", len(QGfilenames)
        if(i == 0):
            Becoords, NScoords, partcoords, partvels, rho, u, num_parts = readQDISK(QGfilenames[i])
            MINRHO = min(rho) * umass / (udist ** 3)
            MAXRHO = max(rho) * umass / (udist ** 3)

        plotAndSaveFig([QGfilenames[i], "{}disk001".format(inputdirectory), "{}IPDISK{}".format(inputdirectory, QGfilenames[i][-5:-2])], 
                        time[i]*utime/(24*3600), name, animType,
                       minRho=MINRHO, maxRho=MAXRHO, orientation=[0,0], show='n', save='y', xylims=[-20, 30, -16, 16])

    """ Create animation """

    os.system("convert -delay 20 -loop 0 {}*.png {}.gif".format(pngdirectory, animationName))
    
    os.system("rm lsQG.txt lsIP.txt".format(inputdirectory))

def plotDiskProfile(partcoords, rho, tol=1e-3):
    """
    Produces a scatter plot of the disk profile for the Be star.

    GEOMETRY
    --------

    z axis is height of disk, y axis is width and x axis is depth.

    INPUT VARIABLES
    ---------------

    partcoords ((npart, 3) float array) - x,y,z coordinates of the disk
                                         particles.

    rhos (npart float array) - density values of the disk particles.

    tol (float) - tolerance either side of x=0 that the function will use particles for
          plotting.

    OUTPUT VARIABLES
    ----------------

    A plot of the disk profile.

    """

    " Initialise arrays "

    ycoordsPLOT = np.array([])
    zcoordsPLOT = np.array([])

    " Get only particles in one slice if the disk i.e. x=0 "

    for i in range(len(partcoords[:,0])):
        if(np.abs(partcoords[i,0]) <= tol):
            ycoordsPLOT = np.append(ycoordsPLOT, partcoords[i,1])
            zcoordsPLOT = np.append(zcoordsPLOT, partcoords[i,2])

    rhomax = max(rho)
    colours = np.zeros((len(rho), 3))
    for i in range(len(rho)):
        rgrad = round(1 - rho[i]/rhomax, 1)
        if(rgrad > 1):
            rgrad = 1
        colours[i] = [rgrad, 0, 0]

    """ Density Scatter """

    plt.scatter(ycoordsPLOT, zcoordsPLOT, color=colours, s=0.5)
    plt.xlim(-0.4, 0.4)
    plt.ylim(-0.1, 0.1)
    plt.show()

"""----------------------------------------------------------------------------

                FINDING PARTICLES

----------------------------------------------------------------------------"""

def findParts(partcoords, boxSize):
    """
    Counts the number of particles inside the given box and returns their
    elements in the particle array.

    GEOMETRY

    z axis is height of disk, y axis is width and x axis is depth.


                         z |
                           |
                           |
                           |_ _ _ _ y
                          /
                         /
                        / x

    INPUT VARIABLES

    partcoords ((npart, 3) float array) - x,y,z coordinates of the disk
                                         particles

    boxSize (6 element float array) - the dimensions of the box that the
                                      program will look for particles within.
                                      xmin, xmax, ymin, ymax, zmin, zmax

    OUTPUT VARIABLES

    elements (array of floats) - the elements of the particles found in the
                                 box.

    """

    elements = np.array([])

    for i in range(len(partcoords)):
        if(partcoords[i, 0] > boxSize[0] and partcoords[i, 0] < boxSize[1]
            and partcoords[i, 1] > boxSize[2] and partcoords[i, 1] < boxSize[3]
            and partcoords[i, 2] > boxSize[4] and partcoords[i, 2] < boxSize[5]
            ):

            elements = np.append(elements, i)

    return elements


"""----------------------------------------------------------------------------

                Halpha Stuff

----------------------------------------------------------------------------"""

def countParts(gridBL, partcoords, sizeOfView, depthOfView):
    """
    Counts the number of particles that an observer will see existing in
    the disk of the Be star when looking at any chosen angle. It then returns
    the array elements of the particle array which are contained in this view.

    GEOMETRY

    z axis is height of disk, y axis is width and x axis is depth.


                         z |
                           |
                           |
                           |_ _ _ _ y
                          /
                         /
                        / x

    INPUT VARIABLES

    gridBL (2 float array) - y,z coords of the bottom left of grid box

    partcoords ((npart, 3) float array) - x,y,z coordinates of the disk
                                         particles

    sizeOfView (2 element float array) - defines the size of the square which
                                        the observer is counting particles in.
                                        (y and z lengths)

    depthOfView (integer) - how far into the disk that can be seen optically
                          (as a number of particles).

    OUTPUT VARIABLES

    elements (len 2 integer array) - elements of particle array that are in
                                     the view of the observer and the number
                                     of particles counted for this element.

    """

    " Replace input array so that it is not altered "

    partcoords_func = np.array(partcoords)

    " Centre on gridTL + 0.5 sizeOfView "

    origin = gridBL + 0.5*sizeOfView
    for i in range(len(partcoords_func)):
        partcoords_func[i][1] -= origin[0]
        partcoords_func[i][2] -= origin[1]

    """ Calculate how far either side of the centre of the grid point the
    viewer sees """

    zlim = 0.5*sizeOfView[1]
    ylim = 0.5*sizeOfView[0]

    " Count particles up to depthOfView within square "

    elements = np.array([])
    collect = np.array([])
    for i in range(len(partcoords_func)):

        " Collect particles within grid "

        if (np.abs(partcoords_func[i][2]) < zlim
             and np.abs(partcoords_func[i][1]) < ylim):
                 collect = np.append(collect, [partcoords_func[i][0], i])
                 #print collect

    """ Split collect into two arrays, one with particle x coordinate and
    and the other with the element it is found in """

    collectCoords = collect[0::2]
    collectElements = collect[1::2]

    """ Find particles with max x coordinate up to depthOfView particles and
    collect their elements """

    j = 0
    while(j < depthOfView):

        " Check collect is not empty "

        if(len(collectCoords) == 0):
            return elements

        " Find max element and register it "
        maxElement = np.argmax(collectCoords)

        elements = np.append(elements, collectElements[maxElement])
        j += 1
        collectElements = np.delete(collectElements, maxElement)
        collectCoords = np.delete(collectCoords, maxElement)

    #print elements
    return elements


def binVels(partcoords, pivot, partvelocities, depthOfView, viewAngle,
            numOfGridPoints, diskMaxY, diskMaxZ):
    """
    Bins velocities of particles into a histogram such that a spectra can be
    made.

    GEOMETRY
    --------

    z axis is height of disk, y axis is width and x axis is depth.

                         z |
                           |
                           |
                           |_ _ _ _ y
                          /
                         /
                        / x

    INPUT VARIABLES
    ---------------

    partcoords ((npart, 3) float array) - x,y,z coordinates of the disk
                                         particles

    pivot (3 float array) - x,y,z coords of the origin to rotate about

    partvels ((npart, 3) float array) - array containing all the x,y,z
                                        velocities of the disk particles.

    depthOfView (integer) - how far into the disk that can be seen optically
                          (as a number of particles).

    viewAngle (degrees) - the angle at which the observer is viewing the disk
                        (90 degrees is considered to be directly above the star)

    numOfGridPointsX (2 integer array) - defines the number of squares which
                                         the disk is split into
                                         (y and z directions).

    diskMaxY, diskMaxZ - length from origin that the disk extends to in y and
                         z coordinates respectively in computational units
                         (plot the system then estimate).

    OUTPUT VARIABLES
    ----------------

    A array of velocities for plotting.

    """

    " Replace input array so that it is not altered "

    partvels_func = partvelocities[:,0]
    partvels_funcPlot = np.array([])

    partcoords_func = partcoords

    " Convert angle to radians "

    viewAngle = np.deg2rad(viewAngle)

    """ Arrange velocities about y-axis so we only require x component"""

    for i in range(len(partvels_func)):
        partvelx = partvelocities[i][0]
        partvelz = partvelocities[i][2]
        partvels_func[i] = partvelx * np.cos(viewAngle) + partvelz * np.sin(viewAngle)

    """ Arrange coordinates such that pivot is at the origin and then
    rotate coordinates about y-axis """

    for i in range(len(partcoords_func)):
        partcoords_func[i] -= pivot
        partx = partcoords_func[i][0]
        partz = partcoords_func[i][2]
        partcoords_func[i][0] = partx * np.cos(viewAngle) + partz * np.sin(viewAngle)
        partcoords_func[i][2] = - partx * np.sin(viewAngle) + partz * np.cos(viewAngle)

    """ Calculate size of disk at observer angle so it can be split up into grid
    cells,
     make sure to cut off particles outside of disk! """

    #diskMaxY = 0.7
    #diskMaxZ = 0.5


    " Make grid squares of size diskMaxY/grid points by diskMaxZ/grid points "

    sizeOfSquare = np.array([2 * diskMaxY / numOfGridPoints[0],
                             2 * diskMaxZ / numOfGridPoints[1]
                             ])

    " Use count parts on a grid splitting up the disk "
    elements = np.array([])
    for i in range(numOfGridPoints[0]):
        for j in range(numOfGridPoints[1]):
            #print i, j
            " Set input up for countParts "

            gridBL = np.array([-diskMaxY + (i * sizeOfSquare[0]),
                               -diskMaxZ + (j * sizeOfSquare[1])
                               ])

            " Use countParts "

            elements = np.append(elements, countParts(gridBL, partcoords_func,
                                           sizeOfSquare, depthOfView))

    for i in range(len(elements)):
        partvels_funcPlot = np.append(partvels_funcPlot,
                                      partvels_func[elements[i]])

    #plt.hist(partvels_funcPlot, weights=rhos, bins=200)
    #plt.hist(partvels_funcPlot, bins=200)
    #plt.xlim(-0.05, 0.05)
    #plt.show()

    return partvels_funcPlot


def showPeakRelation4Spectra(partcoords, pivot, partvelocities, rhos,
                             depthOfView, numOfGridPoints):
    """
    Uses binVels to try to show a cos(theta) relation for the movement of
    the peaks in the produced spectra.

    !! BROKEN !!

    GEOMETRY
    --------

    z axis is height of disk, y axis is width and x axis is depth.

                         z |
                           |
                           |
                           |_ _ _ _ y
                          /
                         /
                        / x

    INPUT VARIABLES
    ---------------

    partcoords ((npart, 3) float array) - x,y,z coordinates of the disk
                                         particles

    pivot (3 float array) - x,y,z coords of the origin to rotate about

    partvels ((npart, 3) float array) - array containing all the x,y,z
                                        velocities of the disk particles.

    rhos (npart float array) - array containing all the densities of the disk
                               particles.

    depthOfView (integer) - how far into the disk that can be seen optically
                          (as a number of particles).

    numOfGridPointsX (2 integer array) - defines the number of squares which
                                         the disk is split into
                                         (y and z directions).

    OUTPUT VARIABLES
    ----------------

    A plot of the movement of the peaks with viewing angle.

    """

    " Create array for the plot at the end "

    plotArray = np.array([])
    partvels_funcPLOT = np.array([])
    weights_func = np.array([])

    " Loop over angles and find the peaks "
    for i in range(20, 91, 5):
        partvels_funcPLOT, weights_func = binVels(partcoords, pivot, partvelocities,
                                        rhos, depthOfView, i,
                                        numOfGridPoints)
        counts, bins = np.histogram(partvels_funcPLOT, weights=weights_func,
                                      bins=1000)
        print counts
        plotArray = np.append(plotArray, 2*bins[np.argmax(counts)])
        print i
        del counts
        del bins
        del partvels_funcPLOT
        del weights_func

    fig2 = plt.figure("fig2")
    fig2.plot(plotArray, range(0, 91))
    #plt.show(block=False)
    fig2.show()
    
def findDiskSize(filenames):
    SIZES = np.array([])
    for filename in filenames:
        diskname = filename[:-26] + "disk001"
        IPname = filename[:-11] + "IPDISK" + filename[-5:-2]
        Becoords, NScoords, partcoords, partvels, rho, u, num_parts, time, umass, udist, utime = readFILES(filename, diskname, IPname)

        rhoCutOff = 1e-13 #max(rho) * 1e-4
        R = np.sqrt(partcoords[:,0]**2 + partcoords[:,1]**2)
        sortOrder = np.argsort(R)
        R = R[sortOrder]
        rho = rho[sortOrder]
        Rad = max(R)
        for i in range(len(R)):
            if(i==len(R)-1):
                Rad = R[i]
#                print filename[-42:-26], Rad * udist / (7 * 6.96e10)
                break
            if(rho[i] < rhoCutOff):
                Rad = R[i]
#                print filename[-42:-26], Rad * udist / (7 * 6.96e10)
                break
        SIZES = np.append(SIZES, Rad * udist / (7 * 6.96e10))
    return SIZES
    
#def findDiskSize(baseName, DISKstart, DISKend, DISKsteps, STEPstart, STEPend,tol=4e-4):
#    discSizes = np.array([])
#    for i in range(DISKstart, DISKend+1):
#        if(i == DISKstart):
#            start = STEPstart
#        else:
#            start = 1
#        if(i == DISKend):
#            end = STEPend
#        else:s
#            end = DISKsteps
#        for j in range(start, end):
#            filename = '/data/rob1g10/SPH/DATA/DISKS/{}/processedDISKs/QGDISK{:03}{:02}'.format(baseName, i, j)
##            print filename, start, end, i, j
#            diskname = filename[:-26] + "disk001"
#            IPname = filename[:-11] + "IPDISK" + filename[-5:-2]
#            Becoords, NScoords, partcoords, partvels, rho, u, num_parts, time, umass, udist, utime = readFILES(filename, diskname, IPname)
#            #COMPLICATED WAY            
##            R = np.sqrt(partcoords[:,0]**2 + partcoords[:,1]**2)
##            R = np.sort(R)
##            Rad = 0
##            for k in range(len(R) - 1):
##                check = np.abs(R[k] - R[k+1])
##                if(check > tol):
##                    Rad = R[k]
##    #                print filename[-42:-26], Rad * udist / (7 * 6.96e10)
##                    break
#            rhoCutOff = max(rho) * 1e-2
#            R = np.sqrt(partcoords[:,0]**2 + partcoords[:,1]**2)
#            sortOrder = np.argsort(R)
#            R = R[sortOrder]
#            rho = rho[sortOrder]
#            Rad = max(R)
#            for i in range(len(R)):
#                if(rho[i] < rhoCutOff):
#                    Rad = R[i]
#                    print filename[-42:-26], Rad * udist / (7 * 6.96e10)
#                    break
#            discSizes = np.append(discSizes, Rad * udist / (7 * 6.96e10))
#    print( max(discSizes), min(discSizes), np.mean(discSizes), np.std(discSizes), len(discSizes))
#    return Rad


"""----------------------------------------------------------------------------

                PERIASTRON DISTANCE

----------------------------------------------------------------------------"""

def calcSemiMajorFromPeriod(P, binaryMass):
    """
    Calculates the radius in units of the semimajor axis.
    
    INPUT VARIABLES
    ---------------
    
    P (float) - period of the orbit in days.

    binaryMass (float) - mass of the binary star system in solar masses.
    
    OUTPUT VARIABLES
    ----------------
    
    period (float) - period in units of days.
    
    """
    
    """ Convert period to seconds, actualRadius to solar radii 
    and binaryMass to kg"""
    
    binaryMass *= 2e30
    P *= (24 * 3600)
    
    " Calculate period "
    
    coefficient = 4 * np.pi * np.pi / (6.67e-11 * binaryMass) 
    
    return ((P ** 2.) / coefficient) ** (1./3)

def calcPeriastronDistance(P, binaryMass, eccentricity):
    a = calcSemiMajorFromPeriod(P, binaryMass)
    perDist = a * (1. - eccentricity)
    return perDist

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
#    print "Be coords:", Becoords
#    partcoords -= Becoords
#    Becoords -= Becoords
#    NScoords -= Becoords
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
    
#    print "Be coords:", np.sqrt(Becoords[0]**2 + Becoords[1]**2 + Becoords[2]**2) * udist / (7e10*7)
#    print "NS coords:", np.sqrt(NScoords[0]**2 + NScoords[1]**2 + NScoords[2]**2) * udist / (7e10*7)
    partcoords -= Becoords
    NScoords -= Becoords
    Becoords -= Becoords
#    print "Be coords:", np.sqrt(Becoords[0]**2 + Becoords[1]**2 + Becoords[2]**2) * udist / (7e10*7)
#    print "NS coords:", np.sqrt(NScoords[0]**2 + NScoords[1]**2 + NScoords[2]**2) * udist / (7e10*7)
    return Becoords, NScoords, partcoords, partvels, rho, u, num_parts, time, umass, udist, utime
#    except:
#        return Becoords, NScoords, partcoords, partvels, rho, u, num_parts, umass, udist, utime
    
    
    

"""----------------------------------------------------------------------------

                RUNNING FUNCTIONS

----------------------------------------------------------------------------"""

"""----------------------------------------------------------------------------

               Mac

----------------------------------------------------------------------------"""

#filename = '../../Documents/PhD/Okazaki SPH/RadialTempSPH/images/QGDISK04002'
#diskname = '../../Documents/PhD/Okazaki SPH/RadialTempSPH/images/disk001'
#IPname = '../../Documents/PhD/Okazaki SPH/RadialTempSPH/images/IPDISK040'

"""----------------------------------------------------------------------------

                PC

----------------------------------------------------------------------------"""

#filename = 'QGDISK09101'
#filename = '/data/rob1g10/SPH/OriginalOkazaki/images/imagep3DSP/MyRead/infoPartsDISK001'
#filename = '/data/rob1g10/SPH/RoryTest2/sph3D/RBDISK1.31000E+02'
#filename = '/data/rob1g10/SPH/RoryTest2/images/imagep3DSP/QGDISK09001'
#filename = '/data/rob1g10/SPH/RoryTest/sph3D/infoPartsDISK  1'
#filename = '/data/rob1g10/SPH/DATA/NSorbit_DISKsize_rel/a1.5e11_isotherm/processedDISKs/QGDISK10925'
#filename = '/data/rob1g10/SPH/DATA/DISKS/justBe_halfInj2/processedDISKs/QGDISK20101'
#filename = '/data/rob1g10/SPH/DATA/DISKS/justBe_quarterInj/processedDISKs/QGDISK59301'
#filename = '/data/rob1g10/SPH/DATA/DISKS/20d_0.4e/processedDISKs/QGDISK28001'
#diskname = '/data/rob1g10/SPH/DATA/DISKS/20d_0.4e/disk001'
#IPname = '/data/rob1g10/SPH/DATA/DISKS/20d_0.4e/processedDISKs/IPDISK280'
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
#filename = '/data/rob1g10/SPH/DATA/DISKS/sph3D_f90_radialTemp/processedDISKs/QGDISK04501'
#diskname = '/data/rob1g10/SPH/DATA/DISKS/sph3D_f90_radialTemp/disk001'
#IPname = '/data/rob1g10/SPH/DATA/DISKS/sph3D_f90_radialTemp/processedDISKs/IPDISK045'
#filename = '/data/rob1g10/SPH/DATA/DISKS/20d_0.4e_eqTest/processedDISKs/QGDISK38722'
#filename = '/data/rob1g10/SPH/DATA/DISKS/20d_0.4e_eqTest/processedDISKs/QGDISK38756'
#filename = '/data/rob1g10/SPH/DATA/DISKS/20d_0.0e_eqTest/processedDISKs/QGDISK38720'
#filename = '/data/rob1g10/SPH/DATA/DISKS/20d_0.2e_eqTest/processedDISKs/QGDISK38725'
#filename = '/data/rob1g10/SPH/DATA/DISKS/80d_0.4e_eqTest/processedDISKs/QGDISK35701'
#filename = '/data/rob1g10/SPH/DATA/DISKS/40d_0.4e_eqTest/processedDISKs/QGDISK38703'
#filename = '/data/rob1g10/SPH/DATA/DISKS/40d_0.4e_eqTest/processedDISKs/QGDISK38907'
#filename = '/data/rob1g10/SPH/DATA/DISKS/isoBe_fromBeginningTenth/processedDISKs/QGDISK41101'
#filename = raw_input('Type the name of the ASCII file to be read')

#filename = '/data/rob1g10/SPH/DATA/DISKS/60d_0.0e/processedDISKs/QGDISK08901'

#filename = '/data/rob1g10/SPH/DATA/DISKS/20d_0.0e/processedDISKs/QGDISK28604'

#filename = '/data/rob1g10/SPH/DATA/DISKS/40d_0.0e/processedDISKs/QGDISK04901'

#filename = '/data/rob1g10/SPH/DATA/DISKS/60d_0.0e/processedDISKs/QGDISK03302'

#filename = '/data/rob1g10/SPH/DATA/DISKS/80d_0.0e/processedDISKs/QGDISK06101'

#filename = '/data/rob1g10/SPH/DATA/DISKS/40d_0.2e_eqTest/processedDISKs/QGDISK38810'


#filename = '/data/rob1g10/SPH/DATA/DISKS/20d_0.4e_eqTest/processedDISKs/QGDISK38701'
#filename = '/data/rob1g10/SPH/DATA/DISKS/20d_0.2e_eqTest/processedDISKs/QGDISK38701'
#filename = '/data/rob1g10/SPH/DATA/DISKS/20d_0.0e_eqTest/processedDISKs/QGDISK38701'
#filename = '/data/rob1g10/SPH/DATA/DISKS/40d_0.4e_eqTest/processedDISKs/QGDISK38701'
#filename = '/data/rob1g10/SPH/DATA/DISKS/40d_0.2e_eqTest/processedDISKs/QGDISK38701'
#filename = '/data/rob1g10/SPH/DATA/DISKS/40d_0.0e_eqTest/processedDISKs/QGDISK38701'

#filename = '/data/rob1g10/SPH/DATA/DISKS/20d_0.2e_eqTest/processedDISKs/QGDISK38731'
#filename = '/data/rob1g10/SPH/DATA/DISKS/BH24d_0.2e/processedDISKs/QGDISK43106'
#filename = '/data/rob1g10/SPH/DATA/DISKS/BH28d_0.2e/processedDISKs/QGDISK43112'
#filename = '/data/rob1g10/SPH/DATA/DISKS/BH20d_0.292e/processedDISKs/QGDISK13018'
#filename = '/data/rob1g10/SPH/DATA/DISKS/BH20d_0.36e/processedDISKs/QGDISK13018'

#filename = '/data/rob1g10/SPH/DATA/DISKS/PSR_J2032+4127/processedDISKs/QGDISK27101'

#filename = ['/data/rob1g10/SPH/DATA/DISKS/20d_0.0e_eqTest/processedDISKs/QGDISK39024',
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
#              '/data/rob1g10/SPH/DATA/DISKS/100d_0.4e_eqTest/processedDISKs/QGDISK38707'][-1]
#filename = '/data/rob1g10/SPH/DATA/DISKS/20d_0.4e_eqTest/processedDISKs/QGDISK39030'
#filename = '/data/rob1g10/SPH/DATA/DISKS/40d_0.4e_eqTest/processedDISKs/QGDISK38907'
#filename = '/data/rob1g10/SPH/DATA/DISKS/80d_0.4e_eqTest/processedDISKs/QGDISK36001'

#filename = '/data/rob1g10/SPH/DATA/DISKS/20d_0.9e/processedDISKs/QGDISK05603'

#filename = '/data/rob1g10/SPH/DATA/DISKS/M1.4_20d_0.0e_5000/processedDISKs/QGDISK12505'
#filename = '/data/rob1g10/SPH/DATA/DISKS/M10_23.7d_0.0e/processedDISKs/QGDISK15517'
#filename = '/data/rob1g10/SPH/DATA/DISKS/20d_0.2e_eqTest/processedDISKs/QGDISK09110'
#filename = '/data/rob1g10/SPH/DATA/DISKS/BH24d_0.2e/processedDISKs/QGDISK09131'
#filename = '/data/rob1g10/SPH/DATA/DISKS/BH24d_0.2e/processedDISKs/QGDISK43106'
#filename = '/data/rob1g10/SPH/DATA/DISKS/BH24.3d_0.2e/processedDISKs/QGDISK09110'
#filename = '/data/rob1g10/SPH/DATA/DISKS/BH25d_0.2e/processedDISKs/QGDISK09128'
#filename = '/data/rob1g10/SPH/DATA/DISKS/BH28d_0.2e/processedDISKs/QGDISK09114'
#filename = '/data/rob1g10/SPH/DATA/DISKS/BH28d_0.2e/processedDISKs/QGDISK43112'
#filename='/data/rob1g10/SPH/DATA/DISKS/M8_35.6d_0.4e/processedDISKs/QGDISK16609'
""" Read File """
#diskname = filename[:-26] + "disk001"
#IPname = filename[:-11] + "IPDISK" + filename[-5:-2]
#filename = 'DATA/QGDISK41101'
#diskname = 'INPUT/disk001'
#IPname = 'INPUT/IPDISK411'
#Becoords, NScoords, partcoords, partvels, rho, u, num_parts, time, umass, udist, utime = readFILES(filename, diskname, IPname)
#print "NS - Be separation:", np.sqrt(NScoords[0]**2 + NScoords[1]**2 + NScoords[2]**2)*udist / (7e10*7)
#plt.scatter(np.sqrt(partcoords[:,0]**2 + partcoords[:,1]**2) * udist / (7 * 7e10), partcoords[:,2] * udist  / (7 * 7e10), s=0.5)
#plt.xlim(1, 2)
#plt.ylim(-0.05, 0.05)
#plt.show()

#Becoords2, NScoords2, partcoords2, partvels2, rho2, u2, num_parts2, time2, umass2, udist2, utime2 = readFILES(filename2)
#print num_parts
#print utime
#print time * utime /(24 * 3600)
#print max(rho * umass / udist**3)
#
#magvels = np.sqrt(partvels[:,0]**2 + partvels[:,1]**2 + partvels[:,2]**2) * (udist / utime)

#Becoords *= udist / (7 * 6.96e10)
#NScoords *= udist / (7 * 6.96e10)
#partcoords *= udist / (7 * 6.96e10)
#Becoords2 *= udist / (7 * 6.96e10)
#NScoords2 *= udist / (7 * 6.96e10)
#partcoords2 *= udist / (7 * 6.96e10)
#R = np.sqrt(partcoords[:,0]**2 + partcoords[:,1]**2) * udist / (7 * 695700e3 * 100)
#print max(R)

#zcoords = plotAndSaveFig([filename, diskname, IPname], time * utime / (24 * 3600), 'outFilename', '2d',
#                   orientation=[0,30], show='y', save='n')
#count = 0                 
#for i in range(len(zcoords)):
#    if(zcoords[i] < 0):
#        count += 1
#print count
#print len(zcoords) - count
#
#print "Apastron distance"
# 1 orbit ---------------------------------------------------------------------
#findDiskSize('20d_0.0e_eqTest', 387, 388, 25, 1, 6)
#findDiskSize('20d_0.2e_eqTest', 387, 387, 35, 1, 33)
#findDiskSize('20d_0.4e_eqTest', 387, 387, 66, 1, 33)
#findDiskSize('40d_0.0e_eqTest', 380, 383, 10, 1, 6)
#findDiskSize('40d_0.2e_eqTest', 380, 383, 11, 1, 3)
#findDiskSize('40d_0.4e_eqTest', 380, 382, 15, 1, 5)
#findDiskSize('60d_0.0e_eqTest', 580, 588, 5, 1, 1)
#findDiskSize('60d_0.2e_eqTest', 580, 586, 5, 1, 3)
#findDiskSize('60d_0.4e_eqTest', 580, 585, 6, 1, 3)
#findDiskSize('80d_0.0e_eqTest', 260, 268, 5, 1, 1)
#findDiskSize('80d_0.2e_eqTest', 260, 268, 5, 1, 1)
#findDiskSize('80d_0.4e_eqTest', 340, 348, 5, 1, 1)
# 60 = [580, 585, 3], [580, 586, 3], [580, 588, 1]
# 80 = [340, 348, 1], [260, 268, 1], [260, 268, 1]
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
DISKstart = np.array([362, 370, 125, 418, 320,
                      300, 290, 350, 290, 360,
                      351, 139,  61, 163, 141])
DISKend = np.array([369, 374, 125, 421, 320,
                    300, 293, 350, 294, 360,
                    354, 141,  62, 166, 142])
STEPstart = np.array([ 5,  4, 13, 32, 14,
                      19, 10, 17, 17, 16,
                       3,  8, 44, 13, 34])
STEPend = np.array([ 2,  2, 45, 16, 47,
                    52, 13, 50,  7, 52,
                     2, 20, 26,  8, 17])
DISKsteps = np.array([ 5, 42, 53, 50, 49, 
                      67, 54, 58, 52, 60,
                      11, 26, 50, 45, 49])
baseName = np.array(['M1.4_20d_0.0e', 'M4_18.4d_0.0e', 'M6_17.4d_0.0e', 'M8_16.6d_0.0e', 'M10_15.8d_0.0e', 
                     'M1.4_28d_0.2e', 'M4_25.7d_0.2e', 'M6_24.3d_0.2e', 'M8_23.1d_0.2e', 'M10_22.1d_0.2e',
                     'M1.4_43d_0.4e', 'M4_39.6d_0.4e', 'M6_37.5d_0.4e', 'M8_35.6d_0.4e', 'M10_34.1d_0.4e'])
for k in range(len(baseName)):
    filenames = np.array([])
    for i in range(DISKstart[k], DISKend[k]+1):
        if(i == DISKstart[k]):
            start = STEPstart[k]
        else:
            start = 1
        if(i == DISKend[k]):
            end = STEPend[k]
        else:
            end = DISKsteps[k]
        for j in range(start, end):
            filenames = np.append(filenames, '/data/rob1g10/SPH/DATA/DISKS/{}/processedDISKs/QGDISK{:03}{:02}'.format(baseName[k], i, j))
    print baseName[k], np.mean(findDiskSize(filenames))
#findDiskSize(['/data/rob1g10/SPH/DATA/DISKS/M1.4_20d_0.0e/processedDISKs/QGDISK36403',
#            '/data/rob1g10/SPH/DATA/DISKS/M6_17.4d_0.0e/processedDISKs/QGDISK12513',
#            '/data/rob1g10/SPH/DATA/DISKS/M10_15.8d_0.0e/processedDISKs/QGDISK32014',
#            
#            '/data/rob1g10/SPH/DATA/DISKS/M1.4_28d_0.2e/processedDISKs/QGDISK30019',
#            '/data/rob1g10/SPH/DATA/DISKS/M6_24.3d_0.2e/processedDISKs/QGDISK35017',
#            '/data/rob1g10/SPH/DATA/DISKS/M10_22.1d_0.2e/processedDISKs/QGDISK36017',
#            
#            '/data/rob1g10/SPH/DATA/DISKS/M1.4_43d_0.4e/processedDISKs/QGDISK35103',
#            '/data/rob1g10/SPH/DATA/DISKS/M6_37.5d_0.4e/processedDISKs/QGDISK06144',
#            '/data/rob1g10/SPH/DATA/DISKS/M10_34.1d_0.4e/processedDISKs/QGDISK14134'])
#print "Periastron Distance"
#
#findDiskSize(['/data/rob1g10/SPH/DATA/DISKS/20d_0.2e_eqTest/processedDISKs/QGDISK39017',
#              '/data/rob1g10/SPH/DATA/DISKS/20d_0.4e_eqTest/processedDISKs/QGDISK39013',
#              '/data/rob1g10/SPH/DATA/DISKS/40d_0.4e_eqTest/processedDISKs/QGDISK38804',
#              '/data/rob1g10/SPH/DATA/DISKS/80d_0.4e_eqTest/processedDISKs/QGDISK35601'])
#              
#print "Mid point"
#
#findDiskSize(['/data/rob1g10/SPH/DATA/DISKS/20d_0.2e_eqTest/processedDISKs/QGDISK39023',
#              '/data/rob1g10/SPH/DATA/DISKS/20d_0.4e_eqTest/processedDISKs/QGDISK39018',
#              '/data/rob1g10/SPH/DATA/DISKS/40d_0.4e_eqTest/processedDISKs/QGDISK38809',
#              '/data/rob1g10/SPH/DATA/DISKS/80d_0.4e_eqTest/processedDISKs/QGDISK35701'])

#print(np.sum(rho))
#plt.show()
#plotDiskProfile(partcoords, rho)
#plotDensityProfileTouhamiandGies(partcoords, rho, u, udist, umass, utime, ylims=['b', '83 day'])
#filename = '/data/rob1g10/SPH/DATA/DISKS/injection_off_JJrun/processedDISKs/QGDISK13605'
#Becoords, NScoords, partcoords, partvels, rho, u, num_parts = readQDISK(filename)
#utime = 2.7502e05
#udist = 5.791e10
#plotDensityProfileTouhamiandGies(partcoords, rho, u, udist, umass, utime, ylims=['r', '24 day'])
#plt.show()

""" binVels(partcoords, pivot, partvelocities, rhos, depthOfView, viewAngle,
            numOfGridPoints) """
#partvelsPLOT = binVels(partcoords, Becoords, partvels,
#                      4e3, i, [10, 10], 0.76, 0.5)
#counts, bins, bars = plt.hist(partvelsPLOT, bins=200)

"""
 for just Be, units are  cm = 6439255109019.17  g = 3.583999100000000E+034  s = 334148.986076665
 frequency for H aplha is 4.57e14
 wavelength is 656.28nm
"""


#unitvel = udist / utime
#
#diskMaxZ = 0.3
#diskMaxY = 1.8
#depth = 1.7e6
#
#partvelsPLOT4 = binVels(partcoords, Becoords, partvels,
#                      depth, 0, [10, 10], diskMaxY, diskMaxZ)
#
#angstroms4 = partvelsPLOT4 * unitvel / 3e8
#angstroms4 += 1
#angstroms4 *= 6562.8
#
#plt.hist(angstroms4, bins=100, normed='True')
#plt.xlabel('Angstroms')
#plt.ylabel('Normalised Counts')

#partvels = 656.28 * 1e-9 * (1 + unitvel * partvels / 3e8)

#for i in [1e2, 5e2, 1e3, 5e3, 1e4, 5e4, 1e5, 5e5]:
#    Becoords, NScoords, partcoords, partvels, rho, num_parts = readASCII(filename)
#    depth = i
#    partvelsPLOT4 = binVels(partcoords, Becoords, partvels,
#                          depth, 0, [10, 10], diskMaxY, diskMaxZ)
#    partvelsPLOT3 = binVels(partcoords, Becoords, partvels,
#                          depth, 30, [10, 10], diskMaxY, diskMaxZ)
#    partvelsPLOT2 = binVels(partcoords, Becoords, partvels,
#                          depth, 60, [10, 10], diskMaxY, diskMaxZ)
#    partvelsPLOT1 = binVels(partcoords, Becoords, partvels,
#                          depth, 90, [10, 10], diskMaxY, diskMaxZ)
#
##    f, (ax2, ax3, ax4) = plt.subplots(1, 3)
#    #ax1.hist(partvelsPLOT1, bins=200)
#    """
#    convert to angstroms
#    """
#
#    angstroms1 = partvelsPLOT1 * unitvel / 3e8
#    angstroms1 += 1
#    angstroms1 *= 6562.8
#
#    angstroms2 = partvelsPLOT2 * unitvel / 3e8
#    angstroms2 += 1
#    angstroms2 *= 6562.8
#
#    angstroms3 = partvelsPLOT3 * unitvel / 3e8
#    angstroms3 += 1
#    angstroms3 *= 6562.8
#
#    angstroms4 = partvelsPLOT4 * unitvel / 3e8
#    angstroms4 += 1
#    angstroms4 *= 6562.8
#
#    counts1, bins1, bars1 = plt.hist(angstroms1, bins=100, normed='True')
#    counts2, bins2, bars2 = plt.hist(angstroms2, bins=100, normed='True')
#    counts3, bins3, bars3 = plt.hist(angstroms3, bins=100, normed='True')
#    counts4, bins4, bars4 = plt.hist(angstroms4, bins=100, normed='True')
#    plt.close()
#    plt.plot(bins1[0:100], counts1, label='90 degrees')
#    plt.plot(bins2[0:100], counts2, label='60 degrees')
#    plt.plot(bins3[0:100], counts3, label='30 degrees')
#    plt.plot(bins4[0:100], counts4, label='0 degrees')
#    plt.legend()
#    plt.ylim(0, max(counts2)+0.02)
#    plt.xlabel('Angstroms')
#    plt.ylabel('Normalised Counts')
##    ax1.set_title('90 degrees')
##    ax2.set_title('60 degrees')
##    ax3.set_title('30 degrees')
##    ax4.set_title('0 degrees')
#    #ax1.set_xlim(-0.05, 0.05)
#    #ax2.set_ylim(0, 600)
#    #ax3.set_ylim(0, 600)
#    #ax4.set_ylim(0, 600)
#    #ax2.set_xlim(-6, 6)
#    #ax3.set_xlim(-6, 6)
#    #ax4.set_xlim(-6, 6)
#    #plt.show()
#    plt.savefig("/data/rob1g10/SPH/DATA/Spectra/justBe/Halpha_{}.png".format(i))
#    plt.close()
#    print "Done", i

#for i in range(0, 61, 30):
#    for depth in [5e2, 1e3, 5e3, 1e4, 5e4, 1e5]:
#        Becoords, NScoords, partcoords, partvels, rho, num_parts = readASCII(filename)
#        partvelsPLOT = binVels(partcoords, Becoords, partvels,
#                              depth, i, [10, 10], diskMaxY, diskMaxZ)
#
#        """
#        convert to angstroms
#        """
#
#        angstroms = partvelsPLOT * unitvel / 3e8
#        angstroms += 1
#        angstroms *= 6562.8
#
#        counts, bins = np.histogram(angstroms, bins=100, normed='True')
#        plt.plot(bins[0:100], counts, label='{} depth'.format(depth))
#    plt.xlabel('Angstroms')
#    plt.ylabel('Normalised Counts')
#    plt.legend(loc=2,prop={'size':6})
#    plt.title('{} degrees'.format(i))
#    plt.savefig("/data/rob1g10/SPH/DATA/Spectra/justBe/Halpha_{}degrees.png".format(i), dpi=100)
#    plt.close()
#    print "Done", i

""" just line plots """

#f, (ax2, ax3, ax4) = plt.subplots(1, 3)
#ax2.plot(np.arange(0,200), counts2)
#ax3.plot(np.arange(0,200), counts3)
#ax4.plot(np.arange(0,200), counts4)
#ax2.set_title('60 degrees')
#ax3.set_title('30 degrees')
#ax4.set_title('0 degrees')
#plt.show()

#plt.xlim(-0.1, 0.1)
#plt.title("0 degrees")
#plt.show()
#print countParts(Becoords, partcoords[2:], -2, [0.01, 0.01])

#plt.hist(np.sqrt(partcoords[113864:,2]**2 + partcoords[113864:,1]**2 + partcoords[113864:,0]**2)/0.08, bins=50)

"""----------------------------------------------------------------------------

                            ALL ROUND PLOTTING

----------------------------------------------------------------------------"""

#plotRotatingFigs([filename, diskname, IPname], "../Plots to Show/All Round Image/20d_0.4e_Spiral/", "20d_0.4e")

"""----------------------------------------------------------------------------

                               ANIMATIONS

----------------------------------------------------------------------------"""

#createAnimation('/data/rob1g10/SPH/DATA/Animations/inclined/inclined',
#                '/data/rob1g10/SPH/DATA/Animations/inclined/',
#                '/home/rob1g10/Dropbox/PhD/processedDISKs/')

#createAnimation('/data/rob1g10/SPH/DATA/Animations/injection_off/40d_0.4e/40d_0.4e',
#                '/data/rob1g10/SPH/DATA/Animations/injection_off/40d_0.4e/',
#                '/data/rob1g10/SPH/DATA/Animations/injection_off/40d_0.4e/processedDISKs/')
                
#createAnimation('/data/rob1g10/SPH/DATA/Animations/injBurst_40d_0.4e/injBurst_40d_0.4e',
#                '/data/rob1g10/SPH/DATA/Animations/injBurst_40d_0.4e/',
#                '/data/rob1g10/SPH/DATA/DISKS/injBurst_40d_0.4e/processedDISKs/',
#                '2d')

#createAnimation('/data/rob1g10/SPH/DATA/Animations/MSSL/MSSL280916',
#                '/data/rob1g10/SPH/DATA/Animations/MSSL/',
#                '/data/rob1g10/SPH/DATA/Animations/MSSL/processedDISKs/',
#                'NSaccretion')

#createAnimation('/data/rob1g10/SPH/DATA/Animations/M10_34.1d_0.4e/M10_34.1d_0.4e',
#                '/data/rob1g10/SPH/DATA/Animations/M10_34.1d_0.4e/',
#                '/data/rob1g10/SPH/DATA/Animations/M10_34.1d_0.4e/processedDISKs/',
#                '2d')
#
#createAnimation('/data/rob1g10/SPH/DATA/Animations/M1.4_43d_0.4e/M1.4_43d_0.4e',
#                '/data/rob1g10/SPH/DATA/Animations/M1.4_43d_0.4e/',
#                '/data/rob1g10/SPH/DATA/Animations/M1.4_43d_0.4e/processedDISKs/',
#                '2d')

#createAnimation('/data/rob1g10/SPH/DATA/Animations/20d_0.2e_eqTest/20d_0.2e_eqTest',
#                '/data/rob1g10/SPH/DATA/Animations/20d_0.2e_eqTest/',
#                '/data/rob1g10/SPH/DATA/Animations/20d_0.2e_eqTest/processedDISKs/',
#                '2d')

#createAnimation('/data/rob1g10/SPH/DATA/Animations/20d_0.9e/20d_0.9e',
#                '/data/rob1g10/SPH/DATA/Animations/20d_0.9e/',
#                '/data/rob1g10/SPH/DATA/Animations/20d_0.9e/processedDISKs/',
#                '2d')
"""----------------------------------------------------------------------------

                               FIND PARTS

----------------------------------------------------------------------------"""
#BOX = [-0.1, 0.1, -25, 25, -1.5, 1.5]
#partsInBox = findParts(partcoords, BOX)
#partsInBox2 = findParts(partcoords2, BOX)
#vels = [partvels[i,0] * udist * 0.01 / utime for i in partsInBox]
#vels2 = [partvels2[i,0] * udist * 0.01 / utime for i in partsInBox2]
##print vels
#plt.hist(vels2, bins = 100)

#partcoords = partcoords[:10]

"""----------------------------------------------------------------------------

                               20d EW

----------------------------------------------------------------------------"""

#angle = range(0, 91, 10)
#
#Ave20d00e = [15.19, 14.34, 13.33, 12.13, 10.85, 10.01, 9.29, 8.42, 7.25, 6.17]
#Ave20d02e = [10.42, 9.54, 8.54, 7.47, 6.56, 6.11, 5.71, 5.24, 4.59, 4.01]
#Ave20d04e = [5.48, 4.58, 3.81, 3.14, 2.78, 2.63, 2.51, 2.36, 2.16, 1.97]
#
#SD20d00e = [0.10, 0.22, 0.18, 0.22, 0.22, 0.13, 0.09, 0.09, 0.11, 0.11]
#SD20d02e = [0.72, 0.24, 0.09, 0.19, 0.23, 0.11, 0.10, 0.16, 0.22, 0.23]
#SD20d04e = [0.39, 0.12, 0.10, 0.11, 0.09, 0.08, 0.09, 0.11, 0.12, 0.12]
#
#plt.plot(angle, Ave20d00e, label="0.0e")
#plt.plot(angle, Ave20d02e, label="0.2e")
#plt.plot(angle, Ave20d04e, label="0.4e")
#plt.legend()
#plt.xlabel("Viewing Angle")
#plt.ylabel("-EW Average")


