import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import pylab as plt
import os
from RorysSPHpackage import readQDISK, readIPDISK, readUnits, readFILES, calcPeriastronDistance, findDiskSize
from RorysSPHpackageOpticalDepth import readOptDepthData4py 


"""----------------------------------------------------------------------------

                Escape Velocity

----------------------------------------------------------------------------"""

def escapeVelPlot(partcoords, partvels, Bemass, Rmax, udist, utime):
    """
    Sorts particles into rings of radius and finds how many are at escape 
    velocity at that radius.

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
                                         
    partvels ((npart, 3) float array) - x,y,z velocities of the disk
                                         particles

    Bemass (float) - mass of the Be star in solar masses.
    
    Rmax (float) - the maximum radius to check up to in stellar radii.
    
    udist, utime (floats) - units for simulation.
    

    OUTPUT VARIABLES

    Plot of number of particles above escape velocity at each radius interval.

    """
    
    magvels = 1e-5 * np.sqrt(partvels[:,0]**2 + partvels[:,1]**2 + partvels[:,2]**2) * 0.01 * udist / utime
    magdists = np.sqrt(partcoords[:,0]**2 + partcoords[:,1]**2 + partcoords[:,2]**2) * 0.01 * udist / (7 * 695700e3)
    
    escapeVel = 1e-5 * np.sqrt(2 * 6.67e-11 * Bemass * 2e30 / (np.arange(1, Rmax) * 7 * 695700e3))
    kepVel = escapeVel / np.sqrt(2)
        
    plt.scatter(magdists, magvels, s=4)
    plt.plot(range(1,Rmax), escapeVel)
    plt.plot(range(1,Rmax), kepVel)
    plt.xlim(1, Rmax)
    plt.xlabel("Distance from centre of Be star in stellar radii")    
    plt.ylabel("Velocity in 100 km/s")
    

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
    
   
def plotSpectralLine(partcoords, pivot, partvelocities, trans, viewAngle,
                     numOfGridPoints, diskMaxY, diskMaxZ, BINS, udist, utime, emissivity, 
                     name):
    """
    Bins velocities of particles into a histogram such that a spectra can be
    made and then creates a line plot out of the histogram.

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

    trans (array of ints) - transmittance of each particle.

    viewAngle (degrees) - the angle at which the observer is viewing the disk
                        (90 degrees is considered to be directly above the star)

    numOfGridPointsX (2 integer array) - defines the number of squares which
                                         the disk is split into
                                         (y and z directions).

    diskMaxY, diskMaxZ - length from origin that the disk extends to in y and
                         z coordinates respectively in computational units
                         (plot the system then estimate)
                         
    BINS (int) - number of bins for the histogram (corresponds to resolution 
                 of line plot) 
                         
    udist, utime (floats) - units of the simulation (distance and time 
                            respectively)
                            
    emissivity - the type of emissivity given by the disk, currenly only
                 "uniform" and "R-3/2"

    OUTPUT VARIABLES
    ----------------

    A array of velocities for plotting.

    """

    " Replace input array so that it is not altered "

    partvels_func = partvelocities[:,0]
    partvels_funcPlot = np.array([])
    R = np.sqrt(partcoords[:,0]**2 + partcoords[:,1]**2)
    
    partcoords_func = partcoords

    " emissivity "
    if(emissivity == "R-3/2"):
        emissScaling = R ** (-1.5)
    if(emissivity == "uniform"):
        emissScaling = np.ones_like(R)
    emissScaling /= min(emissScaling)
    emissScalingPlot = np.array([])
    
    " optical depth using Thomson scattering "
#    opticalDepthScaling = calcOpticalDepth(partcoords * (udist*0.01), rho * umass/((0.01*udist)**3), 
#                            0.1 * max(partcoords[:,0])  * (udist*0.01), 9.1e-31, 5)
#    opticalDepthScalingPlot = opticalDepthScaling / np.min(opticalDepthScaling)

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
                                           sizeOfSquare, 1.7e6))

    for i in range(len(elements)):
        emissScalingPlot = np.append(emissScalingPlot,
                                     emissScaling[elements[i]])        
        partvels_funcPlot = np.append(partvels_funcPlot,
                                          partvels_func[elements[i]])
    
    unitvel = 0.01 * udist / utime
    
    angstroms = partvels_funcPlot * unitvel / 3e8   
    angstroms += 1
    angstroms *= 6562.8
    
    counts, bins = np.histogram(angstroms, bins=BINS)
    binDiff = bins[1] - bins[0]
    
    " Sort angstroms into the bins manually and scale based on emissivity "
    countsPlot = np.zeros_like(counts)
    for i in range(len(angstroms)):
        binNum = (angstroms[i] - bins[0]) / binDiff
        if(int(binNum) < 99):
            countsPlot[int(binNum)] += emissScalingPlot[i]
#            countsPlot[int(binNum)] += trans[i]
            
    countsPlot = countsPlot / float(np.sum(countsPlot))
    
    plt.plot(bins[0:-1], countsPlot, label=name)

    plt.xlabel('Angstroms')
    plt.ylabel('Counts')
    
#    return opticalDepthScaling
    
    
def plotAndSaveSpectralLine(CcodeName, filenames, baseName, thetaRange, phiRange,
                            BINS, emissivity, outFilename, save='n', xlims="n", T=True):
    """
    Bins velocities of particles into a histogram such that a spectra can be
    made and then creates a line plot out of the histogram.

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

    trans (array of ints) - transmittance of each particle.

    viewAngle (degrees) - the angle at which the observer is viewing the disk
                        (90 degrees is considered to be directly above the star)
                         
    BINS (int) - number of bins for the histogram (corresponds to resolution 
                 of line plot) 
                         
    udist, utime (floats) - units of the simulation (distance and time 
                            respectively)
                            
    emissivity - the type of emissivity given by the disk, currenly only
                 "uniform" and "R-3/2"

    OUTPUT VARIABLES
    ----------------

    A array of velocities for plotting.

    """
    for theta1 in range(thetaRange[0], thetaRange[1], thetaRange[2]):
        for phi1 in range(phiRange[0], phiRange[1], phiRange[2]):
            for filename in filenames:
                diskname = filename[:-26] + "disk001"
                IPname = filename[:-11] + "IPDISK" + filename[-5:-2]
                Becoords, NScoords, partcoords, partvelocities, rho, u, num_parts, time, umass, udist, utime = readFILES(filename, diskname, IPname)
                
                R = np.sqrt(partcoords[:,0]**2 + partcoords[:,1]**2) * udist / (7 * 6.96e10)                
                
                " Transmission "
                if(T):
                    OD = readOptDepthData4py("/data/rob1g10/SPH/DATA/DISKS/{}/optDepthData{}{}_theta{}_phi{}_time{:0.1f}.txt".format(baseName, CcodeName, baseName, theta1, phi1, time * utime / (24 * 3600)))
                    trans = (np.exp(-OD)) #* 1e2
                    print np.sum(trans)
                else:
                    trans = np.ones_like(OD)
                    
                " emissivity "
                if(emissivity == "R-3/2"):
                    emissScaling = R ** (-1.5)
                if(emissivity == "uniform"):
                    emissScaling = np.ones_like(partvelocities)
                # emissScaling = emissScaling / min(emissScaling)
            
                " Convert angle to radians "
            
                theta = np.deg2rad(theta1)
                phi = np.deg2rad(phi1)
            
                """ Arrange velocities about y-axis so we only require z component"""
            
                partvelx = partvelocities[:,0] * np.cos(phi) - partvelocities[:, 1] * np.sin(phi)
                partvely = (partvelocities[:, 0] * np.sin(phi) + partvelocities[:, 1] * np.cos(phi)) * np.cos(theta) - partvelocities[:,2] * np.sin(theta)
                partvelz = partvelocities[:, 1] * np.sin(theta) + partvelocities[:, 2] * np.cos(theta)
                
                unitvel = 0.01 * udist / utime
                
                angstroms = partvelz * unitvel / 3e8   
                angstroms += 1
                angstroms *= 6562.8
                
                "calculate keplerian velocity at stellar radius"
                kepVel = np.sqrt(6.67e-11 * 18 * 2e30 / (7 * 695700e3)) * np.sin(theta)
                #print np.sin(theta)
                
                kepAngstroms = kepVel / 3e8   
                kepAngstroms += 1
                kepAngstroms *= 6562.8
                
                counts, bins = np.histogram(angstroms, bins=BINS)
                binDiff = bins[1] - bins[0]
                
                " Sort angstroms into the bins manually and scale based on emissivity and transmissivity "
                countsPlot = np.zeros_like(counts)
                for i in range(len(angstroms)):
                    binNum = (angstroms[i] - bins[0]) / binDiff
                    #print(int(binNum))
                    if(R[i]<10):
                        if(int(binNum) < BINS-1):
                            countsPlot[int(binNum)] += trans[i] #* emissScaling[i]
                        else:
                            countsPlot[BINS-1] += trans[i] #* emissScaling[i]
                        
                countsPlot = countsPlot# / float(np.sum(countsPlot))
                
                plt.plot([kepAngstroms for i in range(20)], np.linspace(0, max(countsPlot), 20), color='black')
                plt.plot(bins[0:-1], countsPlot, label="{:01}degs".format(theta1))
            
                plt.xlabel('Angstroms')
                plt.ylabel('Counts')
                plt.legend()
                if(xlims=='y'):
                    plt.xlim(6540, 6580)
                
    if(save == 'y'):
        plt.savefig(outFilename, figsize=(18, 16), dpi=100)
                    
    plt.show()
            
#                plt.close()
    

def showPeakRelation4Spectra(partcoords, pivot, partvelocities,
                             depthOfView, numOfGridPoints, diskMaxY, diskMaxZ):
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

    numOfGridPoints (2 integer array) - defines the number of squares which
                                         the disk is split into
                                         (y and z directions).

    OUTPUT VARIABLES
    ----------------

    A plot of the movement of the peaks with viewing angle.

    """

    " Create array for the plot at the end "

    plotArray = np.array([])
    partvels_funcPLOT = np.array([])

    " Loop over angles and find the peaks "
    for i in range(0, 91, 5):
        partvels_funcPLOT = binVels(partcoords, pivot, partvelocities,
                                        depthOfView, i,
                                        numOfGridPoints, diskMaxY, diskMaxZ)
        counts, bins = np.histogram(partvels_funcPLOT,
                                      bins=100)
        #print counts
        plotArray = np.append(plotArray, 2*bins[np.argmax(counts)])
#        print plotArray
        del counts
        del bins
        del partvels_funcPLOT
        
    plotArray /= max(plotArray)
    print plotArray
    plt.plot(range(0,91,5), plotArray, label="peak spacing")
    plt.plot(range(0,91), np.cos(np.radians(range(0,91))), label="cosine of viewing angle")
    #plt.show(block=False)
    plt.legend()
    plt.title("Peak spacing relation")
    plt.xlabel("Viewing angle")
    plt.ylabel("Peak spacing")
#    plt.ylim(0,1.1)
    plt.show()
    return plotArray
   
   
def plotHaContributions(partcoords, pivot, partvelocities, viewAngle, BINS,
                        radiusSteps, maxR, udist, utime):
    """
    Bins velocities of particles into a histogram such that a spectra can be
    made and then creates a line plot out of the histogram.

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

    viewAngle (degrees) - the angle at which the observer is viewing the disk
                        (90 degrees is considered to be directly above the star)
                         
    BINS (int) - number of bins for the histogram (corresponds to resolution 
                 of line plot)
                 
    radiusSteps (int) - the number of sections the disk will be split into in
                        radius space.
    
    maxR (float) - radius at which disk should not be included in process.
                         
    udist, utime (floats) - units of the simulation (distance and time 
                            respectively)

    OUTPUT VARIABLES
    ----------------

    A array of velocities for plotting.

    """
    partcoords -= pivot
    viewAngle = np.deg2rad(viewAngle)
    R = np.sqrt(partcoords[:,0]**2 + partcoords[:,1]**2)
    unitvel = 0.01 * udist / utime
    
    maxR /= (udist / (7 * 695700e3 * 100))
   
    if(maxR > max(R)):
        maxR = max(R)
        
    print(maxR)
    
    for i in range(0, radiusSteps):
        
        outerRad = (i + 1) * (maxR / radiusSteps)
        innerRad = i * (maxR / radiusSteps)
        
        " Replace input array so that it is not altered "
    
        partvels_func = np.array([])
    
        partcoords_func = np.array([])
        
        for i in range(len(partcoords)):
            if(R[i] > innerRad and R[i] < outerRad):
                partvelx = partvelocities[i][0]
                partvelz = partvelocities[i][2]
                partvels_func = np.append(partvels_func, partvelx * np.cos(viewAngle) + partvelz * np.sin(viewAngle))
        
        angstroms = partvels_func * unitvel / 3e8   
        angstroms += 1
        angstroms *= 6562.8
        print int(BINS * (float(len(angstroms)) / len(partvels)))
        if(int(BINS * (float(len(angstroms)) / len(partvels))) == 0):
            counts, bins = np.histogram(angstroms, bins=BINS)
        else:
            counts, bins = np.histogram(angstroms, bins=BINS)#int(BINS * (float(len(angstroms)) / len(partvels))))
    
#        plt.plot(bins[0:-1], counts, label="{:.3} to {:.3} ({} counts = {:.3}%)".format(innerRad * udist / (7 * 695700e3 * 100), 
#                                                                                    outerRad * udist / (7 * 695700e3 * 100), 
#                                                                                    len(partvels_func),
#                                                                                    100*len(partvels_func)/float(len(partvels))))
        plt.plot(bins[0:-1], counts, label="{:.3} to {:.3} ({:.3}%)".format(innerRad * udist / (7 * 695700e3 * 100), 
                                                                                    outerRad * udist / (7 * 695700e3 * 100),
                                                                                    100*len(partvels_func)/float(len(partvels))))
                                                                                    
    partvels_func = partvelocities[:, 0] * np.cos(viewAngle) + partvelocities[:, 2] * np.sin(viewAngle)
    angstroms = partvels_func * unitvel / 3e8   
    angstroms += 1
    angstroms *= 6562.8
    
    counts, bins = np.histogram(angstroms, bins=BINS)#, normed='True')
    
    plt.plot(bins[0:-1], counts, label="Complete Halpha")     
    plt.legend(loc="upper left", prop={'size':12})
    plt.xlabel('Angstroms')
    plt.ylabel('Counts')
    
    return


def plotFeatureVels(CcodeName, filenames, baseName, thetaRange, phiRange, binSteps,
                            BINS, emissivity, outFilename, show='n', save='n', xlims="n", T=True):
    """
    Bins velocities of particles into a histogram such that a spectra can be
    made and then creates a line plot out of the histogram.

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

    viewAngle (degrees) - the angle at which the observer is viewing the disk
                        (90 degrees is considered to be directly above the star)
                         
    BINS (int) - number of bins for the histogram (corresponds to resolution 
                 of line plot)
                 
    radiusSteps (int) - the number of sections the disk will be split into in
                        radius space.
    
    maxR (float) - radius at which disk should not be included in process.
                         
    udist, utime (floats) - units of the simulation (distance and time 
                            respectively)

    OUTPUT VARIABLES
    ----------------

    A array of velocities for plotting.

    """
    binSteps = 5
    
    for theta1 in range(thetaRange[0], thetaRange[1], thetaRange[2]):
        for phi1 in range(phiRange[0], phiRange[1], phiRange[2]):
            for filename in filenames:
                diskname = filename[:-26] + "disk001"
                IPname = filename[:-11] + "IPDISK" + filename[-5:-2]
#                diskname = "DATA/disk001"
#                IPname = "DATA/IPDISK" + filename[-5:-2]
                Becoords, NScoords, partcoords, partvelocities, rho, u, num_parts, time, umass, udist, utime = readFILES(filename, diskname, IPname)
                
                R = np.sqrt(partcoords[:,0]**2 + partcoords[:,1]**2) * udist / (7 * 6.96e10) 
                
                Becoords *= udist / (7 * 6.96e10)
                NScoords *= udist / (7 * 6.96e10)
                partcoords *= udist / (7 * 6.96e10)
                
                " Transmission "
                if(T):
                    OD = readOptDepthData4py("/data/rob1g10/SPH/DATA/DISKS/{}/optDepthData{}{}_theta{}_phi{}_time{:0.1f}.txt".format(baseName, CcodeName, baseName, theta1, phi1, time * utime / (24 * 3600)))
#                    OD = readOptDepthData4py(CcodeName)
                    trans = (np.exp(-OD)) * 1e2
                    print np.sum(trans)
                else:
                    trans = np.ones_like(OD)
                    
                " emissivity "
                if(emissivity == "R-3/2"):
                    emissScaling = R ** (-1.5)
                if(emissivity == "uniform"):
                    emissScaling = np.ones_like(partvelocities)
                # emissScaling = emissScaling / min(emissScaling)
            
                " Convert angle to radians "
            
                theta = np.deg2rad(theta1)
                phi = np.deg2rad(phi1)
            
                """ Arrange velocities about y-axis so we only require z component"""
            
                partvelx = partvelocities[:,0] * np.cos(phi) - partvelocities[:, 1] * np.sin(phi)
                partvely = (partvelocities[:, 0] * np.sin(phi) + partvelocities[:, 1] * np.cos(phi)) * np.cos(theta) - partvelocities[:,2] * np.sin(theta)
                partvelz = partvelocities[:, 1] * np.sin(theta) + partvelocities[:, 2] * np.cos(theta)
                
                unitvel = 0.01 * udist / utime
                
                angstroms = partvelz * unitvel / 3e8   
                angstroms += 1
                angstroms *= 6562.8
                
                "calculate keplerian velocity at stellar radius"
                kepVel = np.sqrt(6.67e-11 * 18 * 2e30 / (7 * 695700e3)) * np.sin(theta)
                #print np.sin(theta)
                
                kepAngstroms = kepVel / 3e8   
                kepAngstroms += 1
                kepAngstroms *= 6562.8
                
                counts, bins = np.histogram(angstroms, bins=BINS)
                binDiff = bins[1] - bins[0]
#                
#                " Sort angstroms into the bins manually and scale based on emissivity and transmissivity "
#                countsPlot = np.zeros_like(counts)
#                for i in range(len(angstroms)):
#                    binNum = (angstroms[i] - bins[0]) / binDiff
#                    #print(int(binNum))
#                    if(R[i]<10):
#                        if(int(binNum) < BINS-1):
#                            countsPlot[int(binNum)] += trans[i] #* emissScaling[i]
#                        else:
#                            countsPlot[BINS-1] += trans[i] #* emissScaling[i]
                        
#                countsPlot = countsPlot# / float(np.sum(countsPlot))
                        
                """ SORT VELS """
                
                sortOrder = np.argsort(partvelz)
                
                partvelz = partvelz[sortOrder]
                angstroms = angstroms[sortOrder]
                xcoords = partcoords[:,0][sortOrder]
                ycoords = partcoords[:,1][sortOrder]
                
                """ FIND SEPARATIONS IN ANGSTROMS """

#                seps = np.array([5.8483e4, 7.59810000e+04, 9.75580000e+04, 
#                                 int(0.5*(9.75580000e+04 + 1.15514000e+05)), 1.15514000e+05, -1.00000000e+00])
#                seps = np.array([])                
#                for i in range(binSteps):
#                    thing = np.abs(partvelz) - (i * np.max(partvelz) / float(binSteps))
#                    seps = np.append(seps, np.argmin(np.abs(thing[7e3:])))
##                    print np.abs(partvelz) - (i * np.max(partvelz) / float(binSteps))
#                seps = np.append(seps, -1)
#        
#                print seps
#                print np.argmin(np.abs(partvelz))
                
#                seps = np.array([])
#                thing = np.abs(partvelz)
#                seps = np.append(seps, np.argmin(np.abs(thing)))
#                for i in range(1, binSteps):
#                    thing = np.abs(partvelz) - (i * np.max(partvelz) / float(binSteps))
#                    seps = np.append(seps, np.argmin(np.abs(thing[seps[0]:])))
##                    print np.abs(partvelz) - (i * np.max(partvelz) / float(binSteps))
#                seps = np.append(seps, -1)
#        
#                print seps
#                print np.argmin(np.abs(partvelz))

                thing = np.abs(partvelz)
                cutoff = np.argmin(np.abs(thing))
                partvelzPositive = partvelz
                partvelzPositive = np.delete(partvelzPositive, range(0, int(cutoff)))
                print "HERE:", partvelzPositive
                
                seps = np.array([])
                seps = np.append(seps, 0)
                for i in range(1, binSteps):
                    thing = np.abs(partvelzPositive) - (i * np.max(partvelzPositive) / float(binSteps))
                    seps = np.append(seps, np.argmin(np.abs(thing)))
#                    print np.abs(partvelz) - (i * np.max(partvelz) / float(binSteps))
                seps = np.append(seps, -1 - cutoff)
        
                print seps
                print np.argmin(np.abs(partvelzPositive))
                
                """ PLOTTING """                
                f, axarr = plt.subplots(1, 2, figsize=(16, 8))
                
#                axarr[0].scatter(xcoords, ycoords, s=0.5, color='r')
                
#                axarr[1].plot([kepAngstroms for i in range(20)], np.linspace(0, max(countsPlot), 20), color='black')
                axarr[1].hist(angstroms[0:cutoff], bins=BINS/2)
            
                axarr[1].set_xlabel('Angstroms')
                #axarr[1].set_ylabel('Counts')
                axarr[1].axes.get_yaxis().set_visible(False)

                colours = ['red', 'yellow', 'green', 'orange', 'black'] 
#                binsMax = max(bins)
                
                for i in range(binSteps):
                    axarr[1].hist(angstroms[cutoff+seps[i]:cutoff+seps[i+1]], color=colours[i], 
                                                bins=10)
#                    axarr[1].hist(angstroms[-seps[i]:-seps[i+1]], color=colours[i])
                    axarr[0].scatter(xcoords[cutoff+seps[i]:cutoff+seps[i+1]], ycoords[cutoff+seps[i]:cutoff+seps[i+1]], s=0.5, color=colours[i])
                    axarr[0].scatter(xcoords[0:cutoff], ycoords[0:cutoff], s=0.5, color="blue")
                    print angstroms[seps[i]:seps[i+1]]
                
                axarr[0].scatter(NScoords[0], NScoords[1], color="black")
                axarr[0].set_xlim(-11, 11)
                axarr[0].set_ylim(-11, 11)
                axarr[1].set_ylim(0, 300)
#                axarr[1].legend()
#                if(xlims=='y'):
#                    plt.xlim(6540, 6580)
                
    if(save == 'y'):
        plt.savefig(outFilename, figsize=(18, 16), dpi=100)
    
    if(show=='y'):               
        plt.show()
    if(show=='n'):
        plt.close()
    
    return


def calcEW(partvelocities, T, partFlux, fcMax, fcMin, udist, utime, viewAngles=[0,0], intSteps=200):
    """

    INPUT VARIABLES
    ---------------

    partvelocities (array of floats) - array containing the velocities
                                       of each particle. 

    T (array of floats) - transmission coefficient for all particles.
    
    partFlux (array of floats) - the flux of each particle. 
    
    fcMax, fcMin (floats) - values for the flux continuum on either side of the
                            emission feature.

    intSteps (int) - number of integration steps.

    OUTPUT VARIABLES
    ----------------

    The equivalent width of the given array.

    """
    
    " Convert angle to radians "

    viewAngles[0] = np.deg2rad(viewAngles[0])
    viewAngles[1] = np.deg2rad(viewAngles[1])

    """ Arrange velocities about y-axis so we only require z component"""
    
    phi = viewAngles[0]
    theta = viewAngles[1]

    partvelx = partvelocities[:,0] * np.cos(phi) - partvelocities[:, 1] * np.sin(phi)
    partvely = (partvelocities[:, 0] * np.sin(phi) + partvelocities[:, 1] * np.cos(phi)) * np.cos(theta) - partvelocities[:,2] * np.sin(theta)
    partvelz = partvelocities[:, 1] * np.sin(theta) + partvelocities[:, 2] * np.cos(theta)
    
    unitvel = 0.01 * udist / utime
    
    angstroms = partvelz * unitvel / 3e8   
    angstroms += 1
    angstroms *= 6562.8
    
    """ USING INTEGRATION --------------------------------------------------"""
    
    intSum = 0
    a = min(angstroms)
    b = max(angstroms)
    stepSize = (b-a) / intSteps
    intRegion = np.linspace(a, b, intSteps)
    "approximate continuum emission as a straight line"
    F_c = fcMin + (intRegion * (fcMax - fcMin) / (b-a))
    F_c /= intSteps
    "feature emission"
    F_lambda = np.array(F_c)
    for i in range(len(angstroms)):
        binNum = (angstroms[i] - a) / stepSize
        if(int(binNum) < intSteps-1):
            F_lambda[int(binNum)] += T[i] * partFlux[i]
        else:
            F_lambda[intSteps-1] += T[i] * partFlux[i]
    for i in range(len(intRegion)):
        if(i == 0 or i == len(intRegion)-1):
            intSum += 0.5 * stepSize * (1 - (F_lambda[i] / F_c[i]))
        else:
            intSum += stepSize * (1 - (F_lambda[i] / F_c[i]))
    
    ans = intSum 
            
    """ JUST ADDING FLUX ---------------------------------------------------"""
#    
#    RANGE = (max(angstroms) - min(angstroms))
#    ans = - np.sum(T*partFlux) * 8.0175e-7
    
    return ans
    

def calcEWAverageWithTime(inputdirectory, startDISKs, numDISKs, phi, theta):
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
    IPfilenames = np.array([])
    
    f = open("lsQG.txt", "w")
    for j in range(startDISKs, startDISKs + numDISKs):
        for k in range(1,100):
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
    
    """ IPDISKs """
    
    f = open("lsIP.txt", "r")
    
    lines = f.readlines()
    
    for line in lines:
        IPfilenames = np.append(IPfilenames, line[0:-1])

    f.close()

    times = np.array([])
    EW = 0
    
#    for i in range(len(IPfilenames)):
#        times = np.append(times, readIPDISK(IPfilenames[i]))
        
    for i in range(len(QGfilenames)):
        filename = QGfilenames[i]
        diskname = filename[:-26] + "disk001"
        IPname = filename[:-11] + "IPDISK" + filename[-5:-2]
        Becoords, NScoords, partcoords, partvels, rho, u, num_parts, time, umass, udist, utime = readFILES(filename, diskname, IPname)
        Becoords *= udist / (7 * 6.96e10)
        NScoords *= udist / (7 * 6.96e10)
        partcoords *= udist / (7 * 6.96e10)
#        baseName = QGfilenames[i][29:37]    
#        oD = readOptDepthData4py("/data/rob1g10/SPH/DATA/DISKS/20d_0.4e/optDepthData20d_0.4e_theta{}_phi{}_time{:0.1f}.txt".format(theta, phi, time * utime / (24 * 3600)))
#        T = np.exp(-oD) / 1e-3
        EW += calcEW(partvels, np.ones(len(partvels)), np.ones(len(partvels)), 1.021e6, 1.021e6, udist, utime, intSteps=200, viewAngles=[phi,theta])

    os.system("rm lsQG.txt lsIP.txt")
    
    return EW / 20
    
    
def calcPeakSeparation(partvelocities, trans, viewAngles,
                            BINS, udist, utime, emissivity):
    """
    Bins velocities of particles into a histogram such that a spectra can be
    made and then creates a line plot out of the histogram.

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

    trans (array of ints) - transmittance of each particle.

    viewAngle (degrees) - the angle at which the observer is viewing the disk
                        (90 degrees is considered to be directly above the star)
                         
    BINS (int) - number of bins for the histogram (corresponds to resolution 
                 of line plot) 
                         
    udist, utime (floats) - units of the simulation (distance and time 
                            respectively)
                            
    emissivity - the type of emissivity given by the disk, currenly only
                 "uniform" and "R-3/2"

    OUTPUT VARIABLES
    ----------------

    A array of velocities for plotting.

    """

    " emissivity "
    if(emissivity == "R-3/2"):
        emissScaling = R ** (-1.5)
    if(emissivity == "uniform"):
        emissScaling = np.ones_like(partvelocities)
    # emissScaling = emissScaling / min(emissScaling)

    " Convert angle to radians "

    viewAngles[0] = np.deg2rad(viewAngles[0])
    viewAngles[1] = np.deg2rad(viewAngles[1])

    """ Arrange velocities about y-axis so we only require z component"""
    
    phi = viewAngles[0]
    theta = viewAngles[1]

    partvelx = partvelocities[:,0] * np.cos(phi) - partvelocities[:, 1] * np.sin(phi)
    partvely = (partvelocities[:, 0] * np.sin(phi) + partvelocities[:, 1] * np.cos(phi)) * np.cos(theta) - partvelocities[:,2] * np.sin(theta)
    partvelz = partvelocities[:, 1] * np.sin(theta) + partvelocities[:, 2] * np.cos(theta)
    
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
            countsPlot[int(binNum)] += trans[i] #* emissScaling[i]
        else:
            countsPlot[BINS-1] += trans[i] #* emissScaling[i]
            
    max1 = 0
    max2 = 0     
    peak1 = 0
    peak2 = 0
    for i in range(1, len(countsPlot)/2):
        if(countsPlot[i-1] < countsPlot[i] and countsPlot[i+1] < countsPlot[i]):
            if(countsPlot[i] > max2):
                max2 = countsPlot[i]
                peak2 = bins[i]
    for i in range(len(countsPlot)/2, len(countsPlot) - 1):
        if(countsPlot[i-1] < countsPlot[i] and countsPlot[i+1] < countsPlot[i]):
            if(countsPlot[i] > max1):
                max1 = countsPlot[i]
                peak1 = bins[i]
    
    if(peak1 == 0 or peak2 ==0):
        peak1 = peak2 = 0
    return [peak2, peak1]


def plotEW2dHist(baseName, filename, phi, theta, show='y', save='n', outFilename="figure", cbar=[0,60],
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
    
    EWsX = np.array([])
    EWsY = np.array([])
    
    diskname = filename[:-26] + "disk001"
    IPname = filename[:-11] + "IPDISK" + filename[-5:-2]
    Becoords, NScoords, partcoords, partvels, rho, u, num_parts, time, umass, udist, utime = readFILES(filename, diskname, IPname)
    Becoords *= udist / (7 * 6.96e10)
    NScoords *= udist / (7 * 6.96e10)
    partcoords *= udist / (7 * 6.96e10)
    oD = readOptDepthData4py("/data/rob1g10/SPH/DATA/DISKS/{}/optDepthData{}_theta{}_phi{}_time{:0.1f}.txt".format(baseName, 
                             baseName, theta, phi, time * utime / (24 * 3600)))
    T = np.exp(-oD) * 1e3
    for i in range(len(T)):
        for j in np.linspace(0, T[i], 1):
#            if(partcoords[i,0] < xylims[1] and partcoords[i,0] > xylims[0]
#               and partcoords[i,1] < xylims[3] and partcoords[i,1] > xylims[2]):
           EWsX = np.append(EWsX, partcoords[i,0])
           EWsY = np.append(EWsY, partcoords[i,1])
                   
    print(len(EWsX))
#    plt.hist2d(EWsX, EWsY, bins=200)#, normed=1)
    BINS = 200
    h, xedges, yedges = np.histogram2d(EWsX, EWsY, bins=BINS)
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

     
def plotEW2dHistDiff(baseName, filename, filename2, phi, theta, show='y', save='n', 
                     outFilename="figure", cbar=[0,60], xylims=[-25,25,-25,25]):
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
    
    EWsX = np.array([])
    EWsY = np.array([])
    
    diskname = filename[:-26] + "disk001"
    IPname = filename[:-11] + "IPDISK" + filename[-5:-2]
    Becoords, NScoords, partcoords, partvels, rho, u, num_parts, time, umass, udist, utime = readFILES(filename, diskname, IPname)
    Becoords *= udist / (7 * 6.96e10)
    NScoords *= udist / (7 * 6.96e10)
    partcoords *= udist / (7 * 6.96e10)
    oD = readOptDepthData4py("/data/rob1g10/SPH/DATA/DISKS/{}/optDepthData{}_theta{}_phi{}_time{:0.1f}.txt".format(baseName, 
                             baseName, theta, phi, time * utime / (24 * 3600)))
    T = np.exp(-oD) * 1e3
    for i in range(len(T)):
        for j in np.linspace(0, T[i], 1):
#            if(partcoords[i,0] < xylims[1] and partcoords[i,0] > xylims[0]
#               and partcoords[i,1] < xylims[3] and partcoords[i,1] > xylims[2]):
                   EWsX = np.append(EWsX, partcoords[i,0])
                   EWsY = np.append(EWsY, partcoords[i,1])
                   
    EWsX2 = np.array([])
    EWsY2 = np.array([])
    
    diskname2 = filename2[:-26] + "disk001"
    IPname2 = filename2[:-11] + "IPDISK" + filename2[-5:-2]
    Becoords2, NScoords2, partcoords2, partvels2, rho2, u2, num_parts2, time2, umass2, udist2, utime2 = readFILES(filename2, diskname2, IPname2)
    Becoords2 *= udist2 / (7 * 6.96e10)
    NScoords2 *= udist2 / (7 * 6.96e10)
    partcoords2 *= udist2 / (7 * 6.96e10)
    oD2 = readOptDepthData4py("/data/rob1g10/SPH/DATA/DISKS/{}/optDepthData{}_theta{}_phi{}_time{:0.1f}.txt".format(baseName, 
                             baseName, theta, phi, time2 * utime2 / (24 * 3600)))
    T2 = np.exp(-oD2) * 1e3
    for i in range(len(T2)):
        for j in np.linspace(0, T2[i], 1):
#            if(partcoords2[i,0] < xylims[1] and partcoords2[i,0] > xylims[0]
#               and partcoords2[i,1] < xylims[3] and partcoords2[i,1] > xylims[2]):
                   EWsX2 = np.append(EWsX2, partcoords2[i,0])
                   EWsY2 = np.append(EWsY2, partcoords2[i,1])
    print(len(EWsX))               
    BINS = 200
    h1, xedges, yedges = np.histogram2d(EWsX, EWsY, bins=BINS)
    h2, xedges2, yedges2 = np.histogram2d(EWsX2, EWsY2, bins=BINS)
#    print xedges - xedges2
#    print yedges - yedges2
    h = (h1 - h2) 
    
    plt.pcolor(xedges, yedges, h, cmap='RdBu')
    
#    EWsXplot = np.array([])
#    EWsYplot = np.array([])    
#    
#    for i in range(BINS):
#        print i
#        for j in range(BINS):
#            for k in range(int(h[i][j])):
#                EWsXplot = np.append(EWsXplot, xedges[i])
#                EWsYplot = np.append(EWsYplot, yedges[j])
#   
#    plt.hist2d(EWsXplot, EWsYplot, bins=BINS)
    plt.xlim(xylims[0], xylims[1])
    plt.ylim(xylims[2], xylims[3])
##    plt.clim(cbar[0], cbar[1])
##    plt.ylabel("EW")
##    plt.xlabel("Time / days")
    plt.colorbar()
    
    if(show == 'y'):
        plt.show()
    
    if(save == 'y'):
        plt.savefig(outFilename)#, figsize=(18, 16), dpi=100)

    if(show == 'n'):
        plt.close()
    
    
    
def plotEWwithTime(baseName, inDirs, startDISKs, numDISKs, phi, theta):
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
    
    for i in range(len(inDirs)):
        inputdirectory = inDirs[i]        
        
        os.system("rm lsQG.txt lsIP.txt")
        
        QGfilenames = np.array([])
        IPfilenames = np.array([])
        
        f = open("lsQG.txt", "w")
        for j in range(startDISKs[i], startDISKs[i] + numDISKs[i]):
            for k in range(1,100):
                NAME = "{}QGDISK{:03}{:02}\n".format(inputdirectory, j, k)
                if(os.system("ls {}".format(NAME)) == 0):
                    f.write(NAME)
                else:
                    break
        f.close()
        
        f = open("lsIP.txt", "w")
        for j in range(startDISKs[i], startDISKs[i] + numDISKs[i]):
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
        
        """ IPDISKs """
        
        f = open("lsIP.txt", "r")
        
        lines = f.readlines()
        
        for line in lines:
            IPfilenames = np.append(IPfilenames, line[0:-1])
    
        f.close()
    
        times = np.array([])
        EWs = np.array([])
        num_partsArr = np.array([])
        
        for j in range(len(IPfilenames)):
            times = np.append(times, readIPDISK(IPfilenames[j]))
            
        for j in range(len(QGfilenames)):
            filename = QGfilenames[j]
            diskname = filename[:-26] + "disk001"
            IPname = filename[:-11] + "IPDISK" + filename[-5:-2]
            Becoords, NScoords, partcoords, partvels, rho, u, num_parts, time, umass, udist, utime = readFILES(filename, diskname, IPname)
            Becoords *= udist / (7 * 6.96e10)
            NScoords *= udist / (7 * 6.96e10)
            partcoords *= udist / (7 * 6.96e10)  
            oD = readOptDepthData4py("/data/rob1g10/SPH/DATA/DISKS/{}/optDepthData{}_theta{}_phi{}_time{:0.1f}.txt".format(baseName[i], 
                                     baseName[i], theta, phi, time * utime / (24 * 3600)))
            T = np.exp(-oD)
#            EWs = np.append(EWs, calcEW(partvels, np.ones(len(partvels)), np.ones(len(partvels)), 1.565e5, 1.565e5, udist, utime, intSteps=200, viewAngles=[phi,theta]))
            EWs = np.append(EWs, np.sum(T))            
            num_partsArr = np.append(num_partsArr, num_parts)
             
        EWs /= 1000
        plt.plot(times * utime / (24 *3600), -EWs)
        plt.ylabel("EW")
        plt.xlabel("Time / days")
        MEAN = np.mean(-EWs)
        MAX =  np.max(-EWs)
        MIN =  np.min(-EWs)
        plt.text(1780, -28, "{:.2f}  +/- {:.2f}".format(MEAN, np.max([MAX - MEAN, MEAN - MIN])))
        plt.gca().invert_yaxis()
    plt.show()
#    plt.plot(times * utime / (24 *3600), num_partsArr)
#    plt.ylabel("Number of Particles")
#    plt.xlabel("Time / days")
#    plt.show()
    
    os.system("rm lsQG.txt lsIP.txt")


def calcQDISKsForOrbit(orbitalPeriod, baseName, QDISKstart, QDISKstepsMAX):
    
    tol = 0.025 * orbitalPeriod
    QDISKstepSTART = 1
    
    times = np.array([])
    maxRhos = np.array([])
    nparts = np.array([])
       
    end = [0,1]
    # calculate number of steps required for whole number of orbits
    for i in range(QDISKstart, QDISKstart+10):
        for j in range(QDISKstepSTART, QDISKstepsMAX):
            filename = '/data/rob1g10/SPH/DATA/DISKS/{}/processedDISKs/QGDISK{}{:02}'.format(baseName, i, j)
            diskname = filename[:-26] + "disk001"
            IPname = filename[:-11] + "IPDISK" + filename[-5:-2]
            Becoords, NScoords, partcoords, partvels, rho, u, num_parts, time, umass, udist, utime = readFILES(filename, diskname, IPname)          
            times  = np.append(times, time * utime / (24*3600))
            maxRhos = np.append(maxRhos, max(rho * umass / udist**3))
            nparts = np.append(nparts, num_parts)
#            print filename
            if(np.abs(times[0]-times[-1]) % orbitalPeriod < tol and np.abs(times[0]-times[-1]) / orbitalPeriod >= 1):
                end = [i, j]
                numOrbits = int(np.abs(times[0]-times[-1]) / orbitalPeriod)
                return end, np.mean(maxRhos), np.mean(nparts)
   
def foldEW(orbitalPeriod, baseName, CcodeName, QDISKstart, QDISKnum, QDISKstepsMAX, thetaRange, phi, fold='n'):
    
    tol = 0.025 * orbitalPeriod
    QDISKstepSTART = 1

    MEANs = np.array([])
    SDs = np.array([])
    
    times = np.array([])
        
    end = [0,1]
    # calculate number of steps required for whole number of orbits
    for i in range(QDISKnum):
        for j in range(QDISKstepSTART, QDISKstepsMAX):
            filename = '../../Desktop/40d_0.0e_eqTest/QGDISK{}{:02}'.format(QDISKstart+i, j)
            diskname = "../../Desktop/40d_0.0e_eqTest/disk001"
            IPname = "../../Desktop/40d_0.0e_eqTest/IPDISK" + filename[-5:-2]
            #filename = '/data/rob1g10/SPH/DATA/DISKS/{}/processedDISKs/QGDISK{}{:02}'.format(baseName, QDISKstart+i, j)
            #diskname = filename[:-26] + "disk001"
            #IPname = filename[:-11] + "IPDISK" + filename[-5:-2]
            Becoords, NScoords, partcoords, partvels, rho, u, num_parts, time, umass, udist, utime = readFILES(filename, diskname, IPname)
            try:
                oD = 0.48 * readOptDepthData4py("../../Desktop/40d_0.0e_eqTest/optDepthData{}{}_theta{}_phi{}_time{:0.1f}.txt".format(CcodeName, baseName, 0, phi, time * utime / (24 * 3600)))
                #oD = readOptDepthData4py("/data/rob1g10/SPH/DATA/DISKS/{}/optDepthData{}{}_theta{}_phi{}_time{:0.1f}.txt".format(baseName, CcodeName, baseName, 0, phi, time * utime / (24 * 3600)))
   #         print np.sum(np.exp(-oD))            
                times  = np.append(times, time * utime / (24*3600))
                if(np.abs(times[0]-times[-1]) % orbitalPeriod < tol):
                    end = [i, j]
                    numOrbits = int(np.abs(times[0]-times[-1]) / orbitalPeriod)
            except:
                pass
                print "../../Desktop/40d_0.0e_eqTest/optDepthData{}{}_theta{}_phi{}_time{:0.1f}.txt".format(CcodeName, baseName, 0, phi, time * utime / (24 * 3600)), "does not exist"
                
#    if(end == [0, 1]):
#        print("Tolerance too low to calculate closest whole orbit")
#        return
#    print "END = ", end
#    print "Number of Orbits = ", numOrbits 
    
    if(fold == 'n'):
        end = [QDISKnum, QDISKstepsMAX]
    
    for theta in np.arange(thetaRange[0], thetaRange[1], thetaRange[2]):
        print("THETA = ", theta)
        
        EWs = np.array([])
        times = np.array([])
#        setSize = 0
        # EW calculation        
        for i in range(end[0]+1):
            
            if(end[0] == i):
                QDISKsteps = end[1]
            else:
                QDISKsteps = QDISKstepsMAX
                
            for j in range(QDISKstepSTART, QDISKsteps):
                filename = '../../Desktop/40d_0.0e_eqTest/QGDISK{}{:02}'.format(QDISKstart+i, j)
                diskname = "../../Desktop/40d_0.0e_eqTest/disk001"
                IPname = "../../Desktop/40d_0.0e_eqTest/IPDISK" + filename[-5:-2]
                #print IPname
                #filename = '/data/rob1g10/SPH/DATA/DISKS/{}/processedDISKs/QGDISK{}{:02}'.format(baseName, QDISKstart+i, j)
                #diskname = filename[:-26] + "disk001"
                #IPname = filename[:-11] + "IPDISK" + filename[-5:-2]
#                filename = '../../Desktop/SciamaDATA/NEW/QGDISK{}{:02}'.format(411+i, j)
#                diskname = "../../Desktop/SciamaDATA/NEW/disk001"
#                IPname = "../../Desktop/SciamaDATA/NEW/IPDISK" + filename[-5:-2]
                Becoords, NScoords, partcoords, partvels, rho, u, num_parts, time, umass, udist, utime = readFILES(filename, diskname, IPname)
#                print num_parts                
                Becoords *= udist / (7 * 6.96e10)
                NScoords *= udist / (7 * 6.96e10)
                partcoords *= udist / (7 * 6.96e10)
                #R = np.sqrt(partcoords[:,0]**2 + partcoords[:,1]**2)
#                print j, time * utime / (24*3600)
#                print 'QGDISK{}{:02}'.format(QDISKstart+i, j), time * utime / (24*3600)
                #oD = 0.48 * readOptDepthData4py("/data/rob1g10/SPH/DATA/DISKS/{}/optDepthData{}{}_theta{}_phi{}_time{:0.1f}.txt".format(baseName, CcodeName, baseName, theta, phi, time * utime / (24 * 3600)))
                oD = 0.48 * readOptDepthData4py("../../Desktop/40d_0.0e_eqTest/optDepthData{}{}_theta{}_phi{}_time{:0.1f}.txt".format(CcodeName, baseName, theta, phi, time * utime / (24 * 3600)))
#                if(i == 0):
#                  oD = readOptDepthData4py("../../Desktop/SciamaDATA/NEW/ISOBEtheta{}{}{}.txt".format(411+i, j, theta))
#                else:
#                  oD = readOptDepthData4py("../../Desktop/SciamaDATA/NEW/ISOBEtheta{}{:02}{}.txt".format(411+i, j, theta))
                T = np.exp(-oD)
                #print filename, baseName, int("{:03}".format(387+i))
                EWs = np.append(EWs, np.sum(T))
                times  = np.append(times, time * utime / (24*3600))
#                setSize += 1
        EWs /= 1000
        print times[0]
        times -= times[0]
#        print "SETSIZE = ", setSize
        
        if(fold == 'n'):
            plt.plot(times, EWs, label="{} degs".format(theta))
            plt.xlabel("Time / days")
            plt.ylabel("EW")
            plt.legend(loc='best')
        
        if(fold == 'y'):
            
            orbitArrayLength = int(len(EWs) / numOrbits)
#            print orbitArrayLength, len(EWs), numOrbits
            
            # folding
            foldEWs = EWs[:orbitArrayLength]
            for j in range(1, numOrbits): 
                foldEWs += EWs[j*orbitArrayLength:(j+1)*orbitArrayLength]
            foldEWs /= numOrbits
            #plt.plot(np.linspace(0,20,32), plotEWs, label="{} degs".format(i))
            t = times[:len(foldEWs)]
            plt.plot(np.append(t, t+orbitalPeriod), 
                     np.append(foldEWs, foldEWs), label="{} degs".format(theta))
        #    plt.plot(np.append(times[0:setSize] - 9170, 20 + times[0:setSize] - 9170),
        #                       np.append(plotEWs, plotEWs), label="{} degs".format(theta))
        #    plt.plot([20 for i in range(20)], np.linspace(0, 30, 20), color='black')
        #    plt.plot([40 for i in range(20)], np.linspace(0, 30, 20), color='black')
        #    plt.xlim(0, 65)
        #    plt.ylim(5.5, 16.5)
            plt.xlabel("Time / days")
            plt.ylabel("EW")
            plt.legend(bbox_to_anchor=(1.1, 1.05))#loc='best')
        #text_y = np.array([-4, -24, -27, -30, -46.5, -47, -47.5])
        #for j in range(0,91,10):
#        EWs1 = plotEWs
#        times1 = times[0:setSize] - 9170
#        "polyfit and residuals"
        p = np.polyfit(times, EWs, 1)
        EWRESIDUALS = -EWs - (p[0]*times + p[1])
        MEAN = np.mean(-EWs)
        SD = np.std(EWRESIDUALS)
#    #    plt.text(np.min(times1) + 4, EWs1[0]+ 0.5, "mean = {:.2f}, SD = {:.2f}".format(MEAN, SD))
        print MEAN, SD
        MEANs = np.append(MEANs, MEAN)
        SDs = np.append(SDs, SD)
        
    return times, EWs

#times, EWs = foldEW(20, "20d_0.4e_eqTest", "SCIAMA", 387, 1, 66, [0, 91, 10], 0, fold='y')
#times, EWs = foldEW(20, "20d_0.2e_eqTest", "SCIAMA", 387, 1, 35, [0, 91, 10], 0, fold='y')
#times, EWs = foldEW(20, "BH20d_0.2e", "SCIAMA2", 300, 1, 39, [0, 91, 10], 0, fold='y')
#times, EWs = foldEW(20, "20d_0.0e_eqTest", "SCIAMA", 484, 1, 25, [0, 91, 10], 0, fold='y')
#times, EWs = foldEW(40, "40d_0.4e_eqTest", "SCIAMA", 380, 8, 15, [0, 91, 10], 0, fold='y')
#times, EWs = foldEW(40, "40d_0.2e_eqTest", "5e-12", 380, 4, 11, [0, 91, 10], 0, fold='y')
#times, EWs = foldEW(40, "40d_0.0e_eqTest", "5e-12", 380, 4, 10, [0, 91, 10], 0, fold='y')
#times, EWs = foldEW(60, "60d_0.4e_eqTest", "SCIAMA", 580, 6, 5, [0, 91, 10], 0, fold='y')
#times, EWs = foldEW(80, "80d_0.4e_eqTest", "SCIAMA2", 340, 1, 5, [0, 91, 10], 0, fold='n')
#times, EWs = foldEW(1, "isoBe_fromBeginningTenth", "SCIAMA3", 411, 4, 8, [0, 91, 10], 0, fold='n')
#times, EWs = foldEW(1, "isoBe_fromBeginningTenth", "1e-11", 411, 4, 8, [0, 91, 10], 0, fold='n')
#times, EWs = foldEW(1, "isoBe_fromBeginningTenth", "1e-11UP", 411, 1, 1, [0, 91, 10], 0, fold='n')
#times, EWs = foldEW(1, "isoBe_fromBeginningTenth", "1e-11DOWN", 411, 1, 1, [0, 91, 10], 0, fold='n')
#for i in range(10, times[-1], 20):
#    plt.plot(i * np.ones(10), np.linspace(0, 200, 10), color='orange')

# MEANs AND SDs ----------------------------------------------------------------

MandSD20d_04e = np.array([-11.10, 0.039, -11.05, 0.033, -10.99, 0.033, -10.91, 0.033, -10.78, 0.031, 
                          -10.66, 0.026, -10.46, 0.023, -10.07, 0.025, -9.48, 0.031, -9.25, 0.037])
MandSD20d_02e = np.array([-19.98, 0.049, -19.93, 0.027, -19.85, 0.024, -19.75, 0.044, -19.58, 0.075,
                          -19.44, 0.060, -19.18, 0.040, -18.69, 0.029, -17.72, 0.035, -17.18, 0.068])
MandSD20d_00e = np.array([-27.50, 0.017, -27.45, 0.016, -27.37, 0.028, -27.26, 0.054, -27.10, 0.072,
		                  -26.94, 0.057, -26.67, 0.041, -26.13, 0.035, -24.99, 0.040, -24.12, 0.057])
MeanIsoBe = np.array([-105.90,-104.61,-102.62,-102.61,-103.02,-103.35,-103.29,-102.72,-101.48,-94.07]) 
MeanIsoBeTemp = np.array([-86.79,-82.18,-81.64,-84.30,-86.83,-88.21,-88.10,-85.74,-75.47,-44.22])               
#plt.plot(range(0,91,10), -MandSD20d_04e[0::2], label='20d 0.4e')
#plt.plot(range(0,91,10), -MandSD20d_02e[0::2], label='20d 0.2e')
#plt.plot(range(0,91,10), -MandSD20d_00e[0::2], label='20d 0.0e')
#plt.plot(range(0,91,10), MandSD20d_04e[1::2], label='20d 0.4e')
#plt.plot(range(0,91,10), MandSD20d_02e[1::2], label='20d 0.2e')
#plt.plot(range(0,91,10), MandSD20d_00e[1::2], label='20d 0.0e')
#plt.plot(range(0,91,10), -MeanIsoBe, label='isolated Be star')
#plt.plot(range(0,91,10), -MeanIsoBeTemp, label='isolated Be star')
#plt.legend(loc='best')
#plt.ylabel('-Mean Equivalent Width')
#plt.xlabel('Inclination Angle')
#phi = 0 
#theta = 15
print calcQDISKsForOrbit(40, "40d_0.4e_eqTest", 380, 15)
# 40 = ???, [380, 383, 3], [380, 383, 6]
# 60 = [580, 585, 3], [580, 586, 3], [580, 588, 1]
# 80 = [340, 348, 1], [260, 268, 1], [260, 268, 1]
#[383, 3], [6, 12], 3
#plotEWwithTime(['/data/rob1g10/SPH/DATA/DISKS/20d_0.0e/processedDISKs/',
#                '/data/rob1g10/SPH/DATA/DISKS/20d_0.2e/processedDISKs/',
#                '/data/rob1g10/SPH/DATA/DISKS/20d_0.4e/processedDISKs/'],
#                [266, 250, 207],
#                [20, 20, 20], 
#                phi, theta)

#plotEWwithTime(['/data/rob1g10/SPH/DATA/DISKS/40d_0.0e/processedDISKs/',
#                '/data/rob1g10/SPH/DATA/DISKS/40d_0.2e/processedDISKs/',
#                '/data/rob1g10/SPH/DATA/DISKS/40d_0.4e/processedDISKs/'],
#                [110, 108, 104],
#                [20, 20, 20], 
#                phi, theta)

#plotEWwithTime(['/data/rob1g10/SPH/DATA/DISKS/60d_0.0e/processedDISKs/',
#                '/data/rob1g10/SPH/DATA/DISKS/60d_0.2e/processedDISKs/',
#                '/data/rob1g10/SPH/DATA/DISKS/60d_0.4e/processedDISKs/'],
#                [66, 66, 65],
#                [20, 20, 20], 
#                phi, theta)

#plotEWwithTime(['/data/rob1g10/SPH/DATA/DISKS/80d_0.0e/processedDISKs/',
#                '/data/rob1g10/SPH/DATA/DISKS/80d_0.2e/processedDISKs/',
#                '/data/rob1g10/SPH/DATA/DISKS/80d_0.4e/processedDISKs/'],
#                [41, 41, 40],
#                [20, 20, 20], 
#                phi, theta)

#plotEWwithTime(["justBe_halfInj2"],
#               ['/data/rob1g10/SPH/DATA/DISKS/justBe_halfInj2/processedDISKs/'],
#                [434],
#                [9], 
#                phi, theta)

#plotEWwithTime(["80d_0.0e"],
#               ['/data/rob1g10/SPH/DATA/DISKS/80d_0.0e/processedDISKs/'],
#                [61],
#                [9], 
#                phi, theta)

#print calcEWAverageWithTime('/data/rob1g10/SPH/DATA/DISKS/40d_0.0e/processedDISKs/',
#                      41, 20, phi, theta)

#plotAndSaveSpectralLine("BASIC", ['/data/rob1g10/SPH/DATA/DISKS/20d_0.0e_eqTest/processedDISKs/QGDISK38701'],
#                        "20d_0.0e_eqTest", [0, 91, 30], [0, 9, 10], 100, "uniform", "bleh")
                        
#plotAndSaveSpectralLine("SHEAR", ['/data/rob1g10/SPH/DATA/DISKS/isoBe_fromBeginningTenth/processedDISKs/QGDISK41101'],
#                        "isoBe_fromBeginningTenth", [0, 91, 30], [0, 9, 10], 100, "uniform", "bleh")
                        
#plotFeatureVels("BASIC", ['/data/rob1g10/SPH/DATA/DISKS/isoBe_fromBeginningTenth/processedDISKs/QGDISK41101'],
#                        "isoBe_fromBeginningTenth", [60, 91, 100], [0, 9, 10], 6, 100, "uniform", "bleh")

#plotFeatureVels("OUTPUT/theta1e-1150.txt", ['DATA/QGDISK41101'],
#                        "isoBe_fromBeginningTenth", [50, 91, 100], [0, 9, 10], 6, 100, "uniform", "bleh")

#for i in range(22, 57):
#    plotFeatureVels("BASIC", 
#                    ['/data/rob1g10/SPH/DATA/DISKS/20d_0.4e_eqTest/processedDISKs/QGDISK387{}'.format(i)],
#                        "20d_0.4e_eqTest", [50, 91, 100], [0, 9, 10], 6, 100, "uniform", 
#                        "Talks/Magellanic Clouds Conference/Halpha Mapping/HalphaMapping{}.png".format(i), save='y')

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
#filename = '/data/rob1g10/SPH/DATA/DISKS/justBe_halfInj/processedDISKs/QGDISK50101'
#filename = '/data/rob1g10/SPH/injFiddle/sph3D/processedDISKs/QGDISK13201'
#filename = '/data/rob1g10/SPH/DATA/DISKS/20d_0.4e/processedDISKs/QGDISK28101'
#filename = '/data/rob1g10/SPH/DATA/DISKS/80d_0.4e/processedDISKs/QGDISK06101'
#filename = '/data/rob1g10/SPH/DATA/DISKS/20d_0.4e/processedDISKs/QGDISK28301'
#filename = '/data/rob1g10/SPH/DATA/DISKS/40d_0.4e/processedDISKs/QGDISK13401'
#filename = '/data/rob1g10/SPH/injFiddle/sph3D/processedDISKs/QGDISK32701'
#filename = '/data/rob1g10/SPH/DATA/DISKS/injection_off_JJrun/processedDISKs/QGDISK13601'
#filename = '/data/rob1g10/SPH/DATA/DISKS/injection_off_a1.5e11/processedDISKs/QGDISK14701'
#filename = '/data/rob1g10/SPH/DATA/DISKS/injection_off_justBe/processedDISKs/QGDISK17001'
#filename = '/data/rob1g10/SPH/DATA/DISKS/inclined/processedDISKs/QGDISK04819'
#filename = '/data/rob1g10/SPH/DATA/DISKS/inclined/processedDISKs/QGDISK04901'
#filename = raw_input('Type the name of the ASCII file to be read')
#filename = '/data/rob1g10/SPH/DATA/DISKS/20d_0.4e/processedDISKs/QGDISK28001'

#filename = '/data/rob1g10/SPH/DATA/DISKS/isoBe_close2equilibrium/processedDISKs/QGDISK14001'
#diskname = filename[:-26] + "disk001"
#IPname = filename[:-11] + "IPDISK" + filename[-5:-2]

#filenames = ['/data/rob1g10/SPH/DATA/DISKS/justBe_halfInj2/processedDISKs/QGDISK50301',
#             '/data/rob1g10/SPH/injFiddle/sph3D/processedDISKs/QGDISK13201',
#             '/data/rob1g10/SPH/DATA/DISKS/injFiddle_largebox/processedDISKs/QGDISK32701']
#for filename in filenames:
""" Read File """
# try:
#Becoords, NScoords, partcoords, partvels, rho, u, num_parts, time, umass, udist, utime = readFILES(filename, diskname, IPname)
# except:
#     Becoords, NScoords, partcoords, partvels, rho, u, num_parts, umass, udist, utime = readFILES(filename, diskname, IPname)
#     time = 1
#Becoords2, NScoords2, partcoords2, partvels2, rho2, u2, num_parts2, time2, umass2, udist2, utime2 = readFILES(filename2)

""" binVels(partcoords, pivot, partvelocities, rhos, depthOfView, viewAngle,
            numOfGridPoints) """
#partvelsPLOT = binVels(partcoords, Becoords, partvels,
#                      1.7e6, 0, [10, 10], 2., 0.5)
#counts, bins, bars = plt.hist(partvelsPLOT, bins=100)
# R = np.sqrt(partcoords[:,0]**2 + partcoords[:,1]**2) * udist / (7 * 695700e3 * 100)
# print max(R)
# print min(R)

"""
 for just Be, units are  cm = 6439255109019.17  g = 3.583999100000000E+034  s = 334148.986076665
 frequency for H aplha is 4.57e14
 wavelength is 656.28nm
"""

#
#unitvel = udist / utime
#filenames = ['/data/rob1g10/SPH/DATA/DISKS/20d_0.4e/processedDISKs/QGDISK28101',
#             '/data/rob1g10/SPH/DATA/DISKS/40d_0.4e/processedDISKs/QGDISK13401',
#             '/data/rob1g10/SPH/DATA/DISKS/80d_0.4e/processedDISKs/QGDISK06101']

#for filename in filenames:
#    Becoords, NScoords, partcoords, partvels, rho, u, num_parts, time, umass, udist, utime = readFILES(filename)
#    name = filename[29:37]
#    plotSpectralLine(partcoords, Becoords, partvels, depth, 45,
#                         [10,10], diskMaxY, diskMaxZ, 100, udist, utime, name)
#    plt.legend()
#name = "40"
#for i in range(0, 61, 15):
#optDepth = plotSpectralLine(partcoords, Becoords, partvels, depth, 0,
#                            [10,10], diskMaxY, diskMaxZ, 100, udist, utime, "uniform",
#                            name)

#T = np.exp(-oD) / 1e-4
#print(int(max(T)), int(min(T)))
#for i in range(len(T)):
#    T[i] = T[i]
#plotAndSaveSpectralLine(partvels, T, [0,90], 100, udist, 
#                        utime, "uniform", name)
# print np.max(optDepth)
#plotHaContributions(partcoords, Becoords, partvels, 45,
#                    200, 5, 10.0, udist, utime)
#partvelsPLOT = binVels(partcoords, Becoords, partvels,
#                      depth, 45, [10, 10], diskMaxY, diskMaxZ)
#partvelsPLOT2 = binVels(partcoords2, Becoords2, partvels2,
#                      depth, 0, [10, 10], diskMaxY, diskMaxZ)

#angstroms4 = partvelsPLOT4 * unitvel / 3e8
#angstroms4 += 1
#angstroms4 *= 6562.8
#vels = partvelsPLOT * unitvel / 1e5
#vels2 = partvelsPLOT2 * unitvel / 100
#
#PLOT = np.histogram(vels, bins=100, normed='True') - np.histogram(vels, bins=100, normed='True')

#plt.hist(vels, bins=100)
#plt.xlabel('km/s')
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
#`ax3.plot(np.arange(0,200), counts3)
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

""" Peak spacing relation """

#peakSpacing = showPeakRelation4Spectra(partcoords, Becoords, partvels,
#                             depth, np.array([10, 10]), diskMaxY, diskMaxZ)

"""----------------------------------------------------------------------------

                               FIND PARTS

----------------------------------------------------------------------------"""

#partsInBox = findParts(partcoords, [-0.015, 0.01, 0.46, 0.47, -1.5, 1.5])
#print len(partsInBox)

"""
for DISK090 of John's run, elements 532 and 34847 are close and vary in
density
"""

"""----------------------------------------------------------------------------

                        ESCAPE VELOCITY STUFF

----------------------------------------------------------------------------"""

# escapeVelPlot(partcoords, partvels, 18, 25, udist, utime)
# plt.show()


"""----------------------------------------------------------------------------

                        Equivalent Width 

----------------------------------------------------------------------------"""

#filename = '/data/rob1g10/SPH/DATA/DISKS/40d_0.0e/processedDISKs/QGDISK'
#diskname = filename[:-26] + "disk001"
#IPname = filename[:-11] + "IPDISK" + filename[-5:-2]

"""-----------Relation to Orbital period and eccentricity-------------------"""

#F_c = 7.1e7 # For integration
#theta = 90
#phi = 90
#
#EWs = np.array([])
#pers = np.array([])
#Rs = np.array([])
#peakSeps = np.array([])
#
#orbits_0e = [20,40,60,80]
#DISKs_0e = ["28604", "13001", "08606", "06101"]
#orbits_2e = [20,40,60,80]
#DISKs_2e = ["27004", "12804", "08606", "06101"]
#orbits_4e = [20,40,60,80]
#DISKs_4e = ["22706", "12404", "08504", "06005"]
#
#allOrbits = orbits_0e + orbits_2e + orbits_4e
#allDISKs = DISKs_0e + DISKs_2e + DISKs_4e
#
##for theta in range(0,91,15):
#for i in range(len(allOrbits)):
#    eccentricity = 2 * (i / 4)
#    filename = '/data/rob1g10/SPH/DATA/DISKS/{}/processedDISKs/QGDISK{}'.format(str(allOrbits[i])+"d_0."+str(eccentricity)+"e", allDISKs[i])
#    diskname = filename[:-26] + "disk001"
#    IPname = filename[:-11] + "IPDISK" + filename[-5:-2]
#    Becoords, NScoords, partcoords, partvels, rho, u, num_parts, time, umass, udist, utime = readFILES(filename, diskname, IPname)
##    print filename[29:37], "------------------------"
##    print num_parts
##    print time * utime /(24 * 3600)
#    Becoords *= udist / (7 * 6.96e10)
#    NScoords *= udist / (7 * 6.96e10)
#    partcoords *= udist / (7 * 6.96e10)
#    R = np.sqrt(partcoords[:,0]**2 + partcoords[:,1]**2)
##    print max(R)
##    print(calcEW(partvels, np.ones(len(partvels)), np.ones(len(partvels))*0.01, F_c, F_c, udist, utime, intSteps=200, viewAngles=[phi,theta]))
#    baseName = filename[29:37]    
##    oD = readOptDepthData4py("/data/rob1g10/SPH/DATA/DISKS/{}/optDepthData{}_theta{}_phi{}_time{:0.1f}.txt".format(baseName, baseName, theta, phi, time * utime / (24 * 3600)))
##    T = np.exp(-oD) / 1e-3
##    EWs = np.append(EWs, calcEW(partvels, T, np.ones(len(partvels)), F_c, F_c, udist, utime, intSteps=200, viewAngles=[phi,theta]))
#    print filename, baseName, int("{:03}".format(int(allDISKs[i][:3])))
#    EWs = np.append(EWs, calcEWAverageWithTime('/data/rob1g10/SPH/DATA/DISKS/{}/processedDISKs/'.format(baseName), 
#                                               int("{:03}".format(int(allDISKs[i][:3]))), 20, phi, theta))
#    pers = np.append(pers, calcPeriastronDistance(allOrbits[i], 19.4, 0.1*eccentricity) / (7 * 6.96e8))
#    Rs = np.append(Rs, findDiskSize(partcoords))
#    [peak1, peak2] = calcPeakSeparation(partvels, np.ones(len(partvels)), [phi, theta], 200, udist, utime, "uniform")
#    peakSeps = np.append(peakSeps, np.abs(peak1 - peak2))

#plt.scatter(orbits_0e, -EWs[0:4] / 1770, label='0.0e')
#plt.scatter(orbits_2e, -EWs[4:8] / 1770, color='g', label='0.2e')
#plt.scatter(orbits_4e, -EWs[8:12] / 1770, color='r', label='0.4e')
#plt.legend()
#plt.ylim(0.0, 0.002)
#plt.ylabel("Relative EW growth")
#plt.xlabel("Orbital Period in days")
#plt.show()

"""-----------Relation to periastron distance-------------------------------"""

#plt.scatter(pers[0::4], EWs[0::4], label='20d')
#plt.scatter(pers[1::4], EWs[1::4], color='g', label='40d')
#plt.scatter(pers[2::4], EWs[2::4], color='r', label='60d')
#plt.scatter(pers[3::4], EWs[3::4], color='y', label='80d')
##plt.scatter(pers[0:4], EWs[0:4], label='0.0e')
##plt.scatter(pers[4:8], EWs[4:8], color='g', label='0.2e')
##plt.scatter(pers[8:12], EWs[8:12], color='r', label='0.4e')
#plt.legend()
#plt.ylabel("EW")
#plt.xlabel("Periastron distance in stellar radii")
#plt.gca().invert_yaxis()
#plt.show()

""" With inclinations """
 
#plt.scatter(pers[::4], EWs[::4], label='20d')
#plt.scatter(pers[1::4], EWs[1::4], color='g', label='40d')
#plt.scatter(pers[2::4], EWs[2::4], color='r', label='60d')
#plt.scatter(pers[3::4], EWs[3::4], color='y', label='80d')
#plt.legend()
#plt.ylabel("EW")
#plt.xlabel("Periastron distance in stellar radii")
#plt.gca().invert_yaxis()
#plt.show()

"""-----------Relation to peak separation-----------------------------------"""

#plt.scatter(np.log10(-EWs)[0::4], np.log10(peakSeps)[0::4], label='20d')
#plt.scatter(np.log10(-EWs)[1::4], np.log10(peakSeps)[1::4], color='g', label='40d')
#plt.scatter(np.log10(-EWs)[2::4], np.log10(peakSeps)[2::4], color='r', label='60d')
#plt.scatter(np.log10(-EWs)[3::4], np.log10(peakSeps)[3::4], color='y', label='80d')
##plt.scatter(np.log(-EWs[0:4]), np.log(peakSeps[0:4]), label='0.0e')
##plt.scatter(np.log(-EWs[4:8]), np.log(peakSeps[4:8]), color='g', label='0.2e')
##plt.scatter(np.log(-EWs[8:12]), np.log(peakSeps[8:12]), color='r', label='0.4e')
#plt.legend()
##plt.ylim(np.min(np.append(np.append(-EWs_0e,-EWs_2e),-EWs_4e))-1, 
##         np.max(np.append(np.append(-EWs_0e,-EWs_2e),-EWs_4e))+1)
#plt.xlabel("log(-EW)")
#plt.ylabel("log Peak separation in Angstroms")
#plt.show()

"""-----------Relation to disc radius---------------------------------------"""

#plt.scatter(EWs[0::4], Rs[0::4], label='20d')
#plt.scatter(EWs[1::4], Rs[1::4], color='g', label='40d')
#plt.scatter(EWs[2::4], Rs[2::4], color='r', label='60d')
#plt.scatter(EWs[3::4], Rs[3::4], color='y', label='80d')
##plt.scatter(EWs[0:4], Rs[0:4], label='0.0e')
##plt.scatter(EWs[4:8], Rs[4:8], color='g', label='0.2e')
##plt.scatter(EWs[8:12], Rs[8:12], color='r', label='0.4e')
#plt.legend()
#plt.xlabel("EW")
#plt.ylabel("Disc size in stellar radii")
#plt.gca().invert_xaxis()
#plt.show()


#theta = 90
#phi = 90
##for theta in [0, 90]:
#partx = partcoords[:,0] * np.cos(phi) - partcoords[:, 1] * np.sin(phi)
#party = (partcoords[:, 0] * np.sin(phi) + partcoords[:, 1] * np.cos(phi)) * np.cos(theta) - partcoords[:,2] * np.sin(theta)
#partz = partcoords[:, 1] * np.sin(theta) + partcoords[:, 2] * np.cos(theta)
#R = np.sqrt(partx**2 + party**2) * udist / (7 * 695700e3 * 100)
#baseName = filename[29:37]    
#oD = readOptDepthData4py("/data/rob1g10/SPH/DATA/DISKS/{}/optDepthData{}_theta{}_phi{}_time{:0.1f}.txt".format(baseName, baseName, theta, phi, time * utime / (24 * 3600)))
#T = np.exp(-oD) / 1e-3
#print(calcEW(partvels, T, R**(-1.5) / max(R**(-1.5)), F_c, F_c, udist, utime, intSteps=1000, viewAngles=[90,theta]))
#print(np.sum(T))
#count = 0
#for t in T:
#    if(t == 0):
#        count += 1
#print(count)
#blockCount = 0
#for i in range(len(T)):
#    if(partz[i] < 0 and R[i] < 1): # and T[i] == 0):
#        blockCount += 1
#print(blockCount, "/", num_parts, "=", (float(blockCount) / float(num_parts)) * 100)
#plt.hist(T, 200)
#plt.xlim(0, 1000)
#plt.ylim(0, 100)
#plt.scatter(R, T)

"""----------------------------------------------------------------------------

                        new method with isolated Be 

----------------------------------------------------------------------------"""

phi = 0

EWs1 = np.array([])
times = np.array([])
colours = ["blue", "green", "black", "red", "yellow"]
#
#CcodeNames = [ "TESTRAD", "TESTEXPH", "TESTEXP", "TESTBOTH", "TESTALL"]#, "BASICOLD_"]#["BASICRTEMP2E5"]#, 
#CcodeNames = ["TESTION10", "TESTION6", "TESTION3", "TESTION1", "TESTION06", "TESTION03", "TESTION01"]
#CcodeNames = ["BASICION6", "BASICION3", "BASICION1", "BASICION05", "BASICION01", "BASICION005", "BASICION001"] 
CcodeNames = ["5e-12", "1e-11", "5e-11", "1e-10", "BREAK"]
##CcodeNames = ["Min0.9_R8-15_", "Min0.7_R8-15_", "Min0.5_R8-15_"] #"Min0.9_R5.5-10_", "Min0.7_R5.5-10_", "Min0.5_R5.5-10_"]
#             #["Min0.9_R3-5_", "Min0.7_R3-5_", "Min0.5_R3-5_"] 
#CcodeNames = ["1e-11DOWN"]
for k in range(len(CcodeNames)):
    CcodeName = CcodeNames[k]
    for theta in range(0, 91, 10):
        for i in range(1):
            for j in range(1, 2):
                filename = 'DATA/QGDISK{}{:02}'.format(411+i, j)
                diskname = 'DATA/disk001'
                IPname = 'DATA/IPDISK411'
                #filename = '/data/rob1g10/SPH/DATA/DISKS/isoBe_fromBeginningTenth/processedDISKs/QGDISK{}{:02}'.format(411+i, j)
                #diskname = filename[:-26] + "disk001"
                #IPname = filename[:-11] + "IPDISK" + filename[-5:-2]
                Becoords, NScoords, partcoords, partvels, rho, u, num_parts, time, umass, udist, utime = readFILES(filename, diskname, IPname)
                Becoords *= udist / (7 * 6.96e10)
                NScoords *= udist / (7 * 6.96e10)
                partcoords *= udist / (7 * 6.96e10)
                R = np.sqrt(partcoords[:,0]**2 + partcoords[:,1]**2)
                baseName = "isoBe_fromBeginningTenth"    
                #oD = 0.48 * readOptDepthData4py("/data/rob1g10/SPH/DATA/DISKS/{}/optDepthData{}{}_theta{}_phi{}_time{:0.1f}.txt".format(baseName, CcodeName, baseName, theta, phi, time * utime / (24 * 3600)))
                #oD = 0.48 * readOptDepthData4py("DATA/optDepthData{}{}_theta{}_phi{}_time{:0.1f}.txt".format(CcodeName, baseName, theta, phi, time * utime / (24 * 3600)))
                oD = readOptDepthData4py("../../Desktop/SciamaDATA/NEW/theta{}{}.txt".format(CcodeName, theta)) #* 1e2
                T = np.exp(-oD)
                print max(oD)
#                Tsum = np.sum(T)
                Tsum = 0
                for j in range(len(T)):
#                    if(R[j]<100 and j != 11613 and j != 30074):
                    if(False == np.isnan(T[j])):
                        Tsum += T[j]
                #print filename, baseName, int("{:03}".format(411+i))
                if(CcodeName == "BASICOLD_"):
                    EWs1 = np.append(EWs1, Tsum/1000)
                else:
                    EWs1 = np.append(EWs1, Tsum/1000)
                times  = np.append(times, time * utime / (24*3600))
    ##EWs1 /= 1000
#    print EWs1
    x = range(0, 91, 10)
##    if(CcodeName == "Min0.75_R3-5_Steps5_"):
##        CcodeName = "Min0.95_R3-5_Steps5_"
##    if(CcodeName == "BASICOLD_"):
##        CcodeName = "BASIC"
    plt.plot(x, EWs1[k*len(x):(k+1)*len(x)], label=CcodeName, color=colours[k])
    ##plt.plot(range(0, 91, 10), EWs2, label="10")
    ##plt.plot(range(0, 91, 10), EWs3, label="100")
    ##plt.plot(range(0, 91, 10), EWs4, label="1000")
    plt.xlabel("Inclination Angle in degrees")
    plt.ylabel("-EW")
    plt.legend(loc="best")
    #plt.scatter(R, oD * 1e6, s=0.5)
    #plt.ylim(0, 10)
plt.show()
#setSize = 1#7 * 13
#for j in range(0, 4):
#    i = j * 10
#    plt.plot(times[j*setSize:(j+1)*setSize], EWs[j*setSize:(j+1)*setSize], label="{} degs".format(i))
##plt.xlim(0, 110)
#plt.xlabel("Time / days")
#plt.ylabel("EW")
#plt.legend(loc=5)
#text_y = np.array([-4, -24, -27, -30, -46.5, -47, -47.5])
#for j in range(0,91,10):
#MEANS = np.array([])
#SDs = np.array([])
#for i in range(0, 4):
#    EWs1 = EWs[i*setSize:(i+1)*setSize]
#    times1 = times[i*setSize:(i+1)*setSize] #- 0.9e4
#    "polyfit and residuals"
#    p = np.polyfit(times1, EWs1, 1)
#    EWRESIDUALS = -EWs1 - (p[0]*times1 + p[1])
#    MEAN = np.mean(-EWs1)
#    SD = np.std(EWRESIDUALS)
##    plt.text(np.min(times1) + 4, EWs1[0]+ 0.5, "mean = {:.2f}, SD = {:.2f}".format(MEAN, SD))
#    MEANS = np.append(MEANS, MEAN)
#    print MEAN, SD

#phi = 0
#
#EWs = np.array([])
#times = np.array([])
#
#for theta in range(0, 1, 10):
#    for i in range(13):
#        for j in range(1, 8):
#            if(i==4 and j==1):
#                j = 2
#            filename = '/data/rob1g10/SPH/DATA/DISKS/isoBe_fromBeginningTenth/processedDISKs/QGDISK{}{:02}'.format(411+i, j)
#            diskname = filename[:-26] + "disk001"
#            IPname = filename[:-11] + "IPDISK" + filename[-5:-2]
#            Becoords, NScoords, partcoords, partvels, rho, u, num_parts, time, umass, udist, utime = readFILES(filename, diskname, IPname)
#            Becoords *= udist / (7 * 6.96e10)
#            NScoords *= udist / (7 * 6.96e10)
#            partcoords *= udist / (7 * 6.96e10)
#            R = np.sqrt(partcoords[:,0]**2 + partcoords[:,1]**2)
#            baseName = "isoBe_fromBeginningTenth"    
#            oD = readOptDepthData4py("/data/rob1g10/SPH/DATA/DISKS/{}/optDepthData{}_theta{}_phi{}_time{:0.1f}.txt".format(baseName, baseName, theta, phi, time * utime / (24 * 3600)))
#            T = np.exp(-oD)
#            print filename, baseName, int("{:03}".format(411+i))
#            EWs = np.append(EWs, np.sum(T))
#            times  = np.append(times, time * utime / (24*3600))
#EWs /= 1000
#
#setSize = 7 * 13
#for j in range(0, 1):
#    i = j * 10
#    plt.plot(times[j*setSize:(j+1)*setSize], EWs[j*setSize:(j+1)*setSize], label="{} degs".format(i))
##plt.xlim(0, 120)
#plt.xlabel("Time / days")
#plt.ylabel("EW")
#plt.legend(loc=5)
##text_y = np.array([-4, -24, -27, -30, -46.5, -47, -47.5])
##for j in range(0,91,10):
#for i in range(0, 1):
#    EWs1 = EWs[i*setSize:(i+1)*setSize]
#    times1 = times[i*setSize:(i+1)*setSize]
#    "polyfit and residuals"
#    p = np.polyfit(times1, EWs1, 1)
#    EWRESIDUALS = -EWs1 - (p[0]*times1 + p[1])
#    MEAN = np.mean(-EWs1)
#    SD = np.std(EWRESIDUALS)
##    plt.text(np.min(times1) + 4, EWs1[0]+ 0.5, "mean = {:.2f}, SD = {:.2f}".format(MEAN, SD))
#    print MEAN, SD

# FOLDING -----------------------------------------------------------------------------------------

#phi = 0
#
#for theta in [0]:#, 10, 20, 30, 40, 50, 60, 70, 80, 90]:
#    EWs = np.array([])
#    times = np.array([])
#    for i in range(2):
#        for j in range(2, 66):
#            filename = '/data/rob1g10/SPH/DATA/DISKS/20d_0.4e_eqTest/processedDISKs/QGDISK{}{:02}'.format(387+i, j)
#            diskname = filename[:-26] + "disk001"
#            IPname = filename[:-11] + "IPDISK" + filename[-5:-2]
#            Becoords, NScoords, partcoords, partvels, rho, u, num_parts, time, umass, udist, utime = readFILES(filename, diskname, IPname)
#            Becoords *= udist / (7 * 6.96e10)
#            NScoords *= udist / (7 * 6.96e10)
#            partcoords *= udist / (7 * 6.96e10)
#            #R = np.sqrt(partcoords[:,0]**2 + partcoords[:,1]**2)
#            baseName = "20d_0.4e_eqTest"
#            oD = readOptDepthData4py("/data/rob1g10/SPH/DATA/DISKS/20d_0.4e_eqTest/optDepthData{}_theta{}_phi{}_time{:0.1f}.txt".format(baseName, theta, phi, time * utime / (24 * 3600)))
#            T = np.exp(-oD)
#            #print filename, baseName, int("{:03}".format(387+i))
#            EWs = np.append(EWs, np.sum(T))
#            times  = np.append(times, time * utime / (24*3600))
#    EWs /= 1000
#
#    setSize = 65 * 2 #- 22
#    j = 0
#    plotEWs = EWs[j*setSize:(j+1)*setSize]
##    for j in range(1, 2): 
##        plotEWs += EWs[j*setSize:(j+1)*setSize]
###    plotTimes /= 3
##    plotEWs /= 2
#    #plt.plot(np.linspace(0,20,32), plotEWs, label="{} degs".format(i))
#    plt.plot(times, EWs, label="{} degs".format(i))
##    plt.plot(np.append(times[0:setSize] - 9170, 20 + times[0:setSize] - 9170),
##                       np.append(plotEWs, plotEWs), label="{} degs".format(theta))
##    plt.plot([20 for i in range(20)], np.linspace(0, 30, 20), color='black')
##    plt.plot([40 for i in range(20)], np.linspace(0, 30, 20), color='black')
##    plt.xlim(0, 65)
##    plt.ylim(5.5, 16.5)
#    plt.xlabel("Time / days")
#    plt.ylabel("EW")
#    plt.legend(loc=5)
#    #text_y = np.array([-4, -24, -27, -30, -46.5, -47, -47.5])
#    #for j in range(0,91,10):
#    EWs1 = plotEWs
#    times1 = times[0:setSize] - 9170
#    "polyfit and residuals"
#    p = np.polyfit(times1, EWs1, 1)
#    EWRESIDUALS = -EWs1 - (p[0]*times1 + p[1])
#    MEAN = np.mean(-EWs1)
#    SD = np.std(EWRESIDUALS)
##    plt.text(np.min(times1) + 4, EWs1[0]+ 0.5, "mean = {:.2f}, SD = {:.2f}".format(MEAN, SD))
#    print MEAN, SD

#phi = 0
#
#EWs = np.array([])
#Rs = np.array([])
#peakSeps = np.array([])
#times = np.array([])
#
#theta = 90
#for i in range(100):
#    filename = '/data/rob1g10/SPH/DATA/DISKS/justBe_halfInj2/processedDISKs/QGDISK{}01'.format(325+i)
#    diskname = filename[:-26] + "disk001"
#    IPname = filename[:-11] + "IPDISK" + filename[-5:-2]
#    Becoords, NScoords, partcoords, partvels, rho, u, num_parts, time, umass, udist, utime = readFILES(filename, diskname, IPname)
##        print IPname, "------------------------"
##    print num_parts
##        print time * utime /(24 * 3600)
#    Becoords *= udist / (7 * 6.96e10)
#    NScoords *= udist / (7 * 6.96e10)
#    partcoords *= udist / (7 * 6.96e10)
#    R = np.sqrt(partcoords[:,0]**2 + partcoords[:,1]**2)
##    print max(R)
##    print(calcEW(partvels, np.ones(len(partvels)), np.ones(len(partvels))*0.01, F_c, F_c, udist, utime, intSteps=200, viewAngles=[phi,theta]))
#    baseName = "justBe_halfInj2"    
#    oD = readOptDepthData4py("/data/rob1g10/SPH/DATA/DISKS/{}/optDepthData{}_theta{}_phi{}_time{:0.1f}.txt".format(baseName, baseName, theta, phi, time * utime / (24 * 3600)))
#    T = np.exp(-oD)
##    EWs = np.append(EWs, calcEW(partvels, T, np.ones(len(partvels)), F_c, F_c, udist, utime, intSteps=200, viewAngles=[phi,theta]))
#    print filename, baseName, int("{:03}".format(325+i))
#    EWs = np.append(EWs, np.sum(T))
##    Rs = np.append(Rs, findDiskSize(partcoords))
##    [peak1, peak2] = calcPeakSeparation(partvels, np.ones(len(partvels)), [phi, theta], 200, udist, utime, "uniform")
##    peakSeps = np.append(peakSeps, np.abs(peak1 - peak2))
#    times  = np.append(times, time * utime / (24*3600))
#EWs /= 1000
#
#setSize = 100

#plt.scatter(times, -EWs, label="90 degs")
#"polyfit and residuals"
#p = np.polyfit(times, -EWs, 1)
#EWRESIDUALS = -EWs - (p[0]*times + p[1])
#MEAN = np.mean(-EWs)
#SD = np.std(EWRESIDUALS)
#plt.plot(times, EWRESIDUALS*(0.16 / 0.05), label="90 degs", color='blue')
#plt.plot(times2, EWRESIDUALS2, label="0 degs", color='g')
#plt.text(np.min(times) - 4, -47, "mean = {:.2f}, SD = {:.2f}".format(MEAN, SD))
#plt.scatter(times[setSize:2*setSize], -EWs[setSize:2*setSize], label="15 degs", color='y')
#plt.scatter(times[2*setSize:3*setSize], -EWs[2*setSize:3*setSize], label="30 degs", color='r') 
#plt.scatter(times[3*setSize:4*setSize], -EWs[3*setSize:4*setSize], label="45 degs", color='g')
#plt.scatter(times[4*setSize:5*setSize], -EWs[4*setSize:5*setSize], label="60 degs", color='black') 
#plt.scatter(times[5*setSize:6*setSize], -EWs[5*setSize:6*setSize], label="75 degs", color='orange') 
#plt.scatter(times[6*setSize:7*setSize], -EWs[6*setSize:7*setSize], label="90 degs", color='purple')  
#plt.scatter(times[5:10], -EWs[5:10], label="45 degs", color='y')
#plt.scatter(times[10:15], -EWs[10:15], label="90 degs", color='r')
#plt.legend(loc=4)
#plt.ylabel("EW")
#plt.xlabel("Time in days")
#plt.gca().invert_yaxis()

#text_y = np.array([-33, -38, -45, -46, -46.5, -47, -47.5])
#for i in range(0,7):
#    "polyfit and residuals"
#    p = np.polyfit(times[i*setSize:(i+1)*setSize], -EWs[i*setSize:(i+1)*setSize], 1)
#    EWRESIDUALS = -EWs[i*setSize:(i+1)*setSize] - (p[0]*times[i*setSize:(i+1)*setSize] + p[1])
#    MEAN = np.mean(-EWs[i*setSize:(i+1)*setSize])
#    SD = np.std(EWRESIDUALS)
#    plt.text(np.min(times) - 4, text_y[i], "mean = {:.2f}, SD = {:.2f}".format(MEAN, SD))
plt.show()

"""----------------------------------------------------------------------------

                        new method with 80d 0.0e 

----------------------------------------------------------------------------"""

#phi = 0
#theta = 90
#
#EWs = np.array([])
#Rs = np.array([])
#peakSeps = np.array([])
#times = np.array([])
#
##allDISKs = ["06101", "06102", "06103", "06104", "06105", "06106", "06107",
##            "06201", "06202", "06203", "06204", "06205", "06206", "06207",
##            "06301", "06302", "06303", "06304", "06305", "06306", "06307",
##            "06401", "06402", "06403", "06404", "06405", "06406", "06407",
##            "06501", "06502", "06503", "06504", "06505", "06506", "06507",
##            "06601", "06602", "06603", "06604", "06605", "06606", "06607",
##            "06701", "06702", "06703", "06704", "06705", "06706", "06707",
##            "06801", "06802", "06803", "06804", "06805", "06806", "06807",
##            "06901", "06902", "06903", "06904", "06905", "06906", "06907"]
#
#for baseName in ["80d_0.0e", "80d_0.2e", "80d_0.4e"]:
#    for i in range(100):
#        filename = '/data/rob1g10/SPH/DATA/DISKS/{}/processedDISKs/QGDISK{:03}01'.format(baseName, 61+i)
#        diskname = filename[:-26] + "disk001"
#        IPname = filename[:-11] + "IPDISK" + filename[-5:-2]
#        Becoords, NScoords, partcoords, partvels, rho, u, num_parts, time, umass, udist, utime = readFILES(filename, diskname, IPname)
##        print IPname, "------------------------"
##        print num_parts
##        print time * utime /(24 * 3600)
#        Becoords *= udist / (7 * 6.96e10)
#        NScoords *= udist / (7 * 6.96e10)
#        partcoords *= udist / (7 * 6.96e10)
#        R = np.sqrt(partcoords[:,0]**2 + partcoords[:,1]**2)
##        print max(R)
#    #    print(calcEW(partvels, np.ones(len(partvels)), np.ones(len(partvels))*0.01, F_c, F_c, udist, utime, intSteps=200, viewAngles=[phi,theta]))  
#        oD = readOptDepthData4py("/data/rob1g10/SPH/DATA/DISKS/{}/optDepthData{}_theta{}_phi{}_time{:0.1f}.txt".format(baseName, baseName, theta, phi, time * utime / (24 * 3600)))
#        T = np.exp(-oD)
#    #    EWs = np.append(EWs, calcEW(partvels, T, np.ones(len(partvels)), F_c, F_c, udist, utime, intSteps=200, viewAngles=[phi,theta]))
#        print filename, baseName, int("{:03}".format(61+i))
#        EWs = np.append(EWs, np.sum(T))
#        Rs = np.append(Rs, findDiskSize(partcoords))
#        [peak1, peak2] = calcPeakSeparation(partvels, np.ones(len(partvels)), [phi, theta], 200, udist, utime, "uniform")
#        peakSeps = np.append(peakSeps, np.abs(peak1 - peak2))
#        times  = np.append(times, time * utime / (24*3600))
#EWs /= 1000
#
#setSize = 100

#plt.plot(times[0:setSize], -EWs[0:setSize], label="0.0e")
#plt.plot(times[setSize:2*setSize], -EWs[setSize:2*setSize], label="0.2e", color='y')
#plt.plot(times[2*setSize:3*setSize], -EWs[2*setSize:3*setSize], label="0.4e", color='r')
#plt.scatter(times[3*setSize:4*setSize], -EWs[3*setSize:4*setSize], label="45 degs", color='g')
#plt.scatter(times[4*setSize:5*setSize], -EWs[4*setSize:5*setSize], label="60 degs", color='black') 
#plt.scatter(times[5*setSize:6*setSize], -EWs[5*setSize:6*setSize], label="75 degs", color='orange') 
#plt.scatter(times[6*setSize:7*setSize], -EWs[6*setSize:7*setSize], label="90 degs", color='purple')  
#plt.scatter(times[5:10], -EWs[5:10], label="45 degs", color='y')
#plt.scatter(times[10:15], -EWs[10:15], label="90 degs", color='r')
#plt.legend(loc=4)
#plt.ylabel("EW")
#plt.xlabel("Time in days")
#plt.gca().invert_yaxis()
#
#"lines of periastron"
#for i in range(1800, 3600, 80):
#    plt.plot(i * np.ones(10), np.linspace(np.min(-EWs), np.max(-EWs), 10), color='orange')
#
#text_y = np.array([-35, -30, -25, -46, -46.5, -47, -47.5])
#for i in range(0,3):
#    "polyfit and residuals"
#    p = np.polyfit(times[i*setSize:(i+1)*setSize], -EWs[i*setSize:(i+1)*setSize], 1)
#    EWRESIDUALS = -EWs[i*setSize:(i+1)*setSize] - (p[0]*times[i*setSize:(i+1)*setSize] + p[1])
#    MEAN = np.mean(-EWs[i*setSize:(i+1)*setSize])
#    SD = np.std(EWRESIDUALS)
#    plt.text(np.min(times) - 4, text_y[i], "For {}e: mean = {:.2f}, SD = {:.2f}".format(0.2*i, MEAN, SD))
#plt.show()

"""----------------------------------------------------------------------------

                        new method with 60d

----------------------------------------------------------------------------"""

#phi = theta = 0
#
#EWs = np.array([])
#Rs = np.array([])
#peakSeps = np.array([])
#times = np.array([])
#
#
#for baseName in ["60d_0.0e", "60d_0.2e", "60d_0.4e"]:   
#    allDISKs = np.array([])
#    if(baseName == "60d_0.4e"):
#        NUMS = [85, 94]
#    else:
#        NUMS = [86, 95]
#    for j in range(NUMS[0], NUMS[1]):
#        for k in range(1, 7):
#            NAME = "{:03}{:02}".format(j, k)
#            allDISKs = np.append(allDISKs, NAME)
#
#    for i in range(len(allDISKs)):
#        filename = '/data/rob1g10/SPH/DATA/DISKS/{}/processedDISKs/QGDISK{}'.format(baseName, allDISKs[i])
#        diskname = filename[:-26] + "disk001"
#        IPname = filename[:-11] + "IPDISK" + filename[-5:-2]
#        Becoords, NScoords, partcoords, partvels, rho, u, num_parts, time, umass, udist, utime = readFILES(filename, diskname, IPname)
##        print IPname, "------------------------"
##        print num_parts
##        print time * utime /(24 * 3600)
#        Becoords *= udist / (7 * 6.96e10)
#        NScoords *= udist / (7 * 6.96e10)
#        partcoords *= udist / (7 * 6.96e10)
#        R = np.sqrt(partcoords[:,0]**2 + partcoords[:,1]**2)
##        print max(R)
#    #    print(calcEW(partvels, np.ones(len(partvels)), np.ones(len(partvels))*0.01, F_c, F_c, udist, utime, intSteps=200, viewAngles=[phi,theta]))  
#        oD = readOptDepthData4py("/data/rob1g10/SPH/DATA/DISKS/{}/optDepthData{}_theta{}_phi{}_time{:0.1f}.txt".format(baseName, baseName, theta, phi, time * utime / (24 * 3600)))
#        T = np.exp(-oD)
#    #    EWs = np.append(EWs, calcEW(partvels, T, np.ones(len(partvels)), F_c, F_c, udist, utime, intSteps=200, viewAngles=[phi,theta]))
#        print filename, baseName, int("{:03}".format(int(allDISKs[i][:3])))
#        EWs = np.append(EWs, np.sum(T))
#        Rs = np.append(Rs, findDiskSize(partcoords))
#        [peak1, peak2] = calcPeakSeparation(partvels, np.ones(len(partvels)), [phi, theta], 200, udist, utime, "uniform")
#        peakSeps = np.append(peakSeps, np.abs(peak1 - peak2))
#        times  = np.append(times, time * utime / (24*3600))
#EWs /= 1000

#setSize = 54
#
#plt.plot(times[0:setSize], -EWs[0:setSize], label="0.0e")
#plt.plot(times[setSize:2*setSize], -EWs[setSize:2*setSize], label="0.2e", color='y')
#plt.plot(times[2*setSize:3*setSize], -EWs[2*setSize:3*setSize], label="0.4e", color='r')
#plt.legend(loc=4)
#plt.ylabel("EW")
#plt.xlabel("Time in days")
#plt.gca().invert_yaxis()
#
#"lines of periastron"
#for i in range(1740, 1950, 60):
#    plt.plot(i * np.ones(10), np.linspace(np.min(-EWs), np.max(-EWs), 10), color='orange')
#    
#text_y = np.array([-20, -19.5, -19, -46, -46.5, -47, -47.5])
#for i in range(0,3):
#    "polyfit and residuals"
#    p = np.polyfit(times[i*setSize:(i+1)*setSize], -EWs[i*setSize:(i+1)*setSize], 1)
#    EWRESIDUALS = -EWs[i*setSize:(i+1)*setSize] - (p[0]*times[i*setSize:(i+1)*setSize] + p[1])
#    MEAN = np.mean(-EWs[i*setSize:(i+1)*setSize])
#    SD = np.std(EWRESIDUALS)
#    plt.text(np.max(times) - 4, text_y[i], "For {}e: mean = {:.2f}, SD = {:.2f}".format(0.2*i, MEAN, SD))
#plt.show()

"""-----------Relation to peak separation-----------------------------------"""

#plt.scatter(peakSeps[0:5], -EWs[0:5], label="0 degs")
#plt.scatter(peakSeps[5:10], -EWs[5:10], label="45 degs", color='y')
#plt.scatter(peakSeps[10:15], -EWs[10:15], label="90 degs", color='r')
#plt.legend()
##plt.ylim(np.min(np.append(np.append(-EWs_0e,-EWs_2e),-EWs_4e))-1, 
##         np.max(np.append(np.append(-EWs_0e,-EWs_2e),-EWs_4e))+1)
#plt.xlabel("log(-EW)")
#plt.ylabel("log Peak separation in Angstroms")
#plt.show()

"""-----------Relation to disc radius---------------------------------------"""

#plt.scatter(Rs[0:5], -EWs[0:5], label="0 degs")
#plt.scatter(Rs[5:10], -EWs[5:10], label="45 degs", color='y')
#plt.scatter(Rs[10:15], -EWs[10:15], label="90 degs", color='r')
#plt.legend()
#plt.ylabel("EW")
#plt.xlabel("Disc size in stellar radii")
#plt.gca().invert_yaxis()
#plt.show()

"""----------------------------------------------------------------------------

                        EW HISTOGRAMS

----------------------------------------------------------------------------"""

#pngdirectory = "/data/rob1g10/SPH/DATA/Animations/EW contributions/isoBe5/"
#for i in range(325, 425):
#    for j in range(1, 2):
#        NAME = "{:03}{:02}".format(i, j)
#        plotEW2dHist('justBe_halfInj2', 
#                     '/data/rob1g10/SPH/DATA/DISKS/justBe_halfInj2/processedDISKs/QGDISK{}'.format(NAME), 0, 0, 
#                     show='n', save='y', outFilename='{}{}'.format(pngdirectory, NAME), cbar=[0,8],
#                     xylims=[-5,5,-5,5])

#os.system("convert -delay 20 -loop 0 {}*.png {}.gif".format(pngdirectory, "80_0.0eEW2dHist"))

#plotEW2dHist('justBe_halfInj2', 
#             '/data/rob1g10/SPH/DATA/DISKS/justBe_halfInj2/processedDISKs/QGDISK32701',
#             0, 0, show='y', save='n', outFilename='bleh', cbar=[0,8],
#             xylims=[-15,15,-15,15])
#plotEW2dHist('justBe_halfInj2', 
#             '/data/rob1g10/SPH/DATA/DISKS/justBe_halfInj2/processedDISKs/QGDISK32801',
#             0, 0, show='y', save='n', outFilename='bleh', cbar=[0,8],
#             xylims=[-15,15,-15,15])

#plotEW2dHistDiff('justBe_halfInj2', 
#             '/data/rob1g10/SPH/DATA/DISKS/justBe_halfInj2/processedDISKs/QGDISK32701', 
#             '/data/rob1g10/SPH/DATA/DISKS/justBe_halfInj2/processedDISKs/QGDISK32801', 
#             0, 0, show='y', save='n', outFilename='bleh', cbar=[0,8],
#             xylims=[-15,15,-15,15])
             
#pngdirectory = "/data/rob1g10/SPH/DATA/Animations/EWdiffs/isoBeNegative/"
#for i in range(327, 347):
#    for j in range(1, 2):
#        NAME = "{:03}{:02}".format(i, j)
#        plotEW2dHistDiff('justBe_halfInj2',
#                     '/data/rob1g10/SPH/DATA/DISKS/justBe_halfInj2/processedDISKs/QGDISK32701',
#                     '/data/rob1g10/SPH/DATA/DISKS/justBe_halfInj2/processedDISKs/QGDISK{}'.format(NAME), 0, 0, 
#                     show='n', save='y', outFilename='{}{}'.format(pngdirectory, NAME), cbar=[0,8],
#                     xylims=[-15,15,-15,15])

"""----------------------------------------------------------------------------

                        c_s vs v_r 

----------------------------------------------------------------------------"""

#def calcRadialVel(partvels, partcoords, udist, utime):
#    unitvel = 0.01 * udist / utime
#    vx = partvels[:, 0] * unitvel
#    vy = partvels[:, 1] * unitvel
#    vz = partvels[:, 2] * unitvel
#    x = partcoords[:, 0] * udist * 0.01
#    y = partcoords[:, 1] * udist * 0.01
#    z = partcoords[:, 2] * udist * 0.01
#    xradvel = (vx * x) / np.sqrt(x**2 + y**2 + z**2)
#    yradvel = (vy * y) / np.sqrt(x**2 + y**2 + z**2)
#    zradvel = (vz * z) / np.sqrt(x**2 + y**2 + z**2)
#    radvel = (xradvel + yradvel + zradvel)
#    return radvel
#    
#plt.scatter(R, calcRadialVel(partvels, partcoords, udist, utime))
def plotCswithRadius(partcoords, u, udist, utime):
    R = np.sqrt(partcoords[:,0]**2 + partcoords[:,1]**2) * udist / (7 * 695700e3 * 100)
    p = np.polyfit(np.log10(R), np.log10(np.sqrt(2 * u * (0.01 * udist / utime)**2 / 3)), 1)
    print p
    plt.plot(np.log10(R), np.log10(R)*p[0] + p[1])
    plt.scatter(np.log10(R), np.log10(np.sqrt(2 * u * (0.01 * udist / utime)**2 / 3)), color='g', s=0.5)
    plt.ylabel("log of Sound speed in m/s")
    plt.xlabel("log of Radius from Be star in stellar radii")
