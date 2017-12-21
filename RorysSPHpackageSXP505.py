import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import pylab as plt
import os
from scipy.optimize import leastsq, curve_fit
from RorysSPHpackageStripdisks import createAccretionArray
from RorysSPHpackage import readQDISK, readUnits, readFILES, readQADISK, findPeriastron, calcQDISKsForOrbit, findBinaryPhase, findDiskSize

def gridPlotForDiscSize(filenames):
    
    # find density range for colour bar
    rhoMins = np.array([])
    rhoMaxs = np.array([])
    for i in range(9):
        filename = filenames[i]
        diskname = filename[:-26] + "disk001"
        IPname = filename[:-11] + "IPDISK" + filename[-5:-2]
        Becoords, NScoords, partcoords, partvels, rho, u, num_parts, time, umass, udist, utime = readFILES(filename, diskname, IPname)            
        
        rho = rho*umass / udist**3
        
        rhoMins = np.append(rhoMins, np.min(rho))
        rhoMaxs = np.append(rhoMaxs, np.max(rho))
        
    rhomin = np.min(rhoMins)
    rhomax = np.max(rhoMaxs)
    
    f, axarr = plt.subplots(3, 3, figsize=(27, 27))    
    cm = plt.cm.get_cmap('rainbow')    
    
    Ms = [1.4, 1.4, 1.4, 
          6,   6,   6,
          10,  10,  10]
    es = [0.0, 0.2, 0.4, 
          0.0, 0.2, 0.4,
          0.0, 0.2, 0.4]


    for i in range(9):
        filename = filenames[i]
        diskname = filename[:-26] + "disk001"
        IPname = filename[:-11] + "IPDISK" + filename[-5:-2]
        Becoords, NScoords, partcoords, partvels, rho, u, num_parts, time, umass, udist, utime = readFILES(filename, diskname, IPname)            
        
        partcoords *= udist / (7*7e10)
        NScoords *= udist / (7*7e10)
        rho = rho*umass / udist**3
        

        sc = axarr[i/3, i%3].scatter(partcoords[:,0], partcoords[:,1], 
                              c=rho, s = 1, vmin=rhomin, vmax=rhomax, cmap=cm, edgecolor='')
#        axarr[i/3, i%3].scatter(NScoords[0], NScoords[1], color='b')
        axarr[i/3, i%3].set_xlim(-10, 10)
        axarr[i/3, i%3].set_ylim(-10, 10)
#        axarr[i/3, i%3].set_xlabel(r"R / Be star radius")
        axarr[i/3, i%3].set_ylabel(r"R / Be star radius") 
        axarr[i/3, i%3].set_title(r"M{}, {}e".format(Ms[i], es[i]))
        axarr[i/3, i%3].set(aspect='equal')
                     
#    cb = plt.colorbar(sc, ax=axarr[0, 0])
#    cb.ax.set_ylabel("Density in g per cm cubed")
                     

def findMinRelativeVelocity(baseNames, Ms, orbits):
    binaryPhase = np.array([])
    for k in range(len(baseNames)):
        closestVel = np.array([])
        bPhases = np.array([])
        filenames = np.array([])
        filename = '/data/rob1g10/SPH/DATA/DISKS/'+baseNames[k]
        
        os.system("ls /data/rob1g10/SPH/DATA/DISKS/{}/processedDISKs/QGDISK* > lsDISKSIZE.txt".format(baseNames[k]))
    
        f = open("lsDISKSIZE.txt", 'r')
        lines = f.readlines()
    	
        DISKstart = int(lines[-1][-6:-3])-9
        DISKsteps = int(lines[-1][-3:])
    
        f.close()
        
        DISKend, STEPend = calcQDISKsForOrbit(orbit[k], baseNames[k], DISKstart, DISKsteps)
        
        for i in range(DISKstart, DISKend+1):
            if(i == DISKend):
                end = STEPend
            else:
                end = DISKsteps
            for j in range(1, end):
                filenames = np.append(filenames, '/data/rob1g10/SPH/DATA/DISKS/{}/processedDISKs/QGDISK{:03}{:02}'.format(baseNames[k], i, j))
                diskname = filenames[-1][:-26] + "disk001"
                IPname = filenames[-1][:-11] + "IPDISK" + filenames[-1][-5:-2]
                Becoords, NScoords, partcoords, partvels, rho, u, num_parts, time, umass, udist, utime = readFILES(filenames[-1], diskname, IPname)            
                NSvels = partvels[1]
                Bevels = partvels[0]
#                partvels = np.delete(partvels, [0,1])
#                # find closest particle THAT IS BELOW ESCAPE VELOCITY
#                for pc in range(len(partcoords)):
#                    escapeVel = np.sqrt(2 * 6.67e-11 * Ms[k] / np.sqrt(np.sum( (partcoords[pc]-Becoords)**2)))
#                    if( np.sqrt(np.sum(partvels**2)) > escapeVel):
#                        partcoords = np.delete(partcoords, pc)
#                        partvels = np.delete(partvels, pc)
#                        
#                closestParticle = np.argmin( np.sqrt(np.sum( (partcoords-NScoords)**2)))
#                    
#                closestVel = np.append(closestVel, np.sqrt(np.sum( (NSvels-partvels[closestParticle])**2 )))
                periastrons = np.array([10604.8956481,10876.6435185,9724.82893519,19263.4161806,10703.5259722,22310.5965972,
                                14831.2488889,9698.40474537])
    
                days = (time*utime / (24*3600)) - periastrons[k]
                if( (days % orbits[k])/orbits[k] < 0.5):
                    closestVel = np.append(closestVel, np.sqrt(np.sum( (NSvels-Bevels)**2 )))
                    bPhases = np.append(bPhases, ( days % orbits[k])/orbits[k])
        
        binaryPhase = np.append(binaryPhase, bPhases[np.argmin(closestVel)])
        print baseNames[k], binaryPhase[-1]
        print days
    return binaryPhase


def plotCOaccretion(baseNames, Mx, orbit, colour):
    
    #--------------------------------------------------------------------------
    # FIND BINARY PHASE OF EACH TIMESTEP
    #--------------------------------------------------------------------------
    
    
    Lmax = np.array([])
    Lave = np.array([])    
    maxMaxRhos = np.array([])
    aveMinU = np.array([])
    aveMinV = np.array([])
    aveDiscSize = np.array([])
    
    for BN in range(len(baseNames)):
        NSx = np.array([])
        filenames = np.array([])
        times = np.array([])
        relVel = np.array([])
        Aparts = np.array([])
        maxRhos = np.array([])
        minU = np.array([])
        minV = np.array([])
        discSize = np.array([])
        
        os.system("ls /data/rob1g10/SPH/DATA/DISKS/{}/processedDISKs/QGDISK* > lsDISKSIZE.txt".format(baseNames[BN]))
    
        f = open("lsDISKSIZE.txt", 'r')
        lines = f.readlines()
    	
        DISKstart = int(lines[-1][-6:-3])-9
        DISKsteps = int(lines[-1][-3:])
    
        f.close()
        
        DISKend, STEPend = calcQDISKsForOrbit(orbit[BN], baseNames[BN], DISKstart, DISKsteps)
        
        for i in range(DISKstart, DISKend+1):
            if(i == DISKend):
                end = STEPend
            else:
                end = DISKsteps
            for j in range(1,end):
                filenames = np.append(filenames, '/data/rob1g10/SPH/DATA/DISKS/{}/processedDISKs/QGDISK{:03}{:02}'.format(baseNames[BN], i, j))
                diskname = filenames[-1][:-26] + "disk001"
                IPname = filenames[-1][:-11] + "IPDISK" + filenames[-1][-5:-2]
                Becoords, NScoords, partcoords, partvels, rho, u, num_parts, time, umass, udist, utime = readFILES(filenames[-1], diskname, IPname) 
    
                Acoords = readQADISK('/data/rob1g10/SPH/DATA/DISKS/{}/processedDISKs/QADISK{:03}{:02}'.format(baseNames[BN], i, j))
    
                NSx = np.append(NSx, NScoords[0])
                times = np.append(times, time )
                relVel = np.append(relVel, np.sum( (partvels[0]-partvels[1])**2 )  )
                Aparts = np.append(Aparts, len(Acoords))
                maxRhos = np.append(maxRhos, np.mean(np.sort(rho)))#*umass / udist**3)
                minU = np.append(minU, np.min(u * (udist * 0.01)**2/utime**2))
                minV = np.append(minV, np.min(np.sum(partvels**2, 1)) * udist * 0.01 / utime)
                discSize = np.append(discSize, findDiskSize([filenames[-1]]))
        
        days = times* utime / (24*3600)
        periastron = days[np.argmin(NSx)]
        binaryPhase = ( (periastron - days)  ) / orbit[BN]
#        print times, times*utime / (24*3600)
        
    #--------------------------------------------------------------------------
    # FIND NO. CAPTURED PARTICLES AT EACH TIMESTEP
    #--------------------------------------------------------------------------
        
        os.system("ls /data/rob1g10/SPH/DATA/DISKS/{}/disk* > lsdisks.txt".format(baseNames[BN]))
    
        f = open("lsdisks.txt", 'r')
        lines = f.readlines()
    	
        numdisks = int(lines[-1][-3:])
        
        f.close()

        disks = np.array([])        
        
        for ND in range(1, numdisks+1):
            disks = np.append(disks, '/data/rob1g10/SPH/DATA/DISKS/{}/disk{:03}'.format(baseNames[BN], ND))
            
        accTimes, accArray = createAccretionArray(disks,"NS", None, 'day', tol=0.01) 
#        accArray *= 1e-15 * umass / (24*3600)
        accTimes = accTimes.astype(float)
        
        timeStart = np.argmin(  np.abs(accTimes-times[0]  )  )
        timeEnd = np.argmin(  np.abs(accTimes-times[-1] )  )
        
#        print times[0], times[-1]
        print timeStart, timeEnd
        
        accPhase = np.linspace(binaryPhase[0], binaryPhase[-1], len(accTimes[timeStart:timeEnd]) )
        accPhase = accPhase % 1
        
#        accTimes =  * utime / (24*3600)
#        accPhase = (accTimes % orbit[BN] ) / orbit[BN]
        plotAcc = accArray[timeStart:timeEnd]
        
#        plotAcc = np.array([])
#        for T in range(len(times)):
#            plotAcc = np.append(plotAcc, accArray[  np.argmin(  np.abs(accTimes-times[T])  )  ])
        
#        print len(accArray), accTimes[-10:-1].astype(float)* utime / (24*3600), days, disks
        
        
#        plt.plot(accPhase, plotAcc, color=colour[BN], label='{}'.format(Mx[BN]))
#        plt.plot(range(len(accPhase)), accPhase, color=colour[BN], label='{}'.format(Mx[BN]))
#        plt.scatter(range(len(binaryPhase)), binaryPhase, color=colour[BN], label='{}'.format(Mx[BN]))
        
        Lmax = np.append(Lmax, np.max(plotAcc))
        Lave = np.append(Lave, np.mean(plotAcc))
        maxMaxRhos = np.append(maxMaxRhos, np.max(maxRhos))
        aveMinU = np.append(aveMinU, np.mean(minU) )
        aveMinV = np.append(aveMinV, np.mean(minV))
        aveDiscSize = np.append(aveDiscSize, np.mean(discSize))
        print baseNames[BN]
    plt.legend(loc='best')
#    plt.ylim(0, 3.5e16)
    return Lmax, Lave, maxMaxRhos, aveMinU, aveMinV, aveDiscSize
        
baseNames=np.array([

                    'M1.4_20d_0.0e', 
                    'M3_19d_0.0e',   
                    'M4_18.4d_0.0e', 
                    'M5_17.9d_0.0e', 
                    'M6_17.4d_0.0e', 
                    'M7_17d_0.0e',   
                    'M8_16.6d_0.0e',
                    'M9_16.2d_0.0e', 
                    'M10_15.8d_0.0e', 
                    
                    'M1.4_23.4d_0.1e', 
                    'M3_22.2d_0.1e',   
                    'M4_21.6d_0.1e', 
                    'M5_21d_0.1e', 
                    'M6_20.4d_0.1e', 
                    'M7_19.9d_0.1e',   
                    'M8_19.4d_0.1e',
#                    'M9_19d_0.1e', 
#                    'M10_18.5d_0.1e',                     
                    
                    'M1.4_28d_0.2e', 
                    'M3_26.5d_0.2e', 
                    'M4_25.7d_0.2e', 
                    'M5_25d_0.2e',   
                    'M6_24.3d_0.2e', 
                    'M7_23.7d_0.2e', 
                    'M8_23.1d_0.2e',
                    'M9_22.6d_0.2e',  
                    'M10_22.1d_0.2e',
#                    
#                    'M1.4_43d_0.4e', 
#                    'M3_40.8d_0.4e', 
#                    'M4_39.6d_0.4e', 
#                    'M5_38.5d_0.4e', 
#                    'M6_37.5d_0.4e', 
#                    'M7_36.5d_0.4e', 
#                    'M8_35.6d_0.4e',
#                    'M9_34.8d_0.4e',  
#                    'M10_34.1d_0.4e',
                    
#                    'M1.4_79.1d_0.6e', 
#                    'M3_75d_0.6e', 
#                    'M4_72.8d_0.6e', 
#                    'M5_70.7d_0.6e', 
#                    'M6_68.8d_0.6e', 
#                    'M7_67.1d_0.6e', 
#                    'M8_65.5d_0.6e', 
#                    'M9_64d_0.6e', 
#                    'M10_62.6d_0.6e'
                    
                    
                    ])

#Mx = 2e30*np.array([1.4, 3., 4., 5., 6., 7., 8., 9., 10.,
#                    1.4, 3., 4., 5., 6., 7., 8., 9., 10.,
#                    1.4, 3., 4., 5., 6., 7., 8., 9., 10.,
#                    1.4, 3., 4., 5., 6., 7., 8., 9., 10.])
#orbit = np.array([20.0, 19.0, 18.4, 17.9, 17.4, 17.0, 16.6, 15.8,
#                  28.0, 26.5, 25.7, 25.0, 24.3, 23.7, 23.1, 22.1,
#                  43.0, 40.8, 39.6, 38.5, 37.5, 36.5, 35.6, 34.1])
orbit = np.array([])
Mx = np.array([])
for BN in baseNames:
    underscores = np.array([])
    index = 0
    for char in BN:
#        print char, index
        if(char == 'M'):
            Mplace = index
        if(char == 'd'):
            dplace = index
        if(char == '_'):
            underscores = np.append(underscores, index)
        index += 1
    underscore1 = int(np.min(underscores))
#    print underscore1
    orbit = np.append(orbit, float(BN[underscore1+1:dplace]) )
    Mx = np.append(Mx, float(BN[Mplace+1:underscore1]) )

#print Mx
#print orbit,

bPhases04e = np.array([0.11618648, 0.5976398, 0.38550201, 0.83917929, 0.93665309, 0.75851662, 0.83093399, 0.91999362])

#for k in range(16, len(baseNames)):
#    filenames = np.array([])
#    filename = '/data/rob1g10/SPH/DATA/DISKS/'+baseNames[k]
#    
#    os.system("ls /data/rob1g10/SPH/DATA/DISKS/{}/processedDISKs/QGDISK* > lsDISKSIZE.txt".format(baseNames[k]))
#
#    f = open("lsDISKSIZE.txt", 'r')
#    lines = f.readlines()
#	
#    DISKstart = int(lines[-1][-6:-3])-9
#    DISKsteps = int(lines[-1][-3:])
#
#    f.close()
#    
#    DISKend, STEPend = calcQDISKsForOrbit(orbit[k], baseNames[k], DISKstart, DISKsteps)
#    
##    if(k>0):
##        for num in range(9):
##            DISKend, STEPend = calcQDISKsForOrbit(orbit[k], baseName[k], DISKend, DISKsteps)
##    print DISKend - DISKstart
#    for i in range(DISKstart, DISKend+1):
#        if(i == DISKstart):
#            start = 1#STEPstart[k]
#        else:
#            start = 1
#        if(i == DISKend):
#            end = STEPend
#        else:
#            end = DISKsteps
#        for j in range(start, end):
#            filenames = np.append(filenames, '/data/rob1g10/SPH/DATA/DISKS/{}/processedDISKs/QGDISK{:03}{:02}'.format(baseNames[k], i, j))
#            diskname = filenames[-1][:-26] + "disk001"
#            IPname = filenames[-1][:-11] + "IPDISK" + filenames[-1][-5:-2]
#            Becoords, NScoords, partcoords, partvels, rho, u, num_parts, time, umass, udist, utime = readFILES(filenames[-1], diskname, IPname)            
##            rho0s = np.append(rho0s, np.max(rho)*umass/udist**3)
#    per, perFile = findPeriastron(filenames)
##    print(per*udist/(7*7e10), perFile, baseNames[k])
#    diskname = perFile[:-26] + "disk001"
#    IPname = perFile[:-11] + "IPDISK" + filenames[-1][-5:-2]
#    Becoords, NScoords, partcoords, partvels, rho, u, num_parts, time, umass, udist, utime = readFILES(perFile, diskname, IPname)
#    print (time*utime/(24*3600))
    
#Periastron
#gridPlotForDiscSize(['/data/rob1g10/SPH/DATA/DISKS/M1.4_20d_0.0e/processedDISKs/QGDISK13457',
#                     '/data/rob1g10/SPH/DATA/DISKS/M1.4_28d_0.2e/processedDISKs/QGDISK30902',
#                     '/data/rob1g10/SPH/DATA/DISKS/M1.4_43d_0.4e/processedDISKs/QGDISK51205',
#                     
#                     '/data/rob1g10/SPH/DATA/DISKS/M6_17.4d_0.0e/processedDISKs/QGDISK39517',
#                     '/data/rob1g10/SPH/DATA/DISKS/M6_24.3d_0.2e/processedDISKs/QGDISK36025',
#                     '/data/rob1g10/SPH/DATA/DISKS/M6_37.5d_0.4e/processedDISKs/QGDISK15820',
#                     
#                     '/data/rob1g10/SPH/DATA/DISKS/M10_15.8d_0.0e/processedDISKs/QGDISK32729',
#                     '/data/rob1g10/SPH/DATA/DISKS/M10_22.1d_0.2e/processedDISKs/QGDISK37220',
#                     '/data/rob1g10/SPH/DATA/DISKS/M10_34.1d_0.4e/processedDISKs/QGDISK14102'])
                     
#Apastron
#gridPlotForDiscSize(['/data/rob1g10/SPH/DATA/DISKS/M1.4_20d_0.0e/processedDISKs/QGDISK13515',
#                     '/data/rob1g10/SPH/DATA/DISKS/M1.4_28d_0.2e/processedDISKs/QGDISK30918',
#                     '/data/rob1g10/SPH/DATA/DISKS/M1.4_43d_0.4e/processedDISKs/QGDISK51310',
#                     
#                     '/data/rob1g10/SPH/DATA/DISKS/M6_17.4d_0.0e/processedDISKs/QGDISK39533',
#                     '/data/rob1g10/SPH/DATA/DISKS/M6_24.3d_0.2e/processedDISKs/QGDISK36042',
#                     '/data/rob1g10/SPH/DATA/DISKS/M6_37.5d_0.4e/processedDISKs/QGDISK15836',
#                     
#                     '/data/rob1g10/SPH/DATA/DISKS/M10_15.8d_0.0e/processedDISKs/QGDISK32745',
#                     '/data/rob1g10/SPH/DATA/DISKS/M10_22.1d_0.2e/processedDISKs/QGDISK37237',
#                     '/data/rob1g10/SPH/DATA/DISKS/M10_34.1d_0.4e/processedDISKs/QGDISK14118'])
    
#print findMinRelativeVelocity(baseNames[16:24], Mx[16:24], orbit[16:24])
#for i in range(16,24):
#    phases = findBinaryPhase(baseNames[i], orbit[i] , colour = ['b', 'g', 'r', 'y', 'black', 'brown', 'turquoise', 'pink'][i-16])
#    print phases
    
#Lmax, Lave, maxRhos, minU, minV, discSize = plotCOaccretion(baseNames, Mx, orbit, colour = ['b', 'g', 'r', 'y', 'black', 'brown', 'turquoise', 'pink'])

sec1 = 0
sec2 = sec1 + 9
sec3 = sec2 + 7
sec4 = sec3 + 9


#plt.scatter(Mx[sec4:], np.log10(Lave[sec4:]), marker='H', color='black', label="e = 0.6")
#plt.scatter(Mx[sec3:sec4], np.log10(Lave[sec3:sec4]), marker='^', color='g', label="e = 0.4")
#plt.scatter(Mx[sec2:sec3], np.log10(Lave[sec2:sec3]), marker='x', color='r', label="e = 0.2")
#plt.scatter(Mx[:sec2], np.log10(Lave[:sec2]), marker='o', color='b', label="e = 0.0")
#plt.ylabel("log Average X-ray Luminosity")
#plt.xlabel("Compact object mass")

#plt.scatter(Mx[sec4:], np.log10(Lmax[sec4:]), marker='H', color='black', label="e = 0.6")
#plt.scatter(Mx[sec3:sec4], np.log10(Lmax[sec3:sec4]), marker='^', color='g', label="e = 0.4")
#plt.scatter(Mx[sec2:sec3], np.log10(Lmax[sec2:sec3]), marker='x', color='r', label="e = 0.2")
#plt.scatter(Mx[:sec2], np.log10(Lmax[:sec2]), marker='o', color='b', label="e = 0.0")
#plt.ylabel("log Maximum X-ray Luminosity")
#plt.xlabel("Compact object mass")

#plt.scatter(Mx[sec4:], np.log10(maxRhos[sec4:]), marker='H', color='black', label="e = 0.6")
#p1 = np.polyfit(Mx[sec3:sec4][0:5], np.log10(maxRhos[sec3:sec4][0:5]), 1)
#print p1[0], p1[1]
#plt.plot(Mx[sec3:sec4][0:5], Mx[sec3:sec4][0:5]*p1[0] + p1[1], color='g')
#p1 = np.polyfit(Mx[sec3:sec4][4:], np.log10(maxRhos[sec3:sec4])[4:], 1)
#print p1[0], p1[1]
#plt.plot(Mx[sec3:sec4][4:], Mx[sec3:sec4][4:]*p1[0] + p1[1], color='g')
#plt.scatter(Mx[sec3:sec4], np.log10(maxRhos[sec3:sec4]), marker='^', color='g', label="e = 0.4")
#p2 = np.polyfit(Mx[sec2:sec3], np.log10(maxRhos[sec2:sec3]), 1)
#print p2[0], p2[1]
#plt.plot(Mx[sec2:sec3], Mx[sec2:sec3]*p2[0] + p2[1], color='r')
#plt.scatter(Mx[sec2:sec3], np.log10(maxRhos[sec2:sec3]), marker='x', color='r', label="e = 0.2")
#p3 = np.polyfit(Mx[:sec2], np.log10(maxRhos[:sec2]), 1)
#print p3[0], p3[1]
#plt.plot(Mx[:sec2], Mx[:sec2]*p3[0] + p3[1], color='b')
#plt.scatter(Mx[:sec2], np.log10(maxRhos[:sec2]), marker='o', color='b', label="e = 0.0")
#plt.ylabel("log base gas density")
#plt.xlabel("Compact object mass")

##c_s = np.sqrt(2 * minU / 3)
#BR = 2*6.67e-11*Mx*2e30 / minV**2
#plt.scatter(Mx[sec4:], np.log10(BR[sec4:] / (7*7e10)), marker='H', color='black', label="e = 0.6")
#plt.scatter(Mx[sec3:sec4], np.log10(BR[sec3:sec4] / (7*7e10)), marker='^', color='g', label="e = 0.4")
#plt.scatter(Mx[sec2:sec3], np.log10(BR[sec2:sec3] / (7*7e10)), marker='x', color='r', label="e = 0.2")
#plt.scatter(Mx[:sec2], np.log10(BR[:sec2] / (7*7e10)), marker='o', color='b', label="e = 0.0")
#plt.ylabel("log Bondi accretion radius / stellar radii")
#plt.xlabel("Compact object mass")

#q = 13. / Mx
#rocheLobe = ( 0.49*q**(2./3.) ) / ( 0.6*q**(2./3.) + np.log(1. + q**(1./3.)) )
##rocheLobeSimple = 0.38 + 0.2*np.log10(q) 
##plt.scatter(Mx[sec4:], np.log10(rocheLobe[sec4:]), marker='H', color='black', label="e = 0.6")
##plt.scatter(Mx[sec3:sec4], np.log10(rocheLobe[sec3:sec4]), marker='^', color='g', label="e = 0.4")
##plt.scatter(Mx[sec2:sec3], np.log10(rocheLobe[sec2:sec3]), marker='x', color='r', label="e = 0.2")
#plt.scatter(Mx[:sec2], rocheLobe[:sec2]*10.9, marker='o', color='b')
##plt.scatter(Mx[:sec2], rocheLobeSimple[:sec2], marker='o', color='b')
#plt.ylabel("Roche lobe radius / stellar radii")
#plt.xlabel("Compact object mass")

#plt.scatter(Mx[sec4:], discSize[sec4:], marker='H', color='black', label="e = 0.6")
p1 = np.polyfit(Mx[sec3:sec4], discSize[sec3:sec4], 1)
print p1[0], p1[1]
plt.plot(Mx[sec3:sec4], Mx[sec3:sec4]*p1[0] + p1[1], color='g')
plt.scatter(Mx[sec3:sec4], discSize[sec3:sec4], marker='^', color='g', label="e = 0.4")
p2 = np.polyfit(Mx[sec2:sec3][:-1], discSize[sec2:sec3][:-1], 1)
print p2[0], p2[1]
plt.plot(Mx[sec2:sec3], Mx[sec2:sec3]*p2[0] + p2[1], color='r')
plt.scatter(Mx[sec2:sec3], discSize[sec2:sec3], marker='x', color='r', label="e = 0.2")
p3 = np.polyfit(Mx[:sec2], discSize[:sec2], 1)
print p3[0], p3[1]
plt.plot(Mx[:sec2], Mx[:sec2]*p3[0] + p3[1], color='b')
plt.scatter(Mx[:sec2], discSize[:sec2], marker='o', color='b', label="e = 0.0")
plt.ylabel("Disc size / stellar radii")
plt.xlabel("Compact object mass")
plt.legend(loc='lower left')
