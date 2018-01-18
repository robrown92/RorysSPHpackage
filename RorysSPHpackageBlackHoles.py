import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import pylab as plt
from matplotlib.legend_handler import HandlerErrorbar
import os
from scipy.optimize import leastsq, curve_fit
from RorysSPHpackageStripdisks import createAccretionArray
from RorysSPHpackage import readQDISK, readUnits, readFILES, readQADISK, findPeriastron, calcQDISKsForOrbit, findBinaryPhase, findDiskSize

def gridPlotForDiscSizeSIMPLE(filenames):
    
    orbit = np.array([])
    Ms = np.array([])
    es = np.array([])
    for FN in filenames:
        BN = FN[29:][:-27]
        print BN
        underscores = np.array([])
        index = 0
        for char in BN:
            if(char == 'M'):
                Mplace = index
            if(char == 'd'):
                dplace = index
            if(char == '_'):
                underscores = np.append(underscores, index)
            index += 1
        underscore1 = int(np.min(underscores))

        orbit = np.append(orbit, float(BN[underscore1+1:dplace]) )
        Ms = np.append(Ms, float(BN[Mplace+1:underscore1]) )
        es = np.append(es, float(BN[dplace+2:dplace+5]) )
    
    # find density range for colour bar
    rhoMins = np.array([])
    rhoMaxs = np.array([])
    NSx = np.array([])
    for i in range(len(filenames)):
        filename = filenames[i]
        diskname = filename[:-26] + "disk001"
        IPname = filename[:-11] + "IPDISK" + filename[-5:-2]
        Becoords, NScoords, partcoords, partvels, rho, u, num_parts, time, umass, udist, utime = readFILES(filename, diskname, IPname)            
        
        rho = rho*umass / udist**3
        
        rhoMins = np.append(rhoMins, np.min(rho))
        rhoMaxs = np.append(rhoMaxs, np.max(rho))
        NSx = np.append(NSx, NScoords[0])
        
    rhomin = np.min(np.abs(rhoMins))
    rhomax = np.max(rhoMaxs)
    
    f, axarr = plt.subplots(3, 4, figsize=(27, 36))    
    cm = plt.cm.get_cmap('rainbow')    

    for i in range(len(filenames)):
        filename = filenames[i]
        diskname = filename[:-26] + "disk001"
        IPname = filename[:-11] + "IPDISK" + filename[-5:-2]
        Becoords, NScoords, partcoords, partvels, rho, u, num_parts, time, umass, udist, utime = readFILES(filename, diskname, IPname)            
        
        partcoords *= udist / (7*7e10)
        NScoords *= udist / (7*7e10)
        rho = rho*umass / udist**3
        
        partcoords *= 7

        sc = axarr[i/4, i%4].scatter(partcoords[:,0], partcoords[:,1], 
                              c=rho, s = 1, vmin=rhomin, vmax=rhomax, cmap=cm, edgecolor='')
#        axarr[i/4, i%4].scatter(NScoords[0], NScoords[1], color='b')
        axarr[i/4, i%4].set_xlim(-80, 80)
        axarr[i/4, i%4].set_ylim(-80, 80)
#        axarr[i/3, i%3].set_xlabel(r"R / Be star radius")
#        axarr[i/4, i%4].set_ylabel(r"R (R$_\odot$)") 
#        axarr[i/4, i%4].set_title(r"{}M$_\odot$, e={}".format(Ms[i], es[i]))
        axarr[i/4, i%4].text(-78, 65., r"{}M$_\odot$, e={}".format(Ms[i], es[i]))
        axarr[i/4, i%4].set(aspect='equal')
    f.text(0.5, 0.04, r"X (R$_\odot$)", ha='center', fontsize=15)
    f.text(0.1, 0.5, r"Y (R$_\odot$)", va='center', rotation='vertical', fontsize=15)
#    cb = plt.colorbar(sc, ax=axarr[0, 0])
#    cb.ax.set_ylabel("Density in g per cm cubed")

def gridPlotForDiscSize(baseNames, periastron='n'):
    
    filenames = np.array([])
    filenamesSEARCH = np.array([])
    orbit = np.array([])
    Ms = np.array([])
    es = np.array([])
    for BN in baseNames:
        underscores = np.array([])
        index = 0
        for char in BN:
#            print char, index
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
        Ms = np.append(Ms, float(BN[Mplace+1:underscore1]) )
        es = np.append(es, float(BN[dplace+2:dplace+5]) )
        print BN
        
        # determine periastron/apastron
        os.system("ls /data/rob1g10/SPH/DATA/DISKS/{}/processedDISKs/QGDISK* > lsDISKSIZE.txt".format(BN))
        
        f = open("lsDISKSIZE.txt", 'r')
        lines = f.readlines()
    	
        DISKend = int(lines[-1][-6:-3])
        DISKstart = int(lines[0][-6:-3])
        DISKsteps = int(lines[-1][-3:])
    
        f.close()
        
        NSx = np.array([])
        for DISK in range(DISKstart, DISKend):
            for DISKstep in range(1, DISKsteps):
                filenamesSEARCH = np.append(filenamesSEARCH, '/data/rob1g10/SPH/DATA/DISKS/{}/processedDISKs/QGDISK{:03}{:02}'.format(BN, DISK, DISKstep))
                diskname = filenamesSEARCH[-1][:-26] + "disk001"
                IPname = filenamesSEARCH[-1][:-11] + "IPDISK" + filenamesSEARCH[-1][-5:-2]
                Becoords, NScoords, partcoords, partvels, rho, u, num_parts, time, umass, udist, utime = readFILES(filenamesSEARCH[-1], diskname, IPname) 

                NSx = np.append(NSx, NScoords[0])
        if(periastron=='y'):       
            filenames = np.append(filenames, filenamesSEARCH[np.argmin(NSx)])
        else:       
            filenames = np.append(filenames, filenamesSEARCH[np.argmax(NSx)])
            
    
    
    
    # find density range for colour bar
    rhoMins = np.array([])
    rhoMaxs = np.array([])
    NSx = np.array([])
    for i in range(len(filenames)):
        filename = filenames[i]
        diskname = filename[:-26] + "disk001"
        IPname = filename[:-11] + "IPDISK" + filename[-5:-2]
        Becoords, NScoords, partcoords, partvels, rho, u, num_parts, time, umass, udist, utime = readFILES(filename, diskname, IPname)            
        
        rho = rho*umass / udist**3
        
        rhoMins = np.append(rhoMins, np.min(rho))
        rhoMaxs = np.append(rhoMaxs, np.max(rho))
        NSx = np.append(NSx, NScoords[0])
        
    rhomin = np.min(rhoMins)
    rhomax = np.max(rhoMaxs)
    
    
    f, axarr = plt.subplots(3, 4, figsize=(27, 36))    
    cm = plt.cm.get_cmap('rainbow')    


    for i in range(len(filenames)):
        filename = filenames[i]
        diskname = filename[:-26] + "disk001"
        IPname = filename[:-11] + "IPDISK" + filename[-5:-2]
        Becoords, NScoords, partcoords, partvels, rho, u, num_parts, time, umass, udist, utime = readFILES(filename, diskname, IPname)            
        
        partcoords *= udist / (7*7e10)
        NScoords *= udist / (7*7e10)
        rho = rho*umass / udist**3
        

        sc = axarr[i/4, i%4].scatter(partcoords[:,0], partcoords[:,1], 
                              c=rho, s = 1, vmin=rhomin, vmax=rhomax, cmap=cm, edgecolor='')
#        axarr[i/3, i%3].scatter(NScoords[0], NScoords[1], color='b')
        axarr[i/4, i%4].set_xlim(-10, 10)
        axarr[i/4, i%4].set_ylim(-10, 10)
#        axarr[i/3, i%3].set_xlabel(r"R / Be star radius")
        axarr[i/4, i%4].set_ylabel(r"R / Be star radius") 
        axarr[i/4, i%4].set_title(r"M{}, {}e".format(Ms[i], es[i]))
        axarr[i/4, i%4].set(aspect='equal')
                     
#    cb = plt.colorbar(sc, ax=axarr[0, 0])
#    cb.ax.set_ylabel("Density in g per cm cubed")
                     

def findMinRelativeVelocity(baseNames, Ms, orbits):
    binaryPhase = np.array([])
    for k in range(len(baseNames)):
        closestVel = np.array([])
        bPhases = np.array([])
        filenames = np.array([])
#        filename = '/data/rob1g10/SPH/DATA/DISKS/'+baseNames[k]
        
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
    LaveSQRT = np.array([])
    aveMaxRhos = np.array([])
    minMaxRhos = np.array([])
    maxMaxRhos = np.array([])
    aveMinU = np.array([])
    minMinV = np.array([])
    aveDiscSize = np.array([])
    minDiscSize = np.array([])
    maxDiscSize = np.array([])
    
    # FOR WRITING DATA TO FILE

#    os.system("cp BHdataNEW.txt BHdataOLD.txt")
    
    fBH = open("BHdataOLD.txt", 'r')
    
    linesOLD = fBH.readlines()
        
    fBH.close()
    
    fBH = open("BHdataNEW.txt", 'w')
    fBH.close()
    
    baseNamesFULL=np.array([                                'M1.4_20d_0.0e', 
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
                                                            'M9_19d_0.1e', 
                                                            'M10_18.5d_0.1e',                     
                                                            
                                                            'M1.4_28d_0.2e', 
                                                            'M3_26.5d_0.2e', 
                                                            'M4_25.7d_0.2e', 
                                                            'M5_25d_0.2e',   
                                                            'M6_24.3d_0.2e', 
                                                            'M7_23.7d_0.2e', 
                                                            'M8_23.1d_0.2e',
                                                            'M9_22.6d_0.2e',  
                                                            'M10_22.1d_0.2e',
                                                            
                                                            "M1.4_34.1d_0.3e",
                                                            "M3_32.4d_0.3e",
                                                            "M4_31.4d_0.3e",
                                                            "M5_30.5d_0.3e",
                                                            "M6_29.7d_0.3e",
                                                            "M7_29d_0.3e",
                                                            "M8_28.3d_0.3e",
                                                            "M9_27.6d_0.3e",
                                                            "M10_27d_0.3e",   
                                                            
                                                            'M1.4_43d_0.4e', 
                                                            'M3_40.8d_0.4e', 
                                                            'M4_39.6d_0.4e', 
                                                            'M5_38.5d_0.4e', 
                                                            'M6_37.5d_0.4e', 
                                                            'M7_36.5d_0.4e', 
                                                            'M8_35.6d_0.4e',
                                                            'M9_34.8d_0.4e',  
                                                            'M10_34.1d_0.4e',
                                                            
                                                            'M1.4_56.6d_0.5e', 
                                                            'M3_53.7d_0.5e', 
                                                            'M4_52.1d_0.5e', 
                                                            'M5_50.6d_0.5e', 
                                                            'M6_49.2d_0.5e', 
                                                            'M7_48d_0.5e', 
                                                            'M8_46.8d_0.5e', 
                                                            'M9_45.8d_0.5e', 
                                                            'M10_44.8d_0.5e',
                                                            
                                                            'M1.4_79.1d_0.6e', 
                                                            'M3_75d_0.6e', 
                                                            'M4_72.8d_0.6e', 
                                                            'M5_70.7d_0.6e', 
                                                            'M6_68.8d_0.6e', 
                                                            'M7_67.1d_0.6e', 
                                                            'M8_65.5d_0.6e', 
                                                            'M9_64d_0.6e', 
                                                            'M10_62.6d_0.6e',
                                                            
                                                            'MWC656',])  
                                                            
    B = 0
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
        
        fBH = open("BHdataNEW.txt", 'a')
        
        os.system("ls /data/rob1g10/SPH/DATA/DISKS/{}/processedDISKs/QGDISK* > lsDISKSIZE.txt".format(baseNames[BN]))
        
#        try: 
        f = open("lsDISKSIZE.txt", 'r')
        lines = f.readlines()
    	
        DISKstart = int(lines[0][-6:-3])
        DISKend = int(lines[-1][-6:-3])
        DISKsteps = int(lines[-1][-3:])
    
        f.close()
        
#        print orbit[BN], baseNames[BN], DISKstart, DISKsteps
    
           
    
        DISKend, STEPend = calcQDISKsForOrbit(orbit[BN], baseNames[BN], 
                                              DISKstart, DISKend, DISKsteps, 
                                              numOrbits=2)
        
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
                maxRhos = np.append(maxRhos, np.max(rho)*umass / udist**3)
                minU = np.append(minU, np.min(u * (udist * 0.01)**2/utime**2))
                minV = np.append(minV, np.min( 
                                              np.abs( 
                                              
                                                     np.sqrt(  np.sum(partvels[:,2:]**2, 1)  ) 
                                                     
                                                   - np.sqrt(  np.sum(partvels[1]**2)      )
                                                     
                                                     
                                                    ) * udist * 0.01 / utime)
                                )
                discSize = np.append(discSize, findDiskSize([filenames[-1]]))
        
        days = times* utime / (24*3600)
        periastron = days[np.argmin(NSx)]
        binaryPhase = ( (periastron - days)  ) / orbit[BN]
        numOrbits = int((days[-1] - days[0]) / orbit[BN]) 
        orbitLen = int(len(days) / int(numOrbits))
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
        

        accTimes, accArray = createAccretionArray(disks,"NScapture", None, 'day', tol=0.01)
#        accArray *= 1e-15 * umass / (24*3600)
        accTimes = accTimes.astype(float)
        
#        print len(times)
        
        timeStart = np.argmin(  np.abs(accTimes-times[0]  )  )
        timeEnd = timeStart + np.argmin(  np.abs(accTimes[timeStart+2:]-times[-1] )  )
        
#        print times[0]* utime / (24*3600), times[-1]* utime / (24*3600)
#        plt.plot(range(len(accTimes)), accTimes)
        print timeStart, timeEnd
        
        accPhase = np.linspace(binaryPhase[0], binaryPhase[-1], len(accTimes[timeStart:timeEnd]) )
        accPhase = accPhase % 1
        
#        accTimes =  times * utime / (24*3600)
        accPhase = ( (accTimes * utime) / (24*3600) ) / orbit[BN]
#        accPhase -= accPhase[0]
#        print accTimes
#        print accPhase
#        print len(accArray)
#        plotAcc = accArray[timeStart:timeEnd]
        
        plotAcc = np.array([])
        plotPhase = np.array([])
        for T in range(len(times)):
            plotAcc = np.append(plotAcc, accArray[  np.argmin(  np.abs(accTimes-times[T])  )  ])
            plotPhase = np.append(plotPhase, accTimes[  np.argmin(  np.abs(accTimes-times[T])  )  ] * utime / (24*3600))
         
#        plotPhase = (times * utime / (24*3600))  
        plotPhase -= periastron
        plotPhase = plotPhase / orbit[BN]
        if(BN < 4):
            plotPhase += 0.5
        else:
            plotPhase -= 0.05
        accPhase -= accPhase[-50:][np.argmax(accArray[-50:])]
#        print len(accArray), accTimes[-10:-1].astype(float)* utime / (24*3600), days, disks

#        plotAccFOLD = plotAcc[:orbitLen]
#        for NUMORB in range(1, numOrbits):
#            plotAccFOLD += plotAcc[NUMORB*orbitLen:(NUMORB+1)*orbitLen] 
#        plotAccFOLD /= numOrbits
#        
#        
#        conversion = (1e-15 * 2.867e34 / (24*3600)) * (3e8**2 / 2)
##        plotAccFOLD = np.append(plotAccFOLD, np.append(plotAccFOLD, np.append(plotAccFOLD, np.append(plotAccFOLD, plotAccFOLD))))
##        plotPhaseFOLDstart = (plotPhase - int(np.min(plotPhase)))[0]
##        plotPhaseFOLD = np.linspace(plotPhaseFOLDstart, plotPhaseFOLDstart+5, len(plotAccFOLD))
#        
##        plotAcc = np.log10(plotAcc*conversion)
##        for PA in range(len(plotAcc)):
##            if(plotAcc[PA]<0):
##                plotAcc[PA] = 0
#        accArray = np.log10(accArray*conversion)
#        for AA in range(len(accArray)):
#            if(accArray[AA]<0):
#                accArray[AA] = 31
##        plt.plot(plotPhase - int(np.min(plotPhase)) - 1, plotAcc, color=colour[BN], label='M={}, e={}'.format(Mx[BN], baseNames[BN][-4:-1]))
#        plt.plot(accPhase-int(np.max(accPhase))+1.14, accArray, color=colour[BN], label=r'{}M$_\odot$, e={}'.format(Mx[BN], baseNames[BN][-4:-1]))
##        plt.plot(plotPhaseFOLD, plotAccFOLD, color=colour[BN], label='M={}, e={}'.format(Mx[BN], baseNames[BN][-4:-1]))
##        plt.plot(range(len(accTimes)), accTimes)
##        plt.plot(range(len(accPhase)), accPhase, color=colour[BN], label='{}'.format(Mx[BN]))
##        plt.scatter(range(len(binaryPhase)), binaryPhase, color=colour[BN], label='{}'.format(Mx[BN]))
        
        accSlice = int(5*orbit[BN])
        
        Lmax = np.append(Lmax, np.max(accArray[-accSlice:]))
        Lave = np.append(Lave, np.mean(accArray[-accSlice:]))
        LaveSQRT = np.append(LaveSQRT, np.mean(np.sqrt(accArray[-accSlice:])))
        aveMaxRhos = np.append(aveMaxRhos, np.mean(maxRhos))
        minMaxRhos = np.append(minMaxRhos, np.min(maxRhos))
        maxMaxRhos = np.append(maxMaxRhos, np.max(maxRhos))
        aveMinU = np.append(aveMinU, np.mean(minU) )
        minMinV = np.append(minMinV, np.min(minV))
#        print minMinV
        aveDiscSize = np.append(aveDiscSize, np.mean(discSize))
        minDiscSize = np.append(minDiscSize, np.min(discSize))
        maxDiscSize = np.append(maxDiscSize, np.max(discSize))
        
        for FULL in range(len(baseNamesFULL)):
            if(baseNamesFULL[B] == baseNames[BN]):
                fBH.write('{0:15} {1:15} {2:15} {3:15} {4:15} {5:15} {6:15} {7:15} {8:15} {9:15} {10:15} {11:15}\n'.format(baseNames[BN], Lmax[BN], Lave[BN], LaveSQRT[BN], aveMaxRhos[BN], minMaxRhos[BN], maxMaxRhos[BN], aveMinU[BN], minMinV[BN], aveDiscSize[BN], minDiscSize[BN], maxDiscSize[BN]))
                break
            else:
                fBH.write(linesOLD[B])
                B += 1
        B += 1
        print baseNames[BN], "\n"
        fBH.close()
#        except Exception as e:
#        fBH.write(linesOLD[B])
#        B += 1
#        print baseNames[BN], "failed. Skipping..."
#        print "Error is:"
#        print e, "\n"
#        fBH.close()
    
    fBH = open("BHdataNEW.txt", 'a')
    for FIN in range(B, len(baseNamesFULL)):
        fBH.write(linesOLD[FIN])
    fBH.close()
    
#    plt.xlim(0,1)
#    plt.ylim(31,34)
#    plt.ylabel(r"Average number of particles captured around an orbit")
#    plt.xlabel(r"binary phase")
#    plt.ylim(0, 350)
#    plt.legend(loc='upper right')
    
#    plt.twinx()
#    plt.ylabel(r"log L$_{X}$ (erg s$^{-1}$)")
#    plt.ylim(0, 350 * (1e-15 * umass / (24*3600)) * (3e8**2 / 2))
#    plt.legend(loc='best')
    
    return Lmax, Lave, LaveSQRT, aveMaxRhos, minMaxRhos, maxMaxRhos, aveMinU, minMinV, aveDiscSize, minDiscSize, maxDiscSize



def readBHdata(baseNames):
    f = open('BHdataOLD.txt', 'r')
    
    lines = f.readlines()
    
    f.close()
    
    Lmax = np.array([])
    Lave = np.array([])
    LaveSQRT = np.array([])
    aveRhos = np.array([])
    minRhos = np.array([])
    maxRhos = np.array([])
    aveMinU = np.array([])
    minV = np.array([])
    discSize = np.array([])
    minDiscSize = np.array([])
    maxDiscSize = np.array([])
    
    for line in range(len(lines)):
        if(np.any(lines[line].split()[0] == baseNames)):
            Lmax = np.append(Lmax, float(lines[line].split()[1]))
            Lave = np.append(Lave, float(lines[line].split()[2]))
            LaveSQRT = np.append(LaveSQRT, float(lines[line].split()[3]))
            aveRhos = np.append(aveRhos, float(lines[line].split()[4]))
            minRhos = np.append(minRhos, float(lines[line].split()[5]))
            maxRhos = np.append(maxRhos, float(lines[line].split()[6]))
            aveMinU = np.append(aveMinU, float(lines[line].split()[7]))
            minV = np.append(minV, float(lines[line].split()[8]))
            discSize = np.append(discSize, float(lines[line].split()[9]))
            minDiscSize = np.append(minDiscSize, float(lines[line].split()[10]))
            maxDiscSize = np.append(maxDiscSize, float(lines[line].split()[11]))
            
    return Lmax, Lave, LaveSQRT, aveRhos, minRhos, maxRhos, aveMinU, minV, discSize, minDiscSize, maxDiscSize



baseNames=np.array([

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
                    
#                    "M1.4_34.1d_0.3e",
#                    "M3_32.4d_0.3e",
#                    "M4_31.4d_0.3e",
#                    "M5_30.5d_0.3e",
#                    "M6_29.7d_0.3e",
#                    "M7_29d_0.3e",
#                    "M8_28.3d_0.3e",
#                    "M9_27.6d_0.3e",
#                    "M10_27d_0.3e",   
                    
#                    'M1.4_43d_0.4e', 
#                    'M3_40.8d_0.4e', 
#                    'M4_39.6d_0.4e', 
#                    'M5_38.5d_0.4e', 
#                    'M6_37.5d_0.4e', 
#                    'M7_36.5d_0.4e', 
#                    'M8_35.6d_0.4e',
#                    'M9_34.8d_0.4e',  
#                    'M10_34.1d_0.4e',
#                    
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
                    
                    'MWC656'
                    
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
for BN in baseNames[:-1]:
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
# MWC656
Mx = np.append(Mx, 5.3)
orbit = np.append(orbit, 60.37)

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
#gridPlotForDiscSize(['M1.4_20d_0.0e',
#                     'M1.4_28d_0.2e',
#                     'M1.4_43d_0.4e',
#                     'M1.4_79.1d_0.6e',
#                     
#                     'M6_17.4d_0.0e',
#                     'M6_24.3d_0.2e',
#                     'M6_37.5d_0.4e',
#                     'M6_68.8d_0.6e',
#                     
#                     'M10_15.8d_0.0e',
#                     'M10_22.1d_0.2e',
#                     'M10_34.1d_0.4e',
#                     'M10_62.6d_0.6e'], periastron='n')
#gridPlotForDiscSizeSIMPLE(['/data/rob1g10/SPH/DATA/DISKS/M1.4_20d_0.0e/processedDISKs/QGDISK19255',
#                         '/data/rob1g10/SPH/DATA/DISKS/M1.4_28d_0.2e/processedDISKs/QGDISK30934',
#                         '/data/rob1g10/SPH/DATA/DISKS/M1.4_43d_0.4e/processedDISKs/QGDISK53704',
#                         '/data/rob1g10/SPH/DATA/DISKS/M1.4_79.1d_0.6e/processedDISKs/QGDISK38701',
#                         
#                         '/data/rob1g10/SPH/DATA/DISKS/M6_17.4d_0.0e/processedDISKs/QGDISK39516',
#                         '/data/rob1g10/SPH/DATA/DISKS/M6_24.3d_0.2e/processedDISKs/QGDISK35951',
#                         '/data/rob1g10/SPH/DATA/DISKS/M6_37.5d_0.4e/processedDISKs/QGDISK15820',
#                         '/data/rob1g10/SPH/DATA/DISKS/M6_68.8d_0.6e/processedDISKs/QGDISK27505',
#                         
#                         '/data/rob1g10/SPH/DATA/DISKS/M10_15.8d_0.0e/processedDISKs/QGDISK32728',
#                         '/data/rob1g10/SPH/DATA/DISKS/M10_22.1d_0.2e/processedDISKs/QGDISK37220',
#                         '/data/rob1g10/SPH/DATA/DISKS/M10_34.1d_0.4e/processedDISKs/QGDISK14102',
#                         '/data/rob1g10/SPH/DATA/DISKS/M10_62.6d_0.6e/processedDISKs/QGDISK13509'])   
               
#Apastron
#gridPlotForDiscSizeSIMPLE(['/data/rob1g10/SPH/DATA/DISKS/M1.4_20d_0.0e/processedDISKs/QGDISK19315',
#                         '/data/rob1g10/SPH/DATA/DISKS/M1.4_28d_0.2e/processedDISKs/QGDISK30918',
#                         '/data/rob1g10/SPH/DATA/DISKS/M1.4_43d_0.4e/processedDISKs/QGDISK53903',
#                         '/data/rob1g10/SPH/DATA/DISKS/M1.4_79.1d_0.6e/processedDISKs/QGDISK38204',
#                         
#                         '/data/rob1g10/SPH/DATA/DISKS/M6_17.4d_0.0e/processedDISKs/QGDISK39533',
#                         '/data/rob1g10/SPH/DATA/DISKS/M6_24.3d_0.2e/processedDISKs/QGDISK36042',
#                         '/data/rob1g10/SPH/DATA/DISKS/M6_37.5d_0.4e/processedDISKs/QGDISK15836',
#                         '/data/rob1g10/SPH/DATA/DISKS/M6_68.8d_0.6e/processedDISKs/QGDISK27501',
#                         
#                         '/data/rob1g10/SPH/DATA/DISKS/M10_15.8d_0.0e/processedDISKs/QGDISK32745',
#                         '/data/rob1g10/SPH/DATA/DISKS/M10_22.1d_0.2e/processedDISKs/QGDISK37237',
#                         '/data/rob1g10/SPH/DATA/DISKS/M10_34.1d_0.4e/processedDISKs/QGDISK14118',
#                         '/data/rob1g10/SPH/DATA/DISKS/M10_62.6d_0.6e/processedDISKs/QGDISK13413'])
    
#print findMinRelativeVelocity(baseNames[16:24], Mx[16:24], orbit[16:24])
#for i in range(16,24):
#    phases = findBinaryPhase(baseNames[i], orbit[i] , colour = ['b', 'g', 'r', 'y', 'black', 'brown', 'turquoise', 'pink'][i-16])
#    print phases
    
Lmax, Lave, LaveSQRT, aveRhos, minRhos, maxRhos, aveMinU, minV, discSize, minDiscSize, maxDiscSize = plotCOaccretion(baseNames, Mx, orbit, colour = ['b', 'g', 'r', 'y', 'black', 'brown', 'turquoise', 'pink'])

#Lmax, Lave, LaveSQRT, aveRhos, minRhos, maxRhos, aveMinU, minV, discSize, minDiscSize, maxDiscSize = readBHdata(baseNames)



#LaveSQRT2 = np.delete(LaveSQRT, -1)

#Rx = 2 * Mx * 6.67e-11 / 3e8**2

sec1 = 0
sec2 = sec1 + 9
sec3 = sec2 + 9
sec4 = sec3 + 9
sec5 = sec4 + 8
sec6 = sec5 + 9
sec7 = sec6 + 7
sec8 = sec7 + 9
MWC656 = sec8 + 1

conversion = (1e-15 * 2.867e34 / (24*3600)) * (3e8**2 / 2)

#def plotLave(Mx, Lave):
#    return
#    
#def plotLmax():
#    return
#    
#def plotRhos():
#    return
#    
#def plotDiscSize():
#    return

#e06=plt.errorbar(Mx[sec7:sec8], np.log10(Lave[sec7:sec8]*conversion), 
#             marker='v', ms=8, color='turquoise', label="e = 0.6")
#plt.errorbar(Mx[sec7:sec8], np.log10(Lave[sec7:sec8]*conversion), 
#             yerr= np.log10(Lave[sec7:sec8]*conversion) - np.log10(LaveSQRT[sec7:sec8]*conversion),
#             marker='v', ms=8, color='turquoise')
#e05=plt.errorbar(Mx[sec6:sec7], np.log10(Lave[sec6:sec7]*conversion), 
#             marker='>', ms=8, color='orange', label="e = 0.5")
#plt.errorbar(Mx[sec6:sec7], np.log10(Lave[sec6:sec7]*conversion), 
#             yerr= np.log10(Lave[sec6:sec7]*conversion) - np.log10(LaveSQRT[sec6:sec7]*conversion),
#             marker='>', ms=8, color='orange')
#e04=plt.errorbar(Mx[sec5:sec6], np.log10(Lave[sec5:sec6]*conversion), 
#             marker='<', ms=8, color='purple', label="e = 0.4")
#plt.errorbar(Mx[sec5:sec6], np.log10(Lave[sec5:sec6]*conversion), 
#             yerr= np.log10(Lave[sec5:sec6]*conversion) - np.log10(LaveSQRT[sec5:sec6]*conversion),
#             marker='<', ms=8, color='purple')
#e03=plt.errorbar(Mx[sec4:sec5], np.log10(Lave[sec4:sec5]*conversion),
#             marker='s', ms=8, color='black', label="e = 0.3")
#plt.errorbar(Mx[sec4:sec5], np.log10(Lave[sec4:sec5]*conversion),
#             yerr= np.log10(Lave[sec4:sec5]*conversion) - np.log10(LaveSQRT[sec4:sec5]*conversion),
#             marker='s', ms=8, color='black')
#e02=plt.errorbar(Mx[sec3:sec4], np.log10(Lave[sec3:sec4]*conversion), 
#             marker='^', ms=8, color='g', label="e = 0.2")
#plt.errorbar(Mx[sec3:sec4], np.log10(Lave[sec3:sec4]*conversion), 
#             yerr= np.log10(Lave[sec3:sec4]*conversion) - np.log10(LaveSQRT[sec3:sec4]*conversion),
#             marker='^', ms=8, color='g')
#e01=plt.errorbar(Mx[sec2:sec3], np.log10(Lave[sec2:sec3]*conversion), 
#             marker='x', ms=8, color='r', label="e = 0.1")
#plt.errorbar(Mx[sec2:sec3], np.log10(Lave[sec2:sec3]*conversion), 
#             yerr= np.log10(Lave[sec2:sec3]*conversion) - np.log10(LaveSQRT[sec2:sec3]*conversion),
#             marker='x', ms=8, color='r')
##plt.scatter(Mx[sec2:sec3], Lmax[sec2:sec3]), marker='x', color='r', label="e = 0.1")
#e00=plt.errorbar(Mx[:sec2], np.log10(Lave[:sec2]*conversion), 
#             marker='o', ms=8, color='b', label="e = 0.0")
#plt.errorbar(Mx[:sec2], np.log10(Lave[:sec2]*conversion), 
#             yerr= np.log10(Lave[:sec2]*conversion) - np.log10(LaveSQRT[:sec2]*conversion),
#             marker='o', ms=8, color='b')
##plt.ylabel(r"Peak number of particles captured around an orbit")
#plt.xlabel(r"compact object mass (M$_{\odot}$)")
#plt.xlim(1, 14)
#plt.legend(handler_map={e00: HandlerErrorbar(numpoints=1), 
#                        e01: HandlerErrorbar(numpoints=1), 
#                        e02: HandlerErrorbar(numpoints=1), 
#                        e03: HandlerErrorbar(numpoints=1), 
#                        e04: HandlerErrorbar(numpoints=1),
#                        e05: HandlerErrorbar(numpoints=1), 
#                        e06: HandlerErrorbar(numpoints=1), },loc="upper right")
#plt.ylabel(r"log L$_{ave}$ (erg s$^{-1}$)")


#e06=plt.errorbar(Mx[sec7:sec8], np.log10(Lmax[sec7:sec8]*conversion), 
#             marker='v', ms=8, color='turquoise', label="e = 0.6")
#plt.errorbar(Mx[sec7:sec8], np.log10(Lmax[sec7:sec8]*conversion), 
#             yerr= np.log10(Lmax[sec7:sec8]*conversion) - np.log10(np.sqrt(Lmax[sec7:sec8])*conversion),
#             marker='v', ms=8, color='turquoise')
#e05=plt.errorbar(Mx[sec6:sec7], np.log10(Lmax[sec6:sec7]*conversion), 
#             marker='>', ms=8, color='orange', label="e = 0.5")
#plt.errorbar(Mx[sec6:sec7], np.log10(Lmax[sec6:sec7]*conversion), 
#             yerr= np.log10(Lmax[sec6:sec7]*conversion) - np.log10(np.sqrt(Lmax[sec6:sec7])*conversion),
#             marker='>', ms=8, color='orange')
#e04=plt.errorbar(Mx[sec5:sec6], np.log10(Lmax[sec5:sec6]*conversion), 
#             marker='<', ms=8, color='purple', label="e = 0.4")
#plt.errorbar(Mx[sec5:sec6], np.log10(Lmax[sec5:sec6]*conversion), 
#             yerr= np.log10(Lmax[sec5:sec6]*conversion) - np.log10(np.sqrt(Lmax[sec5:sec6])*conversion),
#             marker='<', ms=8, color='purple')
#e03=plt.errorbar(Mx[sec4:sec5], np.log10(Lmax[sec4:sec5]*conversion),
#             marker='s', ms=8, color='black', label="e = 0.3")
#plt.errorbar(Mx[sec4:sec5], np.log10(Lmax[sec4:sec5]*conversion),
#             yerr= np.log10(Lmax[sec4:sec5]*conversion) - np.log10(np.sqrt(Lmax[sec4:sec5])*conversion),
#             marker='s', ms=8, color='black')
#e02=plt.errorbar(Mx[sec3:sec4], np.log10(Lmax[sec3:sec4]*conversion), 
#             marker='^', ms=8, color='g', label="e = 0.2")
#plt.errorbar(Mx[sec3:sec4], np.log10(Lmax[sec3:sec4]*conversion), 
#             yerr= np.log10(Lmax[sec3:sec4]*conversion) - np.log10(np.sqrt(Lmax[sec3:sec4])*conversion),
#             marker='^', ms=8, color='g')
#e01=plt.errorbar(Mx[sec2:sec3], np.log10(Lmax[sec2:sec3]*conversion), 
#             marker='x', ms=8, color='r', label="e = 0.1")
#plt.errorbar(Mx[sec2:sec3], np.log10(Lmax[sec2:sec3]*conversion), 
#             yerr= np.log10(Lmax[sec2:sec3]*conversion) - np.log10(np.sqrt(Lmax[sec2:sec3])*conversion),
#             marker='x', ms=8, color='r')
##plt.scatter(Mx[sec2:sec3], Lmax[sec2:sec3]), marker='x', color='r', label="e = 0.1")
#e00=plt.errorbar(Mx[:sec2], np.log10(Lmax[:sec2]*conversion), 
#             marker='o', ms=8, color='b', label="e = 0.0")
#plt.errorbar(Mx[:sec2], np.log10(Lmax[:sec2]*conversion), 
#             yerr= np.log10(Lmax[:sec2]*conversion) - np.log10(np.sqrt(Lmax[:sec2])*conversion),
#             marker='o', ms=8, color='b')
##plt.ylabel(r"Peak number of particles captured around an orbit")
#plt.xlabel(r"compact object mass (M$_{\odot}$)")
#plt.xlim(1, 14)
#plt.legend(handler_map={e00: HandlerErrorbar(numpoints=1), 
#                        e01: HandlerErrorbar(numpoints=1), 
#                        e02: HandlerErrorbar(numpoints=1), 
#                        e03: HandlerErrorbar(numpoints=1), 
#                        e04: HandlerErrorbar(numpoints=1),
#                        e05: HandlerErrorbar(numpoints=1), 
#                        e06: HandlerErrorbar(numpoints=1), },loc="upper right")
#plt.ylabel(r"log L$_{max}$ (erg s$^{-1}$)")
##plt.twinx()
#plt.ylabel(r"log L$_{max}$ (erg s$^{-1}$)")
##plt.ylim(0, 350 * (1e-15 * umass / (24*3600)) * (3e8**2 / 2))

#tempRho = np.array(aveRhos)
#aveRhos = np.log(aveRhos)
#tempminRho = np.array(minRhos)
#minRhos = np.log(minRhos)
#tempmaxRho = np.array(maxRhos)
#maxRhos = np.log(maxRhos)


#e06=plt.errorbar(Mx[sec7:sec8], np.log10(aveRhos[sec7:sec8]), marker='v', color='turquoise', label="e = 0.6")
#plt.errorbar(Mx[sec7:sec8], np.log10(aveRhos[sec7:sec8]), yerr=[np.log10(aveRhos[sec7:sec8])-np.log10(minRhos[sec7:sec8]), np.log10(maxRhos[sec7:sec8])-np.log10(aveRhos[sec7:sec8])], color='turquoise', marker='v', ms=8)
#
#e05=plt.errorbar(Mx[sec6:sec7], np.log10(aveRhos[sec6:sec7]), marker='>', color='orange', label="e = 0.5")
#plt.errorbar(Mx[sec6:sec7], np.log10(aveRhos[sec6:sec7]), yerr=[np.log10(aveRhos[sec6:sec7])-np.log10(minRhos[sec6:sec7]), np.log10(maxRhos[sec6:sec7])-np.log10(aveRhos[sec6:sec7])], color='orange', marker='>', ms=8)
#
#p6 = np.polyfit(Mx[sec5+5:], aveRhos[sec5+5:], 1)
#print p6[0], p6[1]
#p6 = np.polyfit(Mx[sec5+2:sec5+6], aveRhos[sec5+2:sec5+6], 1)
#print p6[0], p6[1]
#e04=plt.errorbar(Mx[sec5:sec6], np.log10(aveRhos[sec5:sec6]), marker='<', color='purple', label="e = 0.4")
#plt.errorbar(Mx[sec5:sec6], np.log10(aveRhos[sec5:sec6]), yerr=[np.log10(aveRhos[sec5:sec6])-np.log10(minRhos[sec5:sec6]), np.log10(maxRhos[sec5:sec6])-np.log10(aveRhos[sec5:sec6])], color='purple', marker='<', ms=8)
#
#p5 = np.polyfit(Mx[sec4+4:sec5], np.log10(aveRhos[sec4+4:sec5]), 1)
#print p5[0], p5[1]
#plt.plot(Mx[sec4+4:sec5], Mx[sec4+4:sec5]*p5[0] + p5[1], color='black')
#
#p4 = np.polyfit(Mx[sec4+1:sec4+4], np.log10(aveRhos[sec4+1:sec4+4]), 1)
#print p4[0], p4[1]
#
##plt.plot(Mx[sec4+1:sec4+5], Mx[sec4+1:sec4+5]*p4[0] + p4[1], color='black')
#e03=plt.errorbar(Mx[sec4:sec5], np.log10(aveRhos[sec4:sec5]), marker='s', color='black', label="e = 0.3")
#plt.errorbar(Mx[sec4:sec5], np.log10(aveRhos[sec4:sec5]), yerr=[np.log10(aveRhos[sec4:sec5])-np.log10(minRhos[sec4:sec5]), np.log10(maxRhos[sec4:sec5])-np.log10(aveRhos[sec4:sec5])], marker='s', color='black', ms=8)
#p3 = np.polyfit(Mx[sec3:sec4], np.log10(aveRhos[sec3:sec4]), 1)
#print p3[0], p3[1]
#plt.plot(Mx[sec3:sec4], Mx[sec3:sec4]*p3[0] + p3[1], color='g')
#e02=plt.errorbar(Mx[sec3:sec4], np.log10(aveRhos[sec3:sec4]), marker='^', color='g', label="e = 0.2")
#plt.errorbar(Mx[sec3:sec4], np.log10(aveRhos[sec3:sec4]), yerr=[np.log10(aveRhos[sec3:sec4])-np.log10(minRhos[sec3:sec4]), np.log10(maxRhos[sec3:sec4])-np.log10(aveRhos[sec3:sec4])], marker='^', color='g', ms=8)
#
#p2 = np.polyfit(Mx[sec2:sec3], np.log10(aveRhos[sec2:sec3]), 1)
#print p2[0], p2[1]
#plt.plot(Mx[sec2:sec3], Mx[sec2:sec3]*p2[0] + p2[1], color='r')
#e01=plt.errorbar(Mx[sec2:sec3], np.log10(aveRhos[sec2:sec3]), marker='x', color='r', label="e = 0.1")
#plt.errorbar(Mx[sec2:sec3], np.log10(aveRhos[sec2:sec3]), yerr=[np.log10(aveRhos[sec2:sec3])-np.log10(minRhos[sec2:sec3]), np.log10(maxRhos[sec2:sec3])-np.log10(aveRhos[sec2:sec3])], marker='x', color='r', ms=8)
#
#p1 = np.polyfit(Mx[:sec2], np.log10(aveRhos[:sec2]), 1)
#print p1[0], p1[1]
#plt.plot(Mx[:sec2], Mx[:sec2]*p1[0] + p1[1], color='b')
#e00=plt.errorbar(Mx[:sec2], np.log10(aveRhos[:sec2]), marker='o', color='b', label="e = 0.0")
#plt.errorbar(Mx[:sec2], np.log10(aveRhos[:sec2]), yerr=[np.log10(aveRhos[:sec2])-np.log10(minRhos[:sec2]), np.log10(maxRhos[:sec2])-np.log10(aveRhos[:sec2])], color='b',  marker='o', ms=8)
#
#plt.xlim(1, 14)
#plt.ylabel(r"log base gas density (g cm$^{-3}$)")
#plt.xlabel(r"compact object mass (M$_{\odot}$)")
#plt.legend(handler_map={e00: HandlerErrorbar(numpoints=1), 
#                        e01: HandlerErrorbar(numpoints=1), 
#                        e02: HandlerErrorbar(numpoints=1), 
#                        e03: HandlerErrorbar(numpoints=1), 
#                        e04: HandlerErrorbar(numpoints=1),
#                        e05: HandlerErrorbar(numpoints=1), 
#                        e06: HandlerErrorbar(numpoints=1), },loc="upper right")

#aveRhos = tempRho
#minRhos = tempminRho
#maxRhos = tempmaxRho

##c_s = np.sqrt(2 * minU / 3)
#BR = 2*6.67e-11*Mx*2e30 / minV**2
#plt.plot(Mx[sec5:], BR[sec5:] / (7*7e10), color='purple', label="e = 0.6")
#plt.plot(Mx[sec4:sec5], BR[sec4:sec5] / (7*7e10), color='black', label="e = 0.4")
#plt.plot(Mx[sec3:sec4], BR[sec3:sec4] / (7*7e10), color='g', label="e = 0.2")
##plt.plot(Mx[sec2:sec3], BR[sec2:sec3] / (7*7e10), color='r', label="e = 0.1")
#plt.plot(Mx[:sec2], np.log10(BR[:sec2] / (7*7e10)), color='b', label="e = 0.0")
#plt.ylabel("Bondi accretion radius (stellar radii)")
#plt.xlabel("Compact object mass")
#
#q = 13. / Mx
#rocheLobe = ( 0.49*q**(2./3.) ) / ( 0.6*q**(2./3.) + np.log(1. + q**(1./3.)) )
##rocheLobeSimple = 0.38 + 0.2*np.log10(q) 
##plt.scatter(Mx[sec4:], np.log10(rocheLobe[sec4:]), marker='s', color='black', label="e = 0.6")
##plt.scatter(Mx[sec3:sec4], np.log10(rocheLobe[sec3:sec4]), marker='^', color='g', label="e = 0.4")
##plt.scatter(Mx[sec2:sec3], np.log10(rocheLobe[sec2:sec3]), marker='x', color='r', label="e = 0.2")
#plt.plot(Mx[:sec2], rocheLobe[:sec2]*10.9, color='brown', label='Roche Lobe')
##plt.scatter(Mx[:sec2], rocheLobeSimple[:sec2], marker='o', color='b')
#plt.ylabel("Roche lobe radius (stellar radii)")
#plt.xlabel("Compact object mass")

#plt.errorbar(Mx[sec7:sec8], discSize[sec7:sec8], yerr=[discSize[sec7:sec8]-minDiscSize[sec7:sec8], maxDiscSize[sec7:sec8]-discSize[sec7:sec8]], color='turquoise')
#plt.plot(Mx[sec7:sec8], discSize[sec7:sec8], color='turquoise', label="e = 0.6")
#
#plt.errorbar(Mx[sec6:sec7], discSize[sec6:sec7], yerr=[discSize[sec6:sec7]-minDiscSize[sec6:sec7], maxDiscSize[sec6:sec7]-discSize[sec6:sec7]], color='orange')
#plt.plot(Mx[sec6:sec7], discSize[sec6:sec7], color='orange', label="e = 0.5")
#
##p5 = np.polyfit(Mx[sec5:sec6], discSize[sec5:sec6], 1)
##print p5[0], p5[1]
##plt.scatter(Mx[sec5:], discSize[sec5:], marker='<', color='purple', label="e = 0.6")
#plt.errorbar(Mx[sec5:sec6], discSize[sec5:sec6], yerr=[discSize[sec5:sec6]-minDiscSize[sec5:sec6], maxDiscSize[sec5:sec6]-discSize[sec5:sec6]], color='purple')
#plt.plot(Mx[sec5:sec6], discSize[sec5:sec6], color='purple', label="e = 0.4")
#
##p4 = np.polyfit(Mx[sec4:sec5], discSize[sec4:sec5], 1)
##print p4[0], p4[1]
##plt.scatter(Mx[sec4:sec5], discSize[sec4:sec5], marker='s', color='black', label="e = 0.4")
#plt.errorbar(Mx[sec4:sec5], discSize[sec4:sec5], yerr=[discSize[sec4:sec5]-minDiscSize[sec4:sec5], maxDiscSize[sec4:sec5]-discSize[sec4:sec5]], color='black')
#plt.plot(Mx[sec4:sec5], discSize[sec4:sec5], color='black', label="e = 0.3")
#
##p3 = np.polyfit(Mx[sec3:sec4], discSize[sec3:sec4], 1)
##print p3[0], p3[1]
##plt.scatter(Mx[sec3:sec4], discSize[sec3:sec4], marker='^', color='g', label="e = 0.2")
#plt.errorbar(Mx[sec3:sec4], discSize[sec3:sec4], yerr=[discSize[sec3:sec4]-minDiscSize[sec3:sec4], maxDiscSize[sec3:sec4]-discSize[sec3:sec4]], color='g')
#plt.plot(Mx[sec3:sec4], discSize[sec3:sec4], color='g', label="e = 0.2")
#
##p2 = np.polyfit(np.delete(Mx[sec2:sec3], 2), np.delete(discSize[sec2:sec3], 2), 1)
##print p2[0], p2[1]
##plt.scatter(Mx[sec2:sec3], discSize[sec2:sec3], marker='x', color='r', label="e = 0.1")
#plt.errorbar(Mx[sec2:sec3], discSize[sec2:sec3], yerr=[discSize[sec2:sec3]-minDiscSize[sec2:sec3], maxDiscSize[sec2:sec3]-discSize[sec2:sec3]], color='r')
#plt.plot(Mx[sec2:sec3], discSize[sec2:sec3], color='r', label="e = 0.1")
#
##p1 = np.polyfit(Mx[:sec2], discSize[:sec2], 1)
##print p1[0], p1[1]
##plt.scatter(Mx[:sec2], discSize[:sec2], marker='o', color='b', label="e = 0.0")
#plt.errorbar(Mx[:sec2], discSize[:sec2], yerr=[discSize[:sec2]-minDiscSize[:sec2], maxDiscSize[:sec2]-discSize[:sec2]], color='b')
#plt.plot(Mx[:sec2], discSize[:sec2], color='b', label="e = 0.0")
#
#plt.ylabel(r"disc size (R$_{\odot}$)")
#plt.xlabel(r"compact object mass (M$_{\odot}$)")
#
#plt.legend(loc='upper right')


#------------------------------------------------------------------------------
#         Histogram of real data              
#------------------------------------------------------------------------------

eccGalacticCoe = np.array([0.03, 0.16, 0.28, 0.3, 0.34, 0.42, 0.54])
eccLMCCoe = np.array([])
eccSMCCoe = np.array([0.07, 0.16, 0.26, 0.28, 0.4, 0.41, 0.43])

eccSMCTownsend = np.array([0.07, 0.16, 0.26, 0.41, 0.28, 0.43, 0.4])
LmaxSMCTownsend = np.array([2.1e38, 5e37, 3.3e37, 7.3e37, 1e37, 3e37, 3.5e37])

eccGalactic = np.array([0.03, 0.03, 0.036, 0.09, 0.12, 0.26, 0.3, 0.33, 0.34, 0.4, 0.41, 0.446, 0.47, 0.537, 0.66, 0.88,
                        0.16, 0.42, 0.34, 0.28, 0.03, 0.3
                        ])
LmaxGalactic = np.array([1.1e37, 2.1e37, 1.1e38, 7e36, 1.1e37, 3.9e34, 1.1e38, 5.4e36, 3e37, 3e36, 1e38, 8e36, 2e37, 2e34, 2.9e35, 6e36,
                         4.9e35, 1.4e38, 1.4e36, 2.3e37, 8.8e37, 8.8e34
                         ])
eccLMCRaguzova = np.array([])
LmaxLMCRaguzova = np.array([])
eccSMCRaguzova = np.array([])



#eccGalactic = np.append(eccGalacticCoe, eccGalactic)
#eccLMC = np.append(eccLMCCoe, eccLMCRaguzova)
#eccSMC = np.append(eccSMCCoe, eccSMCRaguzova)

#plt.hist(np.append(eccGalactic, np.append(eccLMC, eccSMC)), bins=np.linspace(0, 1.0, 20), label="SMC")
##plt.hist(np.append(eccGalactic, eccLMC), bins=np.linspace(0, 1.0, 20), label="LMC")
#plt.hist(eccGalactic, bins=np.linspace(0, 1.0, 20), label="Galactic")
#plt.ylabel("Number of Be/neutron star systems")
#plt.xlabel("Eccentricity")
#plt.legend()
#plt.rcParams.update({'font.size': 14})


#SMC = plt.errorbar(eccSMCTownsend, np.log10(LmaxSMCTownsend), label="SMC", color='b', marker='o', ms=10, linestyle='None')
##plt.scatter(eccLMCRaguzova, np.log10(LmaxLMCRaguzova), label="LMC", color='g', marker='^')
#Galactic = plt.errorbar(eccGalactic, np.log10(LmaxGalactic), label="Galactic", color='g', marker='<', ms=10, linestyle='None')
#MWC656 = plt.errorbar(0.1, np.log10(3.7e31), color='r', marker='s', ms=10, linestyle='None')
#plt.text(0.12, np.log10(4.3e31), "MWC 656", fontsize=15)
#plt.xlim(0, 1)
#plt.ylabel(r"log L$_{max}$ (erg s$^{-1}$)")
#plt.xlabel("Eccentricity")
#plt.legend(handler_map={MWC656: HandlerErrorbar(numpoints=1), 
#                        SMC: HandlerErrorbar(numpoints=1), 
#                        Galactic: HandlerErrorbar(numpoints=1)},
#                        loc="lower right")


#------------------------------------------------------------------------------
#         PSR              
#------------------------------------------------------------------------------

# X-ray dip is from aprrox Oct 13th to approx Nov 22nd with periastron at ~13 Nov

# Therefore 31 days before periastron and 9 days after