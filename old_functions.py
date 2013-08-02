
def doCCDMapOld(SEDMap, u, minLambda, maxLambda, deltaLambda, minOrder, maxOrder, deltaOrder, fLength, stheta, intNormalize,Interpolate=False,BackImage='',GaussFit=False):

    
    dataOut=np.zeros(5)

    
    #Loads SEDMap based on selection. 
#    SEDMap = wt.do_SED_map(SEDMode, minLambda, maxLambda, deltaLambda, intNormalize)
    blaze_angle= stheta #Approximately np.arctan(2)
    allFlux=np.array([0])
    allLambdas=np.array([0])
    
    '''Main loop
    Navigates orders within the range given
    For each order navigates the list of wavelenghts in SEDMap constrained by +/- FSP/2'''    
    for nOrder in range(minOrder, maxOrder + deltaOrder, deltaOrder):

        #the wavelength range is from the blaze wavelength of the next order and the blaze wavelength of the previous order
        LambdaBlMin = 2*d*np.sin(blaze_angle)/(nOrder+1) #Was 0.5
        LambdaBlMax = 2*d*np.sin(blaze_angle)/(nOrder-1) #Was 0.5
 
#        #SEDMapLoop is an array of wavelengths (paired with intensity) called SEDMap that are in this particular order.
#        SEDMapLoop=SEDMap.copy()
#        #constrain by +/- FSP (was FSP/2)
#        SEDMapLoop = SEDMapLoop[SEDMapLoop[:,0]>=LambdaBlMin]
#        if SEDMapLoop.shape[0]>0:       
#            SEDMapLoop = SEDMapLoop[SEDMapLoop[:,0]<=LambdaBlMax]     

        #loop lambda for current order
        for Lambda,inI in SEDMap[SEDMap[:,0]>=LambdaBlMin][SEDMap[:,0]<=LambdaBlMax]: 
            
            #refractive indeces of prism and air
            nPrism = wt.n(Lambda,1)
            nAir = wt.n(Lambda,0)
            
            #Computes the unit vector that results from the optical system for a given wavelength and order
            #This is the actual tracing of the ray for each wavelength
#            v, isValid = wt.rayTrace(nAir, nPrism, nOrder, Lambda, d, u, n1, n2, n4, n5, s, l)
            v, isValid = wt.ray_trace_flex(Beam, Lambda, nOrder, Optics)
            
            if isValid: #no errors in calculation and ray makes it back through the prism
                x=v[0]*fLength*1000/pixelSize # x-coord in focal plane in pixels
                z=v[2]*fLength*1000/pixelSize # z-coord in focal plane in pixels
    
                '''Appends data table
                x,z=pojected coordinates
                inI, outI = intensities, inI=from source (file, random,...) 0.0 to 1.0'''
                outI=wt.Intensity(Lambda, minLambda, maxLambda)              
                dataOut= np.vstack((dataOut,np.array([x,z, Lambda, inI*outI ,nOrder]))) 
                
        
        #Order extraction
        if (Interpolate==True and len(dataOut[dataOut[:,4]==nOrder][:,0])>=3):
            xPlot=dataOut[dataOut[:,4]==nOrder][:,0]
            yPlot=dataOut[dataOut[:,4]==nOrder][:,1]
            LambdaPlot=dataOut[dataOut[:,4]==nOrder][:,2]
            
            fLambda = interpolate.interp1d(yPlot, LambdaPlot)
            fX = interpolate.interp1d(yPlot, xPlot, 'quadratic', bounds_error=False)
          
            
            hdulist = pyfits.open(BackImage)
            imWidth = hdulist[0].header['NAXIS1']
            imHeight = hdulist[0].header['NAXIS2']
            
            newY=np.arange(-imHeight/2,imHeight/2)
            
            newX=fX(newY)
            
            nanMap= np.isnan(newX)
            newX=newX[-nanMap]
            newY=newY[-nanMap]
                    
            flux,flux2,flux3 = wt.extractOrder(newX,newY,BackImage)
            
            #read flats
            image ='simple_flat.fits'
            fluxFlat,flux2,flux3 = wt.extractOrder(newX,newY,image)
            
            
            Lambdas=fLambda(newY)
            
            #Blackbody curve to balance flats
            BB=Lambdas**(-4) / (np.exp(14400/Lambdas/3000)- 1)
         
         
            cleanFlux= flux/fluxFlat*BB#/nOrder**2#/np.max(fluxFlat))

            #Write flux to files
#            f = open('Order_'+ str(nOrder) +'.txt', 'w')
#            for k in range(0,len(Lambdas)):
#                outText=str(Lambdas[k])+'\t'+ str(cleanFlux[k])+'\n'
#                f.write(outText)
#            f.close()
#            
#            peakIndex=np.where(normalizedFlux==np.max(normalizedFlux))[0][0]
#            peakIndex=len(normalizedFlux)/2
#            Range=(len(normalizedFlux)- peakIndex)*0.6
            

            # Fit a Gaussian             
            if GaussFit==True: 
#                # Create some sample data
#                known_param = np.array([2.0, .7])
#                xmin,xmax = -1.0, 5.0
#                N = 1000
#                X = np.linspace(xmin,xmax,N)
#                Y = gauss(X, known_param)
#                
#                # Add some noise
#                Y += .10*np.random.random(N)
#                
#                # Renormalize to a proper PDF
#                Y /= ((xmax-xmin)/len(Y))*Y.sum()  
                
                
                         
                X=Lambdas.copy()
                Y=cleanFlux.copy()
                if len(Y)>0:
                    #Y /= ((np.max(X) - np.min(X))/len(Y))*Y.sum()  
                    
                    a,FWHMIndex = wt.find_nearest(Y,np.max(Y)/2)
                    maxIndex=np.where(Y==np.max(Y))[0][0]
                    FWHM=2*(X[maxIndex]-X[FWHMIndex])
                    fit_mu=X[maxIndex]
                    #print 'FWHM',FWHM
#                cleanX=X[np.arange(0,len(X),1)].copy()
#                cleanY=Y[np.arange(0,len(Y),1)].copy()

#                    p0 = [0,1] # Inital guess is a normal distribution
#                    errfunc = lambda p, x, y: gauss(x, p) - y # Distance to the target function
#                    p1, success = leastsq(errfunc, p0[:], args=(X,Y))
#                    
#                    fit_mu, fit_stdev = p1
#                    
#                    p1[1]=5*fit_stdev
#                    
#                    FWHM = 4*np.sqrt(2*np.log(2))*fit_stdev
#                    print 'FWHM2',FWHM
                    
                    R=fit_mu/FWHM
                    print fit_mu ,'            ', R
                    
                    plt.plot(X,Y) 
    #                plt.plot(cleanX,cleanY) 
                    plt.ylabel('Intensity (Relative Units)')
                    plt.xlabel('Wavelength (Micrometers)')
                    plt.title('Spectral Resolution at '+ str("{:0.4f}".format(fit_mu))+' micrometers')
                    plt.annotate('FWHM='+str("{:0.4f}".format(FWHM*1000))+'nm R=' + str("{:0.4f}".format(fit_mu/FWHM)), 
                        xy = (X[0],Y[0]), xytext = (220, 250),
                        textcoords = 'offset points', ha = 'right', va = 'bottom',
                        bbox = dict(boxstyle = 'round,pad=0.5', fc = 'white', alpha = 0.9), size=15)
                    #plt.plot(X, gauss(X,p1),lw=3,alpha=.5, color='r')
                    plt.axvspan(fit_mu-FWHM/2, fit_mu+FWHM/2, facecolor='g', alpha=0.5)
                    plt.show()
#            
#            print nOrder
#            print str("{:0.4f}".format(p1[0]/FWHM))
#            print str("{:0.4f}".format(p1[0]))
            
            if np.sum(allFlux)>0:
                intersectStart=bis.bisect(allLambdas,np.min(Lambdas))   
                intersectEnd=len(allLambdas)   
                bestDistance=1e10     
                bestIndex=0
                
                for k in range(0,intersectEnd-intersectStart):
                    currDistance=np.sqrt((allFlux[intersectStart+k]-cleanFlux[k])**2)
                    if currDistance<bestDistance:
                        bestDistance=currDistance
                        bestIndex=k
                        
                allLambdas=allLambdas[allLambdas<Lambdas[bestIndex]]
                allFlux=allFlux[allLambdas<Lambdas[bestIndex]]
                

                
                allLambdas=np.hstack((allLambdas,Lambdas[bestIndex:]))
                allFlux=np.hstack((allFlux,cleanFlux[bestIndex:]))
            else:
                allLambdas=Lambdas
                allFlux=cleanFlux
             
#            if Interpolate==True:   
#                
#                fig = plt.figure()
#                plt.ylabel('Intensity (Counts)')
#                plt.xlabel('Wavelength (Micrometers)')
#                plt.title('Order '+ str(nOrder))
#                ax1 = fig.add_subplot(111)          
#                ax1.plot(Lambdas,cleanFlux)
##                ax1.plot(Lambdas,fluxFlat)
#                ax1.plot(Lambdas,flux3)
#                #plt.savefig(str(nOrder)+'.png')
#                plt.show()
#    

    if Interpolate==True:   
#        a=np.sort(allFlux,0)
#        params = {'legend.fontsize': 200,
#          'legend.linewidth': 2}
#        plt.rcParams.update(params)
        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        ax1.plot(allLambdas,allFlux)
#        ax1.scatter(allLambdas,allFlux ,s=0.1, color='black' , marker='o', alpha =1)
        plt.title('Sodfium Doublet')
        plt.ylabel('Intensity (Relative Units)')
        plt.xlabel('Wavelength (Micrometers)')
        plt.show() 
       
    CCDMap=dataOut[1:,]



    return CCDMap

#backImage='sky_0deg_2.FIT'

def doPlotOld(CCDMap,CalibPoints=False,Labels=False,BackImage=''):
        x = CCDMap[:,0] 
        z = CCDMap[:,1] 
        Lambda = CCDMap[:,2] 
        Intensity= CCDMap[:,3] 

        colorTable = np.array((wav2RGB(Lambda, Intensity))) 
        
        hdulist = pyfits.open(BackImage)
        imWidth = hdulist[0].header['NAXIS1']
        imHeight = hdulist[0].header['NAXIS2']

        im = pyfits.getdata(BackImage)
        im[im<0]=0
        im /= im.max()
        im = np.sqrt(im) #Remove this line for Hg
#        im = np.sqrt(im) #Remove this line for Hg
#    

#        labels = np.array([0])
#         
#        for line in open('solar2.txt'):
#            Lambda = float(str(line).split()[0]) #Wavelength
#            labels = np.vstack((labels,np.array([Lambda])))
#
#        labels=labels[1:]
        
#        im=mpimg.imread('solar.png')
        
        fig = plt.figure()
        ax1 = fig.add_subplot(111)

        plt.imshow(im,extent=[-imWidth/2 , imWidth/2 , -imHeight/2 , imHeight/2])
        plt.set_cmap(cm.Greys_r)
        ax1.scatter(x, -z ,s=8, color=colorTable , marker='o', alpha =.5)
#        color=colorTable
#        print random.randrange(-30,-10) random()
#        plt.subplots_adjust(bottom = 0.1)
        if Labels==True:
            for label, x, y in zip(Lambda, x, -z):
                plt.annotate(
                    label, 
                    xy = (x, y), xytext = (0,-20),
                    textcoords = 'offset points', ha = 'right', va = 'bottom',
                    bbox = dict(boxstyle = 'round,pad=0.5', fc = 'white', alpha = 0.9),
                    arrowprops = dict(arrowstyle="wedge,tail_width=1.",
                                fc=(0, 0, 1), ec=(1., 1, 1),
                                patchA=None,
                                relpos=(0.2, 0.8),
                                connectionstyle="arc3,rad=-0.1"), size=7)

        plt.ylabel('pixels')
        plt.xlabel('pixels')
        
        
        
        if CalibPoints==True:
            x,y,waveList,xSig,ySig = readCalibrationData(specFile)
            ax1.scatter(x-imWidth/2 , -(y-imHeight/2) ,s=400, color='black', marker='x', alpha=1)



#        plt.plot( x,-z, "o", markersize=7, color=colorTable, markeredgewidth=1,markeredgecolor='g', markerfacecolor='None' )
        
        plt.title('Order Identification')

        plt.axis([-imWidth/2 , imWidth/2 , -imHeight/2 , imHeight/2])
        
        plt.show()

def main(p = [272.31422902, 90.7157937, 59.6543365, 90.21334551, 89.67646101, 89.82098015, 68.0936684,  65.33694031, 1.19265536, 31.50321471, 199.13548823], SEDMode=0, booPlot=False, specFile='c_noFlat_Hg_0deg_10s.txt', intNormalize=1, booDistort=False, booInterpolate=False, booPlotCalibPoints=False, booPlotLabels=False, plotBackImage='c_noFlat_sky_0deg_460_median.fits',booGaussianFit=False):  
    '''
    Computes the projection of n beams of monochromatic light passing through an optical system. 

    Parameters
    ----------
    p : np np.array
        (beam phi, beam theta, prism1 phi, prism1 theta, prism2 phi, prism2 theta, grating phi, grating theta, grating alpha,
         blaze period (microns), focal length(mm), distortion term) <-- optical arrangement
    args : np np.array
        (SEDMode(0=Max, 1=Random, 2=Sun, 3=from specFile, 4=from CalibFile), Plot?, specFile, Normalize intensity? (0=no, #=range), Distort?, Interpolate, PlotCalibPoints) <-- other options
   
    Returns
    -------
    x : np np.array
        x coordinate of the target point 
    y : np np.array
        x coordinate of the target point 
    lambda : np np.array
        wavelength at x,y

    Notes
    -----  
    '''   
    
    global n1, n2, n4, n5, s, l, d, flux
    global allFlux

    #Initial beam
    uiphi = np.radians(p[0])              #'Longitude' with the x axis as 
    uitheta = np.radians(p[1])            #Latitude with the y axis the polar axis
    u=np.array([np.cos(uiphi)*np.sin(uitheta),np.sin(uiphi)*np.sin(uitheta),np.cos(uitheta)])
       
    #Focal length
    fLength = p[10]
    
    #Prism surface 1
    n1phi = np.radians(p[2])   
    n1theta = np.radians(p[3]) 
    n1=np.array([np.cos(n1phi)*np.sin(n1theta),np.sin(n1phi)*np.sin(n1theta),np.cos(n1theta)])
    
    #Prism surface 2
    n2phi = np.radians(p[4])   
    n2theta = np.radians(p[5]) 
    n2=np.array([np.cos(n2phi)*np.sin(n2theta),np.sin(n2phi)*np.sin(n2theta),np.cos(n2theta)])

    #Prism surface 3 (surf #2 on return)
    n4=-n2
    
    #Prism surface 4 (surf #1 on return)
    n5=-n1 
    
    #Grating
    d = p[9]  #blaze period in microns  
    sphi = np.radians(p[6])   
    stheta = np.radians(p[7]) 
    s = np.array([np.cos(sphi)*np.sin(stheta),np.sin(sphi)*np.sin(stheta),np.cos(stheta)]) #component perp to grooves
        
    #Now find two vectors (a and b) perpendicular to s:
    a = np.array([s[1]/np.sqrt(s[0]**2 + s[1]**2), -s[0]/np.sqrt(s[0]**2 + s[1]**2), 0])
    b = np.cross(a,s)
    
    #Create l from given alpha using a and b as basis
    alpha = np.radians(p[8]) 
    l = np.cos(alpha)*a + np.sin(alpha)*b #component along grooves
       
    #Distortion np.array
    K = [0] #p[11]
       
    #Launch grid loop. Creates an array of (x,y,lambda)
    CCDMap = doCCDMapOld(u ,minLambda ,maxLambda ,deltaLambda ,minOrder ,maxOrder ,deltaOrder ,fLength ,stheta, SEDMode, intNormalize,booInterpolate,BackImage=plotBackImage,GaussFit=booGaussianFit) 
    CCDMapX=CCDMap[:,0]
    CCDMapY=CCDMap[:,1]
    CCDMapLambda=CCDMap[:,2]
    CCDMapOrder=CCDMap[:,3]
        
    #Distort
    if K!=0: CCDMapX, CCDMapY = distort(CCDMapX, CCDMapY, K)
    
    #Validates output
#    x=y=Lambda=0
#    if CCDMap.size>0:
#        x=CCDMap[:,0]
#        y=CCDMap[:,1]
#        Lambda=CCDMap[:,2]
    

        
    #Validates new output
#    x=y=Lambda=0
#    if CCDMap.size>0:
#        x=CCDMap[:,0]
#        y=CCDMap[:,1]
#        Lambda=CCDMap[:,2]
        #MJI...
        #p[11]=0
        #newx = x*np.cos(np.pi*p[11]/180) - y*np.sin(np.pi*p[11]/180)
        #newy = y*np.cos(np.pi*p[11]/180) + x*np.sin(np.pi*p[11]/180)
        #CCDMap[:,0]=newx.copy()
        #CCDMap[:,1]=newy.copy()
           
    #Plot
    if booPlot==True:
        doPlot(CCDMap,CalibPoints=booPlotCalibPoints,Labels=booPlotLabels,BackImage=plotBackImage)
    
    return CCDMapX, CCDMapY, CCDMapLambda

def doCCDMapOld(SEDMap, p = [272.31422902, 90.7157937, 59.6543365, 90.21334551, 89.67646101, 89.82098015, 68.0936684,  65.33694031, 1.19265536, 31.50321471, 199.13548823]):
    #Old parameters---->>>>> SEDMode=0, booPlot=False, specFile='c_noFlat_Hg_0deg_10s.txt', intNormalize=1, booDistort=False, booInterpolate=False, booPlotCalibPoints=False, booPlotLabels=False, plotBackImage='c_noFlat_sky_0deg_460_median.fits',booGaussianFit=False):  
    '''
    Computes the projection of n beams of monochromatic light passing through an optical system. 

    Parameters
    ----------
    SEDMap : np.array
        n x 2 np.array with wavelength, Energy
    
    p : np np.array (optional)
        (beam phi, beam theta, prism1 phi, prism1 theta, prism2 phi, prism2 theta, grating phi, grating theta, grating alpha,
         blaze period (microns), focal length(mm), distortion term) <-- optical arrangement
         
    args : np np.array
        (SEDMode(0=Max, 1=Random, 2=Sun, 3=from specFile, 4=from CalibFile), Plot?, specFile, Normalize intensity? (0=no, #=range), Distort?, Interpolate, PlotCalibPoints) <-- other options
   
    Returns
    -------
    x : np np.array
        x coordinate of the target point 
    y : np np.array
        x coordinate of the target point 
    lambda : np np.array
        wavelength at x,y

    Notes
    -----  
    '''   
    
#    global n1, n2, n4, n5, s, l, d, flux
#    global allFlux

    #Initial beam
    uiphi = np.radians(p[0])              #'Longitude' with the x axis as 
    uitheta = np.radians(p[1])            #Latitude with the y axis the polar axis
    Beam=np.array([np.cos(uiphi)*np.sin(uitheta),np.sin(uiphi)*np.sin(uitheta),np.cos(uitheta)])
       
    #Focal length
    fLength = p[10]
    
    #Prism surface 1
    n1phi = np.radians(p[2])   
    n1theta = np.radians(p[3]) 
    n1=np.array([np.cos(n1phi)*np.sin(n1theta),np.sin(n1phi)*np.sin(n1theta),np.cos(n1theta)])
    #Prism surface 2
    n2phi = np.radians(p[4])   
    n2theta = np.radians(p[5]) 
    n2=np.array([np.cos(n2phi)*np.sin(n2theta),np.sin(n2phi)*np.sin(n2theta),np.cos(n2theta)])
    Optics = np.append([[n1,[0,0,0],OpticsBoundary,'nkzfs8',0]], [[n2,[0,0,0],OpticsBoundary,'air',0]], 0)
    
    #Grating
    d = p[9]  #blaze period in microns  
    sphi = np.radians(p[6])   
    stheta = np.radians(p[7]) 
    s = np.array([np.cos(sphi)*np.sin(stheta),np.sin(sphi)*np.sin(stheta),np.cos(stheta)]) #component perp to grooves   
    #Now find two vectors (a and b) perpendicular to s:
    a = np.array([s[1]/np.sqrt(s[0]**2 + s[1]**2), -s[0]/np.sqrt(s[0]**2 + s[1]**2), 0])
    b = np.cross(a,s)
    #Create l from given alpha using a and b as basis
    alpha = np.radians(p[8]) 
    l = np.cos(alpha)*a + np.sin(alpha)*b #component along grooves
    Optics = np.append(Optics, [[s,l,OpticsRGrating,'air',d]], 0)
    
    
    #Prism surface 3 (surf #2 on return)
    n4=-n2
    Optics = np.append(Optics,[[n4,[0,0,0],OpticsBoundary,'nkzfs8',0]], 0)
    
    #Prism surface 4 (surf #1 on return)
    n5=-n1     
    Optics = np.append(Optics, [[n5,[0,0,0],OpticsBoundary,'air',0]], 0)

    #Distortion np.array
    K = [] #p[11]
       
    #Launch grid loop. Creates an array of (x,y,lambda, Intensity, Order)
    CCDX, CCDY, CCDLambda, CCDIntensity, CCDOrder = wt.ccd_loop(SEDMap, Beam , Optics, stheta, fLength) #minLambda ,maxLambda ,deltaLambda ,minOrder ,maxOrder ,deltaOrder ,fLength ,stheta) 
     
    #Distort if any distort data present
    if len(K)!=0: CCDX, CCDY = distort(CCDX, CCDY, K)
        
    return CCDX, CCDY, CCDLambda, CCDIntensity, CCDOrder

def rayTrace(nAir, nPrism, nOrder, Lambda, d, u, n1, n2, n4, n5, s, l, booLog=6):
    ''' Traces a beam through the spectrograph. 
    Spectrograph frame of reference, from the opposite end of the camera looking at the camera
    x=to the right, y=to camera, z=up
    u*=beam, n*=surface normals
    s=grating, perp to the grooves.
    l=grating, parallel to the grooves.
    d=blaze period'''
    
    if booLog==0:
        print 's={'+str("%0.10f"%s[0])+', '+str("%0.10f"%s[1])+', '+str("%0.10f"%s[2])+'}'  
        print 'l={'+str("%0.10f"%l[0])+', '+str("%0.10f"%l[1])+', '+str("%0.10f"%l[2])+'}'  
        print 'n1={'+str("%0.10f"%n1[0])+', '+str("%0.10f"%n1[1])+', '+str("%0.10f"%n1[2])+'}'
        print 'n2={'+str("%0.10f"%n2[0])+', '+str("%0.10f"%n2[1])+', '+str("%0.10f"%n2[2])+'}'
        print 'u={'+str("%0.10f"%u[0])+', '+str("%0.10f"%u[1])+', '+str("%0.10f"%u[2])+'}'  
                           
    """Vector transform due to first surface"""
    u = Snell3D(nAir, nPrism, u, n1)
    
    if booLog==1:
        tempu=u*10
        print '{{'+str("%0.10f"%tempu[0])+', '+str("%0.10f"%tempu[1])+', '+str("%0.10f"%tempu[2])+'},'  
        
        
    """Vector transform due to second surface"""
    u = Snell3D(nPrism, nAir, u, n2)
    
    if booLog==1:
        tempu+=u*50
        print '{'+str("%0.10f"%tempu[0])+', '+str("%0.10f"%tempu[1])+', '+str("%0.10f"%tempu[2])+'},'   
    
    """grating dispersion"""              
    u, isValid = Grating(u, s, l, nOrder, Lambda, d)
    
    if booLog==1:
        tempu+=u*50
        print '{'+str("%0.10f"%tempu[0])+', '+str("%0.10f"%tempu[1])+', '+str("%0.10f"%tempu[2])+'},'      

    
    if isValid:
        """Vector transform due to third surface"""
        u = Snell3D(nAir, nPrism, u, n4)
        
        if booLog==1:
            tempu+=u*10
            print '{'+str("%0.10f"%tempu[0])+', '+str("%0.10f"%tempu[1])+', '+str("%0.10f"%tempu[2])+'},'               
        
        """Vector transform due to fourth surface"""
        u = Snell3D(nPrism, nAir, u, n5)
        
        if booLog==1:
            tempu+=u*150
            print '{'+str("%0.10f"%tempu[0])+', '+str("%0.10f"%tempu[1])+', '+str("%0.10f"%tempu[2])+'}},'     
        elif booLog==2:
            print '{'+str("%0.10f"%u[0])+', '+str("%0.10f"%u[1])+', '+str("%0.10f"%u[2])+'},'   

    return u, isValid

def nkzfs8(Lambda):
    #### absorved by n()
    '''Function that calculates the refractive index of the prism for a given wavelength'''
    x = float(Lambda)
    n = np.sqrt(1 + (1.62693651*x**2)/(x**2-0.010880863) + 0.24369876*x**2/(x**2-0.0494207753) + 1.62007141*x**2/(x**2-131.009163) )
    
    return n

def analyze_image(image_filename='test.fits', image_map_filename='image_map.txt'):
    iraf.noao(_doprint=0)     # load noao
    iraf.digiphot(_doprint=0) # load digiphot
    iraf.apphot(_doprint=0)   # load apphot
    

    
    #iraf.daofind.setParam('image',FitsFileName)        #Set ImageName
    iraf.daofind.setParam('verify','no')            #Don't verify
    iraf.daofind.setParam('interactive','no')        #Interactive
    iraf.daofind.setParam('datamin','10000')        #Min Good Value
    iraf.daofind.setParam('fwhmpsf','2.5')        #Interactive   
#    iraf.daofind.setParam('interactive','no')        #Interactive

    
    check_if_file_exists(image_map_filename)
    try: iraf.daofind(image = image_filename, output = image_map_filename)
    except Exception: return 0
#    brightest_star_info = find_brightest_star(outfile)
    return

def find_brightest_star(readinfile):
    try: starfile = open(readinfile)
    except Exception: return 'ERROR' # <-- change this to returning a number
    startemp = starfile.readlines()
    brighteststar = 50
    xpixel = 0
    ypixel = 0
    for lines in startemp:
        if lines[0][0] != '#': #don't want the comments
            linetemp = str.split(lines)
            #print linetemp
            if 1: #float(linetemp[2]) < brighteststar:
                starmag = 1 #float(linetemp[2])
                xpixel = float(linetemp[0])
                ypixel = float(linetemp[1])
                starsharp = float(linetemp[3])
    return [starmag, xpixel, ypixel, starsharp]

def do_read_calib(outputFileName, image_filename='test.fits', analyze=True):
    '''
    Extracts peaks with daofind
    Imports found points
    Plots found points on image
       
    Parameters
    ----------
    image_filename : string
        Name of the calibration image 

    analyze : boolean
        Run the pyraf analysis (reads the last output otherwise)

    Returns
    -------
    nottin
      
    Notes
    -----
    '''      
    if analyze: ic.analyze_image(image_filename)
    
    #Loads from the iraf output file
    imageMapX, imageMapY = wt.load_image_map()
    
    imageMapLambda_match = imageMapLambda = np.zeros(len(imageMapX))
    
    #Create SEDMap from Mercury emission
    SEDMap = do_sed_map(SEDMode=SED_MODE_CALIB, specFile='Hg_5lines_double.txt')
    
    #Create the model based on default parameters
    CCDX, CCDY, CCDLambda, CCDIntensity, CCDOrder = do_ccd_map(SEDMap)  
    
    #Create list of 
    #todo turn this into a 2D array calculation
    for i in range(len(imageMapX)):
    
        distance_array = np.sqrt((CCDX-imageMapX[i])**2+(CCDY-imageMapY[i])**2)
        closest_point=np.min(distance_array)
        closest_point_index=np.where(distance_array==closest_point)[0][0]       
        imageMapLambda[i] = CCDLambda[closest_point_index]
        
        lambda_distance= abs(SEDMap[SEDMapLambda]-imageMapLambda[i])
        lambda_closest_point=np.min(lambda_distance)
        lambda_closest_point_index=np.where(lambda_distance==lambda_closest_point)[0][0]
        imageMapLambda_match[i] = SEDMap[SEDMapLambda][lambda_closest_point_index]
    
    #Create output file with calibration data
    #todo, add sigmas    
    f = open(outputFileName,'w')
    for i in range(len(imageMapX)):
        out_string=str(imageMapX[i])+' '+str(imageMapY[i])+' '+str(imageMapLambda_match[i])+' 1 1\n'
        f.write(out_string) 
    f.close()
    
    #Plot (probably remove)
    hdulist = pyfits.open(image_filename)
    imWidth = hdulist[0].header['NAXIS1']
    imHeight = hdulist[0].header['NAXIS2']
    im = pyfits.getdata(image_filename)  
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    plt.imshow(im,extent=[-imWidth/2 , imWidth/2 , -imHeight/2 , imHeight/2])
    plt.set_cmap(cm.Greys_r)
    ax1.scatter(imageMapX-imWidth/2, -(imageMapY-imHeight/2) ,s=40, color="red" , marker='o', alpha = 0.3)
    plt.title(str(len(imageMapX))+' point(s) found')
    plt.axis([-imWidth/2 , imWidth/2 , -imHeight/2 , imHeight/2])
    plt.show()
    
    return

def findFit(calibrationFile, p_try=[ 271.92998622,   91.03999719,   59.48997316,   89.83000496,   89.37002499,   89.79999531,   68.03002976,   64.9399939,     1.15998754,   31.52736851,  200.00000005],factor_try=1,diag_try=[1,1,1,1,1,1,1,1,1,.1,1]):
    '''
    Wrapper for reading the calibration file, and launching the fitting function
       
    Parameters
    ----------
    calibrationFile : string
        Name of the file with the data from the spectrograph

    Returns
    -------
    fit : np np.array
        1 x 12 np np.array with fitted arguments (p np.array)
      
    Notes
    -----
    '''  
    #old mainArgs=['4','0',calibrationFile,'0','0']
    #x,y, wavelist are the positions of the peaks in calibrationFile.
    #x,y,waveList,xsig,ysig = readCalibrationData(calibrationFile)
    #fit is the output, which is the ideal p vector.

    fit = leastsq(main_errors,p_try, args=[4,False,calibrationFile,0,False,False,False,True,'c_noFlat_sky_0deg_460_median.fits',False], full_output=True, factor=factor_try, diag=diag_try)


    return fit
