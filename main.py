#Imports
import matplotlib.pyplot as plt     #python/matlab
import pylab
import random                       #random generator package
import pyfits
import os
import numpy as np
import matplotlib.cm as cm

#least square package
from scipy.optimize.minpack import leastsq

#interpolate function
from scipy import interpolate

#Astro Libraries
from astLib import astSED           

#Math package 
from math import cos, sin, acos, asin, pi, atan, degrees, radians

'''
Initial Parameters-------------------------------------------------------------- 
Wavelength'''
minLambda=0.380   #min wavelength
maxLambda=0.780    #max wavelength
deltaLambda=0.0001    #step interval
maxLambda+=deltaLambda

#Can plot orders from  146 to 73 (about 390 to 795nm)
minOrder=90
maxOrder=90
deltaOrder=1
maxOrder+=deltaOrder

# still don't know what this is 
booLog=6 
# Pixel size in microns
pixelSize= 5.4 

def main_errors(p, mainArgs):
    #main will return these vectors in a random order. 
    #We assume that there are no rounding errors (probably a hack?)
    #and will use a floating point == to identify the (x,y) corresponding
    #to each wavelength input.
    #NB A much better idea would be to re-write main without vstack but
    #instead with a neatly created output structure that is defined at the start.
    #i.e. when going from SEDMap to SEDMapLoop the information on which line in
    #the input each wavelength came from was lost.
    #print p
    x,y,waveList,xSig,ySig = readCalibrationData(mainArgs[2])

    hdulist = pyfits.open('../c_noFlat_sky_0deg_460_median.fits')
    imWidth = hdulist[0].header['NAXIS1']
    imHeight = hdulist[0].header['NAXIS2']
    
    x=x-imWidth/2
    y=y-imHeight/2
    
    x_model, y_model, Lambda = main(p, mainArgs)
    
    #Loop through the input wavelengths, and find the closest output.
    #The best model fits in x and y (out of 2 options) is called x_best and y_best
    x_best = x.copy()
    y_best = y.copy()
    for k in range(0,len(waveList)):
        ix, = np.where(waveList[k] == Lambda)
        if (ix.size == 0):
         x_best[k]=0
         y_best[k]=0
        else:
         best = ix[np.argmin(np.abs(y_model[ix] - y[k]))]
         #The following lines have a potential bug because they assume there is only one "best"
         x_best[k] = x_model[best]
         y_best[k] = y_model[best]
#    print x_best + 1663
#    print y_best + 1252
#    print x_best - x
#    print y_best - y
#    print np.hstack([(x_best-x)/xSig,(y_best - y)/ySig])

    return np.hstack([(x_best-x)/xSig,(y_best - y)/ySig])

#[250.6, 64.2, 58.3, 77.9, 89.5, 85.5, 58.5, 61.7, 1.6, 33.7, 193.5]

def main(p = [272.31422902, 90.7157937, 59.6543365, 90.21334551, 89.67646101, 89.82098015, 68.0936684,  65.33694031, 1.19265536, 31.50321471, 199.13548823], 
         args=['0','0','../Hg_5lines_double.txt','1','0']):
    '''
    Compute the projection of n beams of monochromatic light passing through an optical system. 

    Parameters
    ----------
    p : np np.array
        (beam phi, beam theta, prism1 phi, prism1 theta, prism2 phi, prism2 theta, grating phi, grating theta, grating alpha,
         blaze period (microns), focal length(mm), distortion term) <-- optical arrangement
    args : np np.array
        (SEDMode(0=Max, 1=Random, 2=Sun, 3=from specFile, 4=from CalibFile), Plot?, specFile, Normalize intensity? (0=no, #=range), Distort?) <-- other options
   
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
    
    global n1, n2, n4, n5, s, l, d, flux, booInterpolate
    global booPlot, specFile, booPlotCalibPoints
    
    booInterpolate=1
    booPlotCalibPoints=0
    
    #Args breakdown
    SEDMode = int(args[0])
    booPlot = int(args[1])
    specFile = args[2]
    intNormalize = int(args[3])
    booDistort = int(args[4])
    
    #Initial beam
    uiphi = radians(p[0])              #'Longitude' with the x axis as 
    uitheta = radians(p[1])            #Latitude with the y axis the polar axis
    u=np.array([cos(uiphi)*sin(uitheta),sin(uiphi)*sin(uitheta),cos(uitheta)])
       
    #Focal length
    fLength = p[10]
    
    #Prism surface 1
    n1phi = radians(p[2])
    n1theta = radians(p[3])
    n1=np.array([cos(n1phi)*sin(n1theta),sin(n1phi)*sin(n1theta),cos(n1theta)])
    
    #Prism surface 2
    n2phi = radians(p[4])
    n2theta = radians(p[5]) 
    n2=np.array([cos(n2phi)*sin(n2theta),sin(n2phi)*sin(n2theta),cos(n2theta)])

    #Prism surface 3 (surf #2 on return)
    n4=-n2
    
    #Prism surface 4 (surf #1 on return)
    n5=-n1 
    
    #Grating
    d = p[9]  #blaze period in microns  
    sphi = radians(p[6])   
    stheta = radians(p[7]) 
    s = np.array([cos(sphi)*sin(stheta),sin(sphi)*sin(stheta),cos(stheta)]) #component perp to grooves
        
    #Now find two vectors perpendicular to s:
    a = np.array([s[1]/np.sqrt(s[0]**2 + s[1]**2), -s[0]/np.sqrt(s[0]**2 + s[1]**2), 0])
    b = np.cross(a,s)
    
    #Create l from given alpha using a and b as basis
    alpha = radians(p[8])
    l = cos(alpha)*a + sin(alpha)*b #component along grooves
       
    #Distortion np.array
    K = [0] #p[11]
       
    #Launch grid loop. Creates an array of (x,y,lambda)
    CCDMap = doCCDMap(u ,minLambda ,maxLambda ,deltaLambda ,minOrder ,maxOrder ,deltaOrder ,fLength ,stheta, SEDMode, intNormalize) 
        
        
    #Distort
    if booDistort==1:
        x=CCDMap[:,0]
        y=CCDMap[:,1]
        CCDMap[:,0], CCDMap[:,1] = distort(x, y, K)
    
    #Validates output
    x=y=Lambda=0
    if CCDMap.size>0:
        x=CCDMap[:,0]
        y=CCDMap[:,1]
        Lambda=CCDMap[:,2]
    

        
    #Validates new output
    x=y=Lambda=0
    if CCDMap.size>0:
        x=CCDMap[:,0]
        y=CCDMap[:,1]
        Lambda=CCDMap[:,2]
        #MJI...
        #p[11]=0
        #newx = x*cos(radians(p[11])) - y*sin(radians(p[11]))
        #newy = y*cos(radians(p[11])) + x*sin(radians(p[11]))
        #CCDMap[:,0]=newx.copy()
        #CCDMap[:,1]=newy.copy()
           
    #Plot
    if booPlot==1:
        doPlot(CCDMap,'../mflat.fits')#,'../c_noFlat_Hg_0deg_10s.fits')
    
    return x, y, Lambda

def extractOrder(x,y,image):
    
    
    flux=np.zeros(len(y))
    flux2=np.zeros(len(y))
    flux3=np.zeros(len(y))

    tempIm=np.zeros((len(y),60))

    hdulist = pyfits.open(image)
    imWidth = hdulist[0].header['NAXIS1']
    imHeight = hdulist[0].header['NAXIS2']
    
    im = pyfits.getdata(image)  

    x=x+imWidth/2
    y=y+imHeight/2
    width=6
    
    for k in range(0,len(y)):
        ###mikes'
#        x_int = round(x[k])
#        #if e.g. x = 10.4, this goes from 5 to 15, and xv at the center is 0.4
#        in_image_temp = im[y[k],x_int-5:x_int+6]
#        xv = np.arange(-5,6) - x[k] + x_int
#        flux[k] =  np.sum(in_image_temp * np.exp(-(xv/3.5)**4))
#        flux2[k] = np.sum(in_image_temp)
#        flux3[k] = np.sum(in_image_temp)- np.sum(in_image_temp * np.exp(-(xv/3.5)**4))
        
        #Carlos' trial
        x_int = int(x[k])
        res=x[k]-x_int- 0.5
        sum_in_image_temp = np.sum(im[y[k],x_int-width/2:x_int+width/2])
        flux[k] = sum_in_image_temp - ((res*im[y[k],x[k]+width/2+1])+ ((1-res)*im[y[k],x[k]-width/2]))
        flux2[k] = ((res*im[y[k],x[k]+width/2+2])+ ((1-res)*im[y[k],x[k]-width/2]))
        flux3[k] = sum_in_image_temp 
        
#    print xv #Mike's check
#       in_image=np.reshape(np.tile(in_image_temp,in_image_temp.shape[0]),(in_image_temp.shape[0],in_image_temp.shape[0]))
#        shift=(0,x[k]-x[0])
#        out_image = fftshift(in_image,shift)
        
#        print in_image[0]
#        print out_image[0]
        
#        flux[k]=np.sum(out_image[0])  
        
#        tempIm[k,:]=im[y[k]+imHeight/2,x[k]-30:x[k]+30]

#    plt.imshow(tempIm)
#    plt.set_cmap(plt.cm.Greys_r)
#    plt.show()
             
    return flux, flux2 +np.average(flux), flux3

def fftshift1D(inImage, shift):
    '''
    This program shifts an image by sub-pixel amounts.
       
    Parameters
    ----------
    inImage :  image
        Input image
    shift : array
        (x,y) pixel shift
    Returns
    -------
    outImage : Image
        Shifted Image
      
    Notes
    -----
    '''  
    
## An example - shift a rough Gaussian by 0.5 pixels.
#in_image = np.array([[0,0,0.06,0,0],[0,0.25,0.5,0.25,0],[0.06,0.5,1,0.5,0.06],[0,0.25,0.5,0.25,0],[0,0,0.06,0,0]])
#out_image = fftshift(in_image,(0,0))
#print out_image

    ftin = np.fft.fft(inImage)
    sh = inImage.shape
    
    #The following line makes a meshgrid np.array as floats. Not sure if there is a neater way.
    xy = np.mgrid[0:sh[0]] + 0.0
    xy[:] = (((xy[0,:] + sh[0]/2) % sh[0]) - sh[0]/2)/float(sh[0])
    #xy[1,:,:] = (((xy[1,:,:] + sh[1]/2) % sh[1]) - sh[1]/2)/float(sh[1])

    return np.real(np.fft.ifft(ftin*np.exp(np.complex(0,-2*np.pi)*(xy[0,:,:]*shift[0] + xy[1,:,:]*shift[1])))) 

def fftshift(inImage, shift):
    '''
    This program shifts an image by sub-pixel amounts.
       
    Parameters
    ----------
    inImage :  image
        Input image
    shift : array
        (x,y) pixel shift
    Returns
    -------
    outImage : Image
        Shifted Image
      
    Notes
    -----
    '''  
    
## An example - shift a rough Gaussian by 0.5 pixels.
#in_image = np.array([[0,0,0.06,0,0],[0,0.25,0.5,0.25,0],[0.06,0.5,1,0.5,0.06],[0,0.25,0.5,0.25,0],[0,0,0.06,0,0]])
#out_image = fftshift(in_image,(0,0))
#print out_image

    ftin = np.fft.fft2(inImage)
    sh = inImage.shape
    
    #The following line makes a meshgrid np.array as floats. Not sure if there is a neater way.
    xy = np.mgrid[0:sh[0],0:sh[1]] + 0.0
    xy[0,:,:] = (((xy[0,:,:] + sh[0]/2) % sh[0]) - sh[0]/2)/float(sh[0])
    xy[1,:,:] = (((xy[1,:,:] + sh[1]/2) % sh[1]) - sh[1]/2)/float(sh[1])
    
    return np.real(np.fft.ifft2(ftin*np.exp(np.complex(0,-2*np.pi)*(xy[0,:,:]*shift[0] + xy[1,:,:]*shift[1])))) 
    
def doCCDMap(u, minLambda, maxLambda, deltaLambda, minOrder, maxOrder, deltaOrder, fLength, stheta, SEDMode, intNormalize):
    
    dataOut=np.zeros((1,5))
    
    #Loads SEDMap based on selection. 
    SEDMap = doSEDMap(SEDMode, minLambda, maxLambda, deltaLambda, intNormalize)
    blaze_angle= stheta #Approximately atan(2)

    '''Main loop
    Navigates orders within the range given
    For each order navigates the list of wavelenghts in SEDMap constrained by +/- FSP/2'''    
    for nOrder in range(minOrder, maxOrder, deltaOrder):
            
        LambdaBlMin = LambdaBlMax=0
        #LambdaBlMin = 2*d*sin(blaze_angle)/(nOrder+0.5) #Was 0.5
        LambdaBlMin = 2*d*sin(blaze_angle)/(nOrder+1) #Was 0.5
#        print LambdaBlMin       
        #LambdaBlMax = 2*d*sin(blaze_angle)/(nOrder-1.25) #Was 0.5
        LambdaBlMax = 2*d*sin(blaze_angle)/(nOrder-1) #Was 0.5
#        print LambdaBlMax
 
        #SEDMapLoop is an array of wavelengths (paired with intensity) called SEDMap that are in this particular order.
        SEDMapLoop=SEDMap.copy()
        
        #constrain by +/- FSP (was FSP/2)
        SEDMapLoop = SEDMapLoop[SEDMapLoop[:,0]>=LambdaBlMin]
        if SEDMapLoop.shape[0]>0:       
            SEDMapLoop = SEDMapLoop[SEDMapLoop[:,0]<=LambdaBlMax]     

        #loop lambda for current order
        for Lambda,inI in SEDMapLoop: 
                                    
            nPrism = nkzfs8(Lambda)
            nAir = n(Lambda)
            
            #Computes the unit vector that results from the optical system for a given wavelength and order
            v, isValid = rayTrace(nAir, nPrism, nOrder, Lambda, d, u, n1, n2, n4, n5, s, l)
            
            if isValid: #no errors in calculation
                x=v[0]*fLength*1000/pixelSize # x-coord in focal plane in pixels
                z=v[2]*fLength*1000/pixelSize # z-coord in focal plane in pixels
    
                '''Appends data table
                x,z=pojected coordinates
                inI, outI = intensities, inI=from source (file, random,...) 0.0 to 1.0'''
                outI=Intensity(Lambda, minLambda, maxLambda)              
                dataOut= np.vstack((dataOut,np.array([x,z, Lambda, inI*outI ,nOrder]))) 
        
        #Order extraction        
        if booInterpolate==1:
            xPlot=dataOut[dataOut[:,4]==nOrder][:,0]
            yPlot=dataOut[dataOut[:,4]==nOrder][:,1]
            LambdaPlot=dataOut[dataOut[:,4]==nOrder][:,2]
            
            fLambda = interpolate.interp1d(yPlot, LambdaPlot)
            fX = interpolate.interp1d(yPlot, xPlot, 'quadratic', bounds_error=False)
            #print y,f(np.trunc(y)) 
            
    #        image ='../c_noFlat_Hg_0deg_10s.fits'
    #        image ='../c_noFlat_sky_0deg_460_median.fits'
#            image ='../simple_flat.fits'
            image ='../mflat.fits'
    
            hdulist = pyfits.open(image)
            imWidth = hdulist[0].header['NAXIS1']
            imHeight = hdulist[0].header['NAXIS2']
            
            newY=np.arange(-imHeight/2,imHeight/2)
            
            newX=fX(newY)
            
            nanMap= np.isnan(newX)
            newX=newX[-nanMap]
            newY=newY[-nanMap]
                    
            flux,flux2,flux3 = extractOrder(newX,newY,image)
            
            fig = plt.figure()
            ax1 = fig.add_subplot(111)
        
            ax1.plot(fLambda(newY),flux)
            ax1.plot(fLambda(newY),flux2)
            ax1.plot(fLambda(newY),flux3)

            #plt.plot(fLambda(newY),flux)
            plt.show()
        
    CCDMap=dataOut[1:,]
#    if dataOut.size>5: #temporary hack for leastsq fit       
#        CCDMap=dataOut[1:,]
#    else:
#        CCDMap=np.array([[0,0,0,0,0],[0,0,0,0,0],[0,0,0,0,0],[0,0,0,0,0],[0,0,0,0,0],[0,0,0,0,0],[0,0,0,0,0]])
        
    return CCDMap

#backImage='../sky_0deg_2.FIT'
   
def x(Lambda, p = [272.31422902, 90.7157937, 59.6543365, 90.21334551, 89.67646101, 89.82098015, 68.0936684,  65.33694031, 1.19265536, 31.50321471, 199.13548823]): #,args=['2','1','../Hg_5lines_double.txt','1','0']):
    
    #Args breakdown
#    SEDMode = int(args[0])
#    booPlot = int(args[1])
#    specFile = args[2]
#    intNormalize = int(args[3])
#    booDistort = int(args[4])
    
    #Initial beam
    uiphi = radians(p[0])              #'Longitude' with the x axis as 
    uitheta = radians(p[1])            #Latitude with the y axis the polar axis
    u=np.array([cos(uiphi)*sin(uitheta),sin(uiphi)*sin(uitheta),cos(uitheta)])
       
    #Focal length
    fLength = p[10]
    
    #Prism surface 1
    n1phi = radians(p[2])
    n1theta = radians(p[3])
    n1=np.array([cos(n1phi)*sin(n1theta),sin(n1phi)*sin(n1theta),cos(n1theta)])
    
    #Prism surface 2
    n2phi = radians(p[4])
    n2theta = radians(p[5])
    n2=np.array([cos(n2phi)*sin(n2theta),sin(n2phi)*sin(n2theta),cos(n2theta)])

    #Prism surface 3 (surf #2 on return)
    n4=-n2
    
    #Prism surface 4 (surf #1 on return)
    n5=-n1 
    
    #Grating
    d = p[9]  #blaze period in microns  
    sphi = radians(p[6])
    stheta = radians(p[7])
    s = np.array([cos(sphi)*sin(stheta),sin(sphi)*sin(stheta),cos(stheta)]) #component perp to grooves
        
    #Now find two vectors perpendicular to s:
    a = np.array([s[1]/np.sqrt(s[0]**2 + s[1]**2), -s[0]/np.sqrt(s[0]**2 + s[1]**2), 0])
    b = np.cross(a,s)
    
    #Create l from given alpha using a and b as basis
    alpha = radians(p[8])
    l = cos(alpha)*a + sin(alpha)*b #component along grooves
       
    #Distortion np.array
    K = [0] #p[11]
    
    nPrism = nkzfs8(Lambda)
    nAir = n(Lambda)

    blaze_angle= stheta #Approximately atan(2)

    dataOut=np.zeros((1,2))
 
    for nOrder in range(minOrder, maxOrder, deltaOrder):
        LambdaBlMin = LambdaBlMax=0
        LambdaBlMin = 2*d*sin(blaze_angle)/(nOrder+1) #Was 0.5
        LambdaBlMax = 2*d*sin(blaze_angle)/(nOrder-1) #Was 0.5
        
        if (LambdaBlMin<=Lambda<=LambdaBlMax):
            
            v, isValid = rayTrace(nAir, nPrism, nOrder, Lambda, d, u, n1, n2, n4, n5, s, l)
        
            tempx=v[0]*fLength*1000/pixelSize
            dataOut= np.vstack((dataOut,np.array([tempx,nOrder]))) 

    CCDMap=dataOut[1:,]
    
    #z=v[2]*fLength*1000/pixelSize # z-coord in focal plane in pixels
    """Remove rows outside the CCD range"""
    CCDMap = CCDMap[CCDMap[:,0]>=-3352/2]     
    CCDMap = CCDMap[CCDMap[:,0]<=3352/2]   
     
    return CCDMap

def doPlot(CCDMap, backImage='../c_noFlat_sky_0deg_460_median.fits'):
        x = CCDMap[:,0] 
        z = CCDMap[:,1] 
        Lambda = CCDMap[:,2] 
        Intensity= CCDMap[:,3] 

        colorTable = np.array((wav2RGB(Lambda, Intensity))) 
        
        hdulist = pyfits.open(backImage)
        imWidth = hdulist[0].header['NAXIS1']
        imHeight = hdulist[0].header['NAXIS2']

        im = pyfits.getdata(backImage)
        im /= im.max()
        #im = np.sqrt(im) #Remove this line for Hg
    

        
        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        
        plt.imshow(im,extent=[-imWidth/2 , imWidth/2 , -imHeight/2 , imHeight/2])
        plt.set_cmap(cm.Greys_r)
        ax1.scatter(x, -z ,s=4, color=colorTable, marker='o')
        
        if booPlotCalibPoints==1:
            x,y,waveList,xSig,ySig = readCalibrationData('../c_noFlat_Hg_0deg_10s.txt')
            ax1.scatter(x-imWidth/2 , -(y-imHeight/2) ,s=10, color='r', marker='s')



#        plt.plot( x,-z, "o", markersize=7, color=colorTable, markeredgewidth=1,markeredgecolor='g', markerfacecolor='None' )

        plt.axis([-imWidth/2 , imWidth/2 , -imHeight/2 , imHeight/2])
        
        plt.show()

def rayTrace(nAir, nPrism, nOrder, Lambda, d, u, n1, n2, n4, n5, s, l):
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
    u, isValid = Grating(u, l, s, nOrder, Lambda, d)
    
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

def Snell3D(n_i,n_r,u,n):
    """Computes the new direction of a vector whne changing medium. 
    n_i, n_r = incident and refractive indices"""
 
    u_p = u - np.dot(u,n)*n
    u_p /= np.linalg.norm(u_p)
    
    theta_i = acos(np.dot(u,n))
    
    if n_i*sin(theta_i)/n_r<=1:
        theta_f = asin(n_i*sin(theta_i)/n_r)
        u = u_p*sin(pi-theta_f) + n*cos(pi-theta_f)    

    return u       

def Grating(u, l, s, nOrder, Lambda, d):
    """Computes the new direction of a vector when hitting a grating."""
    
    isValid=False

    n = np.cross(s, l) 
    u_l = np.dot(u, l)
    u_s = np.dot(u, s) + nOrder*Lambda/d  
    if (1-u_l**2 -u_s**2)>=0: 
        u_n = np.sqrt(1- u_l**2 - u_s**2)
        u = u_l*l + u_s*s + u_n*n
        isValid=True

    return u, isValid     

def doSEDMap(SEDMode, minLambda, maxLambda, deltaLambda, intNormalize):
    '''
    Loads the SED map
       
    Parameters
    ----------
    SEDMode : int
        Mode for the creation of the SEDMap
        0=Max, 1=Random, 2=Sun, 3=from specFile, 4=from Calibration file
        
    minLambda : np.float32
        Lower limit for SEDMap.

    maxLambda : np.float32
        Higher limit for SEDMap.

    deltaLambda : np.float32
        Step between wavelengths.
    
    intNormalize : integer
        if !=0, it normalizes to intNormalize value
        
    Returns
    -------
    SEDMap : np np.array
        n x 2 np np.array with wavelength, Energy
      
    Notes
    -----
    '''  
    if SEDMode==0: #Flat
        SEDMap = np.column_stack((np.arange(minLambda, maxLambda, deltaLambda),np.ones(np.arange(minLambda, maxLambda, deltaLambda).size)))

    elif SEDMode==1: #Random
        SEDMap = np.array([0,0])
        
        for Lambda in range(minLambda, maxLambda + deltaLambda, deltaLambda):
#            Intensity=int()
#            newItem=np.np.array([Lambda,random.random(0.0,1.0)])
            SEDMap = np.vstack((SEDMap,np.array([Lambda,random.random(0.0,1.0)])))
            
        SEDMap = SEDMap[1:,]
                 
    elif SEDMode==2: #Sun
        
        sol = astSED.SOL        
        
        tempA=sol.wavelength.transpose()*1e-4
        tempB=sol.flux.transpose()            
        SEDMap = np.column_stack((tempA, tempB))    
        
        """Remove rows outside the wavelength range"""
        SEDMap = SEDMap[SEDMap[:,0]>=minLambda]     
        SEDMap = SEDMap[SEDMap[:,0]<=maxLambda]     
                                  
    elif SEDMode==3: #From file
        SEDMap = np.array([0,0])
         
        for line in open(specFile):
            Lambda = float(str(line).split()[0]) #Wavelength
            I = float(str(line).split()[1])       #Intensity
            SEDMap = np.vstack((SEDMap,np.array([Lambda,I])))

        SEDMap=SEDMap[1:,]  
         
    elif SEDMode==4: #From calibration file
        SEDMap = np.array([0,0])
         
        for line in open(specFile):
            Lambda = float(str(line).split()[2]) #Wavelength
            I = 1      #Intensity
            SEDMap = np.vstack((SEDMap,np.array([Lambda,I])))
            
        SEDMap=SEDMap[1:,]  
                
    """Normalize the intensity"""   
    if intNormalize!=0:    
        fluxRange=(max(SEDMap[:,1])-min(SEDMap[:,1]))
        
        if fluxRange==0:
            SEDMap = np.column_stack((SEDMap[:,0], np.ones(SEDMap[:,0].size)))  
        else:
            SEDMap = np.column_stack((SEDMap[:,0], (SEDMap[:,1]-min(SEDMap[:,1]))/(fluxRange+1) ))
       
    return SEDMap

def readCalibrationData(calibrationFile):
    
    CalibrationMap = np.array([0,0,0,0,0])
     
    for line in open(calibrationFile):
        x = float(str(line).split()[0]) #x
        y = float(str(line).split()[1]) #y
        Lambda = float(str(line).split()[2]) #Wavelength
        xSig = float(str(line).split()[3]) #xSig
        ySig = float(str(line).split()[4]) #ySig
        CalibrationMap = np.vstack((CalibrationMap,np.array([x,y,Lambda,xSig,ySig])))
    
    CalibrationMap=CalibrationMap[1:,]  
    

    
    return CalibrationMap[:,0], CalibrationMap[:,1], CalibrationMap[:,2], CalibrationMap[:,3], CalibrationMap[:,4]

def nkzfs8(Lambda):
    
    x = float(Lambda)
    n = np.sqrt(1 + (1.62693651*x**2)/(x**2-0.010880863) + 0.24369876*x**2/(x**2-0.0494207753) + 1.62007141*x**2/(x**2-131.009163) )
    
    return n
    
def column(matrix, i):
    
    a=[row[i] for row in matrix]
    
    return a    

def wav2RGB(Lambda, Intensity):
    
    """Converts Lambda into RGB"""
    out=np.array([0,0,0])
    
    for i in range(Lambda.size):
    
        w = Lambda[i]
        I = Intensity[i]
    
        # colour
        if w >= .380 and w < .440:
            R = -(w - .440) / (.440 - .350)
            G = 0.0
            B = 1.0
        elif w >= .440 and w < .490:
            R = 0.0
            G = (w - .440) / (.490 - .440)
            B = 1.0
        elif w >= .490 and w < .510:
            R = 0.0
            G = 1.0
            B = -(w - .510) / (.510 - .490)
        elif w >= .510 and w < .580:
            R = (w - .510) / (.580 - .510)
            G = 1.0
            B = 0.0
        elif w >= .580 and w < .645:
            R = 1.0
            G = -(w - .645) / (.645 - .580)
            B = 0.0
        elif w >= .645 and w <= .780:
            R = 1.0
            G = 0.0
            B = 0.0
        else:
            R = 1.0
            G = 1.0
            B = 1.0
    
        # intensity correction
        if w >= .3800 and w < .4200:
            SSS = 0.3 + 0.7*(w - .3500) / (.4200 - .3500)
        elif w >= .4200 and w <= .7000:
            SSS = 1.0
        elif w > .7000 and w <= .7800:
            SSS = 0.3 + 0.7*(.7800 - w) / (.7800 - .7000)
        else:
            SSS = 1.0
        SSS *= (I)

        out=np.vstack((out,np.array( [float(SSS*R), float(SSS*G), float(SSS*B)]))) 
        
    return out[1:,]

def Intensity(Lambda, minLambda, maxLambda):
    '''
    Retrieves or calculates the expected relative intensity based on distance from the central lambda value
       
    Parameters
    ----------
    Lambda : np.float32
        Wavelength.
    minLambda : np.float32
        Lower end of wavelength range.
    maxLambda : np.float32
        Higher end of wavelength range. 
   
    Returns
    -------
    z : np.float32
        result 0 to 1.
      
    Notes
    -----
    
    '''

#    x = 0.5*(float(Lambda) - 0.5*(float(maxLambda) - float(minLambda)))/(float(maxLambda) - float(minLambda))
    x = (((float(Lambda) - float(minLambda))/(float(maxLambda) - float(minLambda)))-0.5)*2
    if x!=0:
        z=sin(x*pi)/(x*pi)
     

#    print x,z
                
#        """Look for wavelength"""
#        EMax = max(spectrum_dict.iteritems(), key=operator.itemgetter(1))[1]   
#        a= spectrum_dict[str(Lambda)]            
#        if (1==2):
#            E=random.random() * (EMax-EMin)+EMin            
#        else:
#            j = spectrum_dict
#            topLambda = max([i for i in j if i >= 5])
#            topLambda = max(spectrum_dict.iteritems(), key=operator.itemgetter(1))[1]
#            bottomLambda=max(spectrum_dict.iteritems(), key=operator.itemgetter(1))[1]
#            topE=max(spectrum_dict.iteritems(), key=operator.itemgetter(1))[1]
#            bottomE=max(spectrum_dict.iteritems(), key=operator.itemgetter(1))[1]
#            E=((topE-bottomE)((Lambda-bottomLambda) /(topLambda-bottomLambda)))+bottomE            
#            Intensity=(E-EMin)/(EMax-EMin)*255
       
    return z

def n(Lambda, t=18, p=101325):

    n = (0.0472326 * (173.3 - (1/Lambda)**2)**(-1))+1
    
    return n

def distort( inX, inY, K=[0], Xc=0 , Yc=0):
    '''
    Transforms coordinate values based on Brown's 1966 model.
       
    Parameters
    ----------
    inX : np np.array
        1 x n np.array of X-coordinates to be transformed.
    inY : np np.array
        1 x n np.array of Y-coordinates to be transformed.
    K : np np.array
        1 x n np.array of distortion factors.
        The length of the np.array determines the order of the polynomial. 
    Xc : integer
        X-coord of distortion center.
    Yc : integer
        Y-coord of distortion center.

    Returns
    -------
    outX : np np.array
        1 x n np.array of X-coordinates transformed.
    outY : np np.array
        1 x n np.array of Y-coordinates transformed.
      
    Notes
    -----
    
    '''
    
#    Denom=1
#    index=1
#    r = np.sqrt((inX-np.ones(len(inX))*Xc)**2+(inY-np.ones(len(inY))*Yc)**2)
#    
#    for Kn in K:
#        Denom += Kn * r**(index*2)
#        index += 1
#        
#    outX = inX / Denom + Xc   
#    outY = inY / Denom + Yc

    outX=inX+(inX-Xc)/160
    outY=inY
    
    return outX, outY
 
def findFit(calibrationFile):
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
      
    mainArgs=['4','0',calibrationFile,'0','0']
    #x,y, wavelist are the positions of the peaks in calibrationFile.
    #x,y,waveList,xsig,ysig = readCalibrationData(calibrationFile)
    #fit is the output, which is the ideal p vector.
#    fit = leastsq(main_errors, [279, 90.6, 59, 90, 89, 90, 70.75, 65, 0, 31.64, 210,0], args = (mainArgs), full_output=False, factor=0.05)
#    fit = leastsq(main_errors, [279, 90.6, 59, 90, 89, 90, 70.75, 63.9, 0, 31.62, 210], args = (mainArgs), full_output=False, factor=0.03)#Obviously bad
#    fit = leastsq(main_errors, [279, 90.6, 59, 90, 89, 90, 70.75, 66.1, 0, 31.62, 210], args = (mainArgs), full_output=False, factor=0.03)#Obviously bad
#    fit = leastsq(main_errors, [278.6, 90.6, 61, 90, 91, 90, 70.75, 63.95, 0, 31.9, 210], args = (mainArgs), full_output=True, factor=0.06)
#    fit = leastsq(main_errors, [278.2,90.8,59.0,90.4,88.9,89.7,70.8,65.1,0.758,31.67,203.1], args=mainArgs, full_output=True, diag=[1,1,1,1,1,1,1,1,1,.1,1])

    fit = leastsq(main_errors,[ 271.92998622,   91.03999719,   59.48997316,   89.83000496,   89.37002499,   89.79999531,   68.03002976,   64.9399939,     1.15998754,   31.52736851,  200.00000005], args=mainArgs, full_output=True, diag=[1,1,1,1,1,1,1,1,1,.1,1])
#    fit = leastsq(main_errors, [278.2,90.8,59.0,90.4,88.9,89.7,68.8,65.1,0.758,31.67,203.1], args=mainArgs, full_output=True, diag=[1,1,1,1,1,1,1,1,1,.1,1])

    return fit
   
def medianCombine(inFiles, outFilename):
    
    for k in range(0,len(inFiles)):
        im = pyfits.getdata(inFiles[k])
        if (k == 0):
            cube=np.empty([im.shape[0],im.shape[1],len(inFiles)])
        cube[:,:,k] = im
    med_im=np.median(cube,axis=2)
    pyfits.writeto(outFilename,med_im)
    
def averageCombine(inFiles, outFilename):
    
    for k in range(0,len(inFiles)):
        im = pyfits.getdata(inFiles[k])
        if (k == 0):
            cube=np.empty([im.shape[0],im.shape[1],len(inFiles)])
        cube[:,:,k] = im
    med_im=np.average(cube,axis=2)
    pyfits.writeto(outFilename,med_im)
    
def batchMedian():
    
    start = range(1,721,5)
    stop = range(5,721,5)
    os.chdir('/media/win-desktop/6hrs/Darks')
    
#    mink=10
#    maxk=15
#    files = 720/5
    for j in range(720/5-1):
        inFiles = ['Hg_20s_%03d.FIT' % k for k in range(start[j],start[j+1])]
        outFilename = 'med' + str(start[j]) + '_' + str(start[j+1]) + '.fits' 
        medianCombine(inFiles, outFilename)

def batchDark():
    
    files = range(1,721)

    os.chdir('/media/win-desktop/6hrs/Darks')

#    mink=10
#    maxk=15
#    files = 720/5
    for j in len(files):
        
        inFile = pyfits.getdata(files[j],0)
        inFiles = ['Hg_20s_%03d.FIT' % k for k in range(start[j],start[j+1])]
        outFilename = 'med' + str(start[j]) + '_' + str(start[j+1]) + '.fits' 
        medianCombine(inFile, darkFile, outFilename)

def calibrateImage(inFileName, biasFileName, darkFileName, flatFileName, outFileName, pixelmaskFileName=''):
    
    hdulist = pyfits.open(inFileName)
    inExpTime=hdulist[0].header['EXPOSURE']

    
    inFile = pyfits.getdata(inFileName)  

    
    
    #Bias subtract
    biasFile = pyfits.getdata(biasFileName) 
    outFile=inFile-biasFile
    
    #Darks
    hdulist = pyfits.open(darkFileName)
    darkExpTime = hdulist[0].header['EXPOSURE']
    
    darkFile=pyfits.getdata(darkFileName)
    darkFile=darkFile/darkExpTime*inExpTime
    
    outFile=outFile-darkFile
    
    #Flats
    flatFile=pyfits.getdata(flatFileName)
    
    outFile=outFile/flatFile
    
    #Bad Pixels

    
    #write result
    pyfits.writeto(outFileName, outFile)

def subtractDark(inFile, darkFile, outFilename):
    
    outFile = inFile - darkFile
    pyfits.writeto(outFilename, outFile)





#os.chdir('/media/win-desktop/6hrs/Darks')
#inFiles = ['Hg_6hrs_20s_D%01d.FIT' % k for k in range(1,6)]
#outFilename = 'masterDark.fits' 
#medianCombine(inFiles, outFilename)

#
#p=findFit('../c_noFlat_Hg_0deg_10s.txt')
#print p[0]
#main(p[0])

main()

#batchMedian()

#inFiles = ['../sky_0deg_460_%01d.fits' % k for k in range(1,4)]
#medianCombine(inFiles, '../sky_0deg_460_median.fits')
#inFiles = ['../sky_0deg_1550_%01d.fits' % k for k in range(1,4)]
#medianCombine(inFiles, '../sky_0deg_1550_median.fits')

#inFiles = ['../sky_0deg_460_%01d.fits' % k for k in range(1,4)]
#averageCombine(inFiles, '../sky_0deg_460_average.fits')
#inFiles = ['../sky_0deg_1550_%01d.fits' % k for k in range(1,4)]
#averageCombine(inFiles, '../sky_0deg_1550_average.fits')

#calibrateImage('../sky_0deg_460_median.fits','../mbias.fits','../mdark.fits','../mflat.fits','../c_noFlat_sky_0deg_460_median.fits')
#calibrateImage('../sky_0deg_1550_median.fits','../mbias.fits','../mdark.fits','../mflat.fits','../c_noFlat_sky_0deg_1550_median.fits')

#calibrateImage('../Hg_0deg_10s.fits','../mbias.fits','../mdark.fits','../mflat.fits','../c_Hg_0deg_10s.fits')



#print x(0.502)

#p=main_errors(p = [ 271.92998622,   91.03999719,   59.48997316,   89.83000496,   89.37002499,   89.79999531,   68.03002976,   64.9399939,     1.15998754,   31.52736851,  200.00000005], ['4','1','../c_noFlat_Hg_0deg_10s.txt','1','0'])

#main(p = [279, 90.6, 59, 90, 89, 90, 70.75, 65, 0, 31.64, 210])
#main(p = [ 271.92998622,   91.03999719,   59.48997316,   89.83000496,   89.37002499,   89.79999531,   68.03002976,   64.9399939,     1.15998754,   31.52736851,  200.00000005]) ####best fit so far....
##
##print arrayOut
##
#xi = np.arange(0,9)
#A = np.array([ xi, np.ones(9)])
## linearly generated sequence
#y = [19, 20, 20.5, 21.5, 22, 23, 23, 25.5, 24]
#w = np.linalg.lstsq(A.T,y)[0] # obtaining the parameters
#s=1


#best p for previous fitting:
#p = [ 271.92998622,   91.03999719,   59.48997316,   89.83000496,   89.37002499,   89.79999531,   68.03002976,   64.9399939,     1.15998754,   31.52736851,  200.00000005]
