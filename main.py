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
from math import cos, sin, acos, asin, pi, atan, degrees, sqrt, radians

import bisect as bis

import matplotlib.image as mpimg



'''
Initial Parameters-------------------------------------------------------------- 
Wavelength'''
minLambda=0.5886 #min wavelength
maxLambda=0.59    #max wavelength
deltaLambda=0.0001    #step interval
maxLambda+=deltaLambda

#Can plot orders from  146 to 73 (about 390 to 795nm)
minOrder=146
maxOrder=73
deltaOrder=-1
maxOrder+=deltaOrder
booLog=6 
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

    return np.hstack([(x_best-x)/xSig,(y_best - y)/ySig]), waveList



def main(p = [272.31422902, 90.7157937, 59.6543365, 90.21334551, 89.67646101, 89.82098015, 68.0936684,  65.33694031, 1.19265536, 31.50321471, 199.13548823], 
         args=['0','1','../c_noFlat_Hg_0deg_10s.txt','1','0','0','0','0','../solar.png']):#'../solar.txt'

    
    '''
    Compute the projection of n beams of monochromatic light passing through an optical system. 

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
    
    global n1, n2, n4, n5, s, l, d, flux, booInterpolate, plotBackImage
    global booPlot, specFile, booPlotCalibPoints, allFlux,booPlotLabels, specFile
    global booGaussianFit
    

    
    #Args breakdown
    SEDMode = int(args[0])
    booPlot = int(args[1])
    specFile = args[2]
    intNormalize = int(args[3])
    booDistort = int(args[4])
    booInterpolate=int(args[5])
    booPlotCalibPoints=int(args[6])
    booPlotLabels=int(args[7])
    plotBackImage=args[8]
    booGaussianFit=int(args[9])
    
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
        #newx = x*cos(pi*p[11]/180) - y*sin(pi*p[11]/180)
        #newy = y*cos(pi*p[11]/180) + x*sin(pi*p[11]/180)
        #CCDMap[:,0]=newx.copy()
        #CCDMap[:,1]=newy.copy()
           
    #Plot
    if booPlot==1:
        doPlot(CCDMap)
    
    return x, y, Lambda

def extractOrder(x,y,image):
    
    #Don't know what these fluxes are yet!
    flux=np.zeros(len(y))
    flux2=np.zeros(len(y))
    flux3=np.zeros(len(y))

    #Grab image and sizes
    hdulist = pyfits.open(image)
    imWidth = hdulist[0].header['NAXIS1']
    imHeight = hdulist[0].header['NAXIS2']
    
    im = pyfits.getdata(image)  

    x=x+imWidth/2
    y=y+imHeight/2

    
    for k in range(0,len(y)):
        ###mikes'
        x_int = round(x[k])
        #if e.g. x = 10.4, this goes from 5 to 15, and xv at the center is 0.4

        
        in_image_temp = im[y[k],x_int-5:x_int+6]

#        in_image_temp=im[y[k],x_int-5]        
#        for i in range(-4,6):
#            in_image_temp = np.hstack((in_image_temp,im[y[k+i-4],x_int+i]))
            
        in_image_temp[in_image_temp < 0] = 0
        xv = np.arange(-5,6) - x[k] + x_int
        flux[k] =  np.sum(in_image_temp * np.exp(-(xv/3.5)**4))
#        flux2[k] = np.sum(in_image_temp)
#        flux3[k] = np.sum(in_image_temp)- np.sum(in_image_temp * np.exp(-(xv/3.5)**4))
        
        #Carlos' trial
#        x_int = int(x[k])
#        res=x[k]-x_int- 0.5
#        sum_in_image_temp = np.sum(im[y[k],x_int-width/2:x_int+width/2])
#        flux[k] = sum_in_image_temp - ((res*im[y[k],x[k]+width/2+1])+ ((1-res)*im[y[k],x[k]-width/2]))
#        flux2[k] = ((res*im[y[k],x[k]+width/2+2])+ ((1-res)*im[y[k],x[k]-width/2]))
#        flux3[k] = sum_in_image_temp 
        
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
    
    dataOut=np.zeros(5)

    
    #Loads SEDMap based on selection. 
    SEDMap = doSEDMap(SEDMode, minLambda, maxLambda, deltaLambda, intNormalize)
    blaze_angle= stheta #Approximately atan(2)
    allFlux=np.array([0])
    allLambdas=np.array([0])
    
    '''Main loop
    Navigates orders within the range given
    For each order navigates the list of wavelenghts in SEDMap constrained by +/- FSP/2'''    
    for nOrder in range(minOrder, maxOrder, deltaOrder):

        #the wavelength range is from the blaze wavelength of the next order and the blaze wavelength of the previous order
        LambdaBlMin = 2*d*sin(blaze_angle)/(nOrder+1) #Was 0.5
        LambdaBlMax = 2*d*sin(blaze_angle)/(nOrder-1) #Was 0.5
 
        #SEDMapLoop is an array of wavelengths (paired with intensity) called SEDMap that are in this particular order.
        SEDMapLoop=SEDMap.copy()
        
        #constrain by +/- FSP (was FSP/2)
        SEDMapLoop = SEDMapLoop[SEDMapLoop[:,0]>=LambdaBlMin]
        if SEDMapLoop.shape[0]>0:       
            SEDMapLoop = SEDMapLoop[SEDMapLoop[:,0]<=LambdaBlMax]     

        #loop lambda for current order
        for Lambda,inI in SEDMapLoop: 
            
            #refractive indeces of prism and air
            nPrism = nkzfs8(Lambda)
            nAir = n(Lambda)
            
            #Computes the unit vector that results from the optical system for a given wavelength and order
            #This is the actual tracing of the ray for each wavelength
            v, isValid = rayTrace(nAir, nPrism, nOrder, Lambda, d, u, n1, n2, n4, n5, s, l)
            
            if isValid: #no errors in calculation and ray makes it back through the prism
                x=v[0]*fLength*1000/pixelSize # x-coord in focal plane in pixels
                z=v[2]*fLength*1000/pixelSize # z-coord in focal plane in pixels
    
                '''Appends data table
                x,z=pojected coordinates
                inI, outI = intensities, inI=from source (file, random,...) 0.0 to 1.0'''
                outI=Intensity(Lambda, minLambda, maxLambda)              
                dataOut= np.vstack((dataOut,np.array([x,z, Lambda, inI*outI ,nOrder]))) 
        
        #Order extraction        
        if (booInterpolate==1 and len(dataOut[dataOut[:,4]==nOrder][:,0])>=3):
            xPlot=dataOut[dataOut[:,4]==nOrder][:,0]
            yPlot=dataOut[dataOut[:,4]==nOrder][:,1]
            LambdaPlot=dataOut[dataOut[:,4]==nOrder][:,2]
            
            fLambda = interpolate.interp1d(yPlot, LambdaPlot)
            fX = interpolate.interp1d(yPlot, xPlot, 'quadratic', bounds_error=False)
          
            
            hdulist = pyfits.open(plotBackImage)
            imWidth = hdulist[0].header['NAXIS1']
            imHeight = hdulist[0].header['NAXIS2']
            
            newY=np.arange(-imHeight/2,imHeight/2)
            
            newX=fX(newY)
            
            nanMap= np.isnan(newX)
            newX=newX[-nanMap]
            newY=newY[-nanMap]
                    
            flux,flux2,flux3 = extractOrder(newX,newY,plotBackImage)
            
            #read flats
            image ='../simple_flat.fits'
            fluxFlat,flux2,flux3 = extractOrder(newX,newY,image)
            
            
            Lambdas=fLambda(newY)
            
            #Blackbody curve to balance flats
            BB=Lambdas**(-4) / (np.exp(14400/Lambdas/3000)- 1)
         
         
            cleanFlux= flux/fluxFlat*BB#/nOrder**2#/np.max(fluxFlat))

            #Write flux to files
#            f = open('../Order_'+ str(nOrder) +'.txt', 'w')
#            for k in range(0,len(Lambdas)):
#                outText=str(Lambdas[k])+'\t'+ str(cleanFlux[k])+'\n'
#                f.write(outText)
#            f.close()
#            
#            peakIndex=np.where(normalizedFlux==np.max(normalizedFlux))[0][0]
#            peakIndex=len(normalizedFlux)/2
#            Range=(len(normalizedFlux)- peakIndex)*0.6
            

            # Fit a Gaussian             
            if booGaussianFit==1: 
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
                    
                    a,FWHMIndex = find_nearest(Y,np.max(Y)/2)
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
                    currDistance=sqrt((allFlux[intersectStart+k]-cleanFlux[k])**2)
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
             
#            if booInterpolate==1:   
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

    if booInterpolate==1:   
#        a=np.sort(allFlux,0)
#        params = {'legend.fontsize': 200,
#          'legend.linewidth': 2}
#        plt.rcParams.update(params)
        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        ax1.plot(allLambdas,allFlux)
#        ax1.scatter(allLambdas,allFlux ,s=0.1, color='black' , marker='o', alpha =1)
        plt.title('Sodium Doublet')
        plt.ylabel('Intensity (Relative Units)')
        plt.xlabel('Wavelength (Micrometers)')
        plt.show() 
       
    CCDMap=dataOut[1:,]



    return CCDMap

#backImage='../sky_0deg_2.FIT'

def gauss(x, p): # p[0]==mean, p[1]==stdev
    return 1.0/(p[1]*np.sqrt(2*np.pi))*np.exp(-(x-p[0])**2/(2*p[1]**2))

def find_nearest(array,value):
    idx=(np.abs(array-value)).argmin()
    return array[idx], idx

def doPlot(CCDMap):
        x = CCDMap[:,0] 
        z = CCDMap[:,1] 
        Lambda = CCDMap[:,2] 
        Intensity= CCDMap[:,3] 

        colorTable = np.array((wav2RGB(Lambda, Intensity))) 
        
        hdulist = pyfits.open(plotBackImage)
        imWidth = hdulist[0].header['NAXIS1']
        imHeight = hdulist[0].header['NAXIS2']

        im = pyfits.getdata(plotBackImage)
        im[im<0]=0
        im /= im.max()
        im = np.sqrt(im) #Remove this line for Hg
#        im = np.sqrt(im) #Remove this line for Hg
#    

#        labels = np.array([0])
#         
#        for line in open('../solar2.txt'):
#            Lambda = float(str(line).split()[0]) #Wavelength
#            labels = np.vstack((labels,np.array([Lambda])))
#
#        labels=labels[1:]
        
#        im=mpimg.imread('../solar.png')
        
        fig = plt.figure()
        ax1 = fig.add_subplot(111)

        plt.imshow(im,extent=[-imWidth/2 , imWidth/2 , -imHeight/2 , imHeight/2])
        plt.set_cmap(cm.Greys_r)
        ax1.scatter(x, -z ,s=8, color=colorTable , marker='o', alpha =.5)
#        color=colorTable
#        print random.randrange(-30,-10) random()
#        plt.subplots_adjust(bottom = 0.1)
        if booPlotLabels==1:
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
        
        
        
        if booPlotCalibPoints==1:
            x,y,waveList,xSig,ySig = readCalibrationData(specFile)
            ax1.scatter(x-imWidth/2 , -(y-imHeight/2) ,s=400, color='black', marker='x', alpha=1)



#        plt.plot( x,-z, "o", markersize=7, color=colorTable, markeredgewidth=1,markeredgecolor='g', markerfacecolor='None' )
        
        plt.title('Order Identification')

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
    z=1

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
    mainArgs=['4','0',calibrationFile,'0','0','0','0','1','../c_noFlat_sky_0deg_460_median.fits']   
    #old mainArgs=['4','0',calibrationFile,'0','0']
    #x,y, wavelist are the positions of the peaks in calibrationFile.
    #x,y,waveList,xsig,ysig = readCalibrationData(calibrationFile)
    #fit is the output, which is the ideal p vector.

    fit = leastsq(main_errors,p_try, args=mainArgs, full_output=True, factor=factor_try, diag=diag_try)
#    fit = leastsq(main_errors, [278.2,90.8,59.0,90.4,88.9,89.7,68.8,65.1,0.758,31.67,203.1], args=mainArgs, full_output=True,factor,diag )

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
    
    #outFile=outFile/flatFile
    
    #Bad Pixels

    
    #write result
    pyfits.writeto(outFileName, outFile)

def subtractDark(inFileName, darkFileName, outFilename):
    
    inFile = pyfits.getdata(inFileName) 
    darkFile = pyfits.getdata(darkFileName) 
    
    outFile = inFile - darkFile
    
    pyfits.writeto(outFilename, outFile)
                              
                                
def compareTemp(inFileName1,inFileName2, outFilename):


    im1 = pyfits.getdata(inFileName1)
    im2 = pyfits.getdata(inFileName2)

    outFile=im2-im1

    pyfits.writeto(outFilename, outFile)


   
   
   
   
   
   
   
   
   
   
   







#os.chdir('/media/win-desktop/6hrs/Darks')
#inFiles = ['Hg_6hrs_20s_D%01d.FIT' % k for k in range(1,6)]
#outFilename = 'masterDark.fits' 
#medianCombine(inFiles, outFilename)

#
#p=findFit('../c_noFlat_Hg_0deg_10s.txt')
#print p[0]
#main(p[0])

#main()

#compareTemp('../6hrs/Medians/med1_5.fits','../6hrs/Medians/med711_715.fits','../6hrs/Results/out1.fits')

#batchAverage2()

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

#calibrateImage('../sky_0deg_460_1.fits','../mbias.fits','../mdark.fits','../mflat.fits','../c_noFlat_sky_0deg_460_1.fits')
#calibrateImage('../sky_0deg_460_2.fits','../mbias.fits','../mdark.fits','../mflat.fits','../c_noFlat_sky_0deg_460_2.fits')
#calibrateImage('../sky_0deg_460_3.fits','../mbias.fits','../mdark.fits','../mflat.fits','../c_noFlat_sky_0deg_460_3.fits')

#calibrateImage('../sky_0deg_1550_median.fits','../mbias.fits','../mdark.fits','../mflat.fits','../c_noFlat_sky_0deg_1550_median.fits')




#inFiles = ['../dark_a%01d.fits' % k for k in range(1,4)]
#medianCombine(inFiles, '../dark_a_median.fits')
#
#inFiles = ['../dark_a%01d.fits' % k for k in range(1,4)]
#averageCombine(inFiles,  '../dark_a_average.fits')
#
#
#subtractDark('../arcturus_1min_1.fits','../dark_a_average.fits','../c_arcturus_1min_1.fits')
#subtractDark('../arcturus_1min_2.fits','../dark_a_average.fits','../c_arcturus_1min_2.fits')
#subtractDark('../arcturus_1min_3.fits','../dark_a_average.fits','../c_arcturus_1min_3.fits')
#subtractDark('../arcturus_1min_4.fits','../dark_a_average.f#inFiles = ['../dark_a%01d.fits' % k for k in range(1,4)]
#medianCombine(inFiles, '../dark_a_median.fits')
#
#inFiles = ['../dark_a%01d.fits' % k for k in range(1,4)]
#averageCombine(inFiles,  '../dark_a_average.fits')
#
#
#subtractDark('../arcturus_1min_1.fits','../dark_a_average.fits','../c_arcturus_1min_1.fits')
#subtractDark('../arcturus_1min_2.fits','../dark_a_average.fits','../c_arcturus_1min_2.fits')
#subtractDark('../arcturus_1min_3.fits','../dark_a_average.fits','../c_arcturus_1min_3.fits')
#subtractDark('../arcturus_1min_4.fits','../dark_a_average.fits','../c_arcturus_1min_4.fits')
#
#inFiles = ['../c_arcturus_1min_%01d.fits' % k for k in range(1,5)]
#medianCombine(inFiles, '../c_arcturus_1min_median.fits')
#
#inFiles = ['../c_arcturus_1min_%01d.fits' % k for k in range(1,5)]
#averageCombine(inFiles, '../c_arcturus_1min_average.fits')its','../c_arcturus_1min_4.fits')
#
#inFiles = ['../c_arcturus_1min_%01d.fits' % k for k in range(1,5)]
#medianCombine(inFiles, '../c_arcturus_1min_median.fits')
#
#inFiles = ['../c_arcturus_1min_%01d.fits' % k for k in range(1,5)]
#averageCombine(inFiles, '../c_arcturus_1min_average.fits')


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
#'../c_noFlat_Hg_0deg_10s.fits')#,'../solar.png')#,'../mdark.fits')#,'../simple_flat.fits')#


#Fraunhofer Lines-CCDMap with labels***********************************************************************************
#identify fraunhofer
#SED fraunhofer
#(SEDMode(0=Max, 1=Random, 2=Sun, 3=from specFile, 4=from CalibFile), Plot?, specFile/calibfile, 
#Normalize intensity? (0=no, #=range), Distort?, Interpolate, PlotCalibPoints, booPlotLabels, plotbackfile, 
#Gaussian) <-- other options#main(args=['4','1','../solar.txt','1','0','0','0','1','../c_noFlat_sky_0deg_460_median.fits'])
#p = [272.31422902, 90.7157937, 59.6543365, 90.21334551, 89.67646101, 89.82098015, 68.0936684,  65.33694031, 1.19265536, 31.50321471, 199.13548823]
#args=['4','1','../solar.txt','1','0','0','0','1','../c_noFlat_sky_0deg_460_median.fits','0']
#main(p,args)
#********************************************************************************************************


#Arcturus callibration**********************************************************************************
#p_try=[ 272.45928478  , 91.17527895 ,  59.25339245 ,  89.61147631 ,  89.63791054,   89.81178189 ,  68.1669475 ,   63.3555271 ,   1.10221798 ,  31.8848023,  199.70165672] #best
##p_try=[250.6, 64.2, 58.3, 77.9, 89.5, 85.5, 58.5, 61.7, 1.6, 33.7, 193.5] #doeesn't work
## ok p_try=[279, 90.6, 59, 90, 89, 90, 70.75, 65, 0, 31.64, 210]
##p_try = [278.6, 90.6, 61, 90, 91, 90, 70.75, 63.95, 0, 31.9, 210]
##p_try = [278.2,90.8,59.0,90.4,88.9,89.7,70.8,65.1,0.758,31.67,203.1]
#only h-alpha
#p_try = [ 272.46778143,   91.14782003 ,  59.2663218 ,   89.65582564 ,  89.62544383,   89.76937348  , 68.19750843 ,  63.29873297  ,  1.10998689  , 31.92112407,  199.70153725]
#
##next
#p_try = [ 272.45928478,   91.17527895  , 59.25339245 ,  89.61147631 ,  89.63791054,  89.81178189 ,  68.1669475 ,   63.3555271,     1.10221798 ,  31.8848023,  199.70165672]
#
#diag_try=[1,1,1,1,1,1,1,1,1,0.1,1]
#factor_try=0.3
#p=findFit('../arcturus.txt',p_try, factor_try, diag_try)
#print p[0]
#main(p[0],args=['4','1','../arcturus.txt','1','0','0','1','0','../c_arcturus_1min_4.fits','0'])
#********************************************************************************************************


#Arcturus plot********************************************************************************************************
#(beam phi, beam theta, prism1 phi, prism1 theta, prism2 phi, prism2 theta, grating phi, grating theta, grating alpha,
#         blaze period (microns), focal length(mm), distortion term)
#(SEDMode(0=Max, 1=Random, 2=Sun, 3=from specFile, 4=from CalibFile), Plot?, specFile/calibfile, 
#Normalize intensity? (0=no, #=range), Distort?, Interpolate, PlotCalibPoints, booPlotLabels, plotbackfile, 
#Gaussian) <-- other options
#p=[ 272.45928478  , 91.5 ,  59.25339245 ,  89.61147631 ,  89.63791054,   89.81178189 ,  68.1669475 ,   63.3555271 ,   1.10221798 ,  31.8848023,  200] 
#p=[ 272.35233865 ,  91.17544429  , 59.24907075 ,  90.19912466  , 89.63148374,   89.2621423  ,  68.1372743   , 62.95098288   , 1.0992334  ,  32.24603515,  199.65078345]
#p = [ 271.9473,   91.137 ,  59.2663218 ,   89.65582564 ,  89.62544383,   89.76937348  , 68.19750843 ,  63.29873297  ,  1.65  , 31.92112407,  199.70153725]
#main(p,args=['0','1','../arcturus.txt','1','0','1','0','0','../c_arcturus_1min_4.fits','0'])
#*********************************************************************************************************************

#Spectral Resolution-CCDMap mercury with labels******************************************************************************************************************
#(beam phi, beam theta, prism1 phi, prism1 theta, prism2 phi, prism2 theta, grating phi, grating theta, grating alpha,
#         blaze period (microns), focal length(mm))
#(SEDMode(0=Max, 1=Random, 2=Sun, 3=from specFile, 4=from CalibFile), Plot?, specFile/calibfile, 
#Normalize intensity? (0=no, #=range), Distort?, Interpolate, PlotCalibPoints, booPlotLabels, plotbackfile,
# Gaussian) <-- other options
#p = [272.31422902, 90.7157937, 59.6543365, 90.21334551, 89.67646101, 89.82098015, 68.0936684,  65.33694031, 1.19265536, 31.50321471, 199.13548823]
#args=['4','1','../c_noFlat_Hg_0deg_10s.txt','1','0','0','0','1','../c_noFlat_Hg_0deg_10s.fits','0']
#main(p,args)
#*********************************************************************************************************************


#Spectral Resolution-Gaussian Fit in Orders******************************************************************************************************************
#(beam phi, beam theta, prism1 phi, prism1 theta, prism2 phi, prism2 theta, grating phi, grating theta, grating alpha,
#         blaze period (microns), focal length(mm), distortion term)
#(SEDMode(0=Max, 1=Random, 2=Sun, 3=from specFile, 4=from CalibFile), Plot?, specFile/calibfile, 
#Normalize intensity? (0=no, #=range), Distort?, Interpolate, PlotCalibPoints, booPlotLabels, plotbackfile,
# Gaussian) <-- other options
#p = [272.31422902, 90.7157937, 59.6543365, 90.21334551, 89.67646101, 89.82098015, 68.0936684,  65.33694031, 1.19265536, 31.50321471, 199.13548823]
#args=['0','0','../c_noFlat_Hg_0deg_10s.txt','1','0','1','0','0','../c_noFlat_Hg_0deg_10s.fits','1']
#main(p,args)
#*********************************************************************************************************************


#Complete Solar Spectrum***********************************************************************************
#(SEDMode(0=Max, 1=Random, 2=Sun, 3=from specFile, 4=from CalibFile), Plot?, specFile/calibfile, 
#Normalize intensity? (0=no, #=range), Distort?, Interpolate, PlotCalibPoints, booPlotLabels, plotbackfile, 
#Gaussian) <-- other options#main(args=['4','1','../solar.txt','1','0','0','0','1','../c_noFlat_sky_0deg_460_median.fits'])
p = [272.31422902, 90.7157937, 59.6543365, 90.21334551, 89.67646101, 89.82098015, 68.0936684,  65.33694031, 1.19265536, 31.50321471, 199.13548823]
args=['0','0','../solar.txt','1','0','1','0','0','../c_noFlat_sky_0deg_460_median.fits','0']
main(p,args)
#********************************************************************************************************


#Spectral Resolution-CCDMap mercury without labels and with calibpoints******************************************************************************************************************
#(beam phi, beam theta, prism1 phi, prism1 theta, prism2 phi, prism2 theta, grating phi, grating theta, grating alpha,
#         blaze period (microns), focal length(mm))
#(SEDMode(0=Max, 1=Random, 2=Sun, 3=from specFile, 4=from CalibFile), Plot?, specFile/calibfile, 
#Normalize intensity? (0=no, #=range), Distort?, Interpolate, PlotCalibPoints, booPlotLabels, plotbackfile,
# Gaussian) <-- other options
#p = [272.31422902, 90.7157937, 59.6543365, 90.21334551, 89.67646101, 89.82098015, 68.0936684,  65.33694031, 1.19265536, 31.50321471, 199.13548823]
#args=['4','1','../c_noFlat_Hg_0deg_10s.txt','1','0','0','1','0','../c_noFlat_Hg_0deg_10s.fits','0']
#main(p,args)
#*********************************************************************************************************************

#Plot Errors Example**********************************************************************************
#p_try=[ 272.45928478  , 91.17527895 ,  59.25339245 ,  89.61147631 ,  89.63791054,   89.81178189 ,  68.1669475 ,   63.3555271 ,   1.10221798 ,  31.8848023,  199.70165672] #best
##p_try=[250.6, 64.2, 58.3, 77.9, 89.5, 85.5, 58.5, 61.7, 1.6, 33.7, 193.5] #doeesn't work
## ok p_try=[279, 90.6, 59, 90, 89, 90, 70.75, 65, 0, 31.64, 210]
##p_try = [278.6, 90.6, 61, 90, 91, 90, 70.75, 63.95, 0, 31.9, 210]
##p_try = [278.2,90.8,59.0,90.4,88.9,89.7,70.8,65.1,0.758,31.67,203.1]
#only h-alpha
#p_try = [ 272.46778143,   91.14782003 ,  59.2663218 ,   89.65582564 ,  89.62544383,   89.76937348  , 68.19750843 ,  63.29873297  ,  1.10998689  , 31.92112407,  199.70153725]
#
###next
#p = [272.31422902, 90.7157937, 59.6543365, 90.21334551, 89.67646101, 89.82098015, 68.0936684,  65.33694031, 1.19265536, 31.50321471, 199.13548823]
#mainArgs=['4','0','../c_noFlat_Hg_0deg_10s.txt','1','0','0','1','0','../c_noFlat_Hg_0deg_10s.fits','0']
#temp , waves = main_errors(p, mainArgs)
#x2=temp[0]
#y2=temp[1]
#
#d=np.sqrt(x2**2+y2**2)
#plt.ylabel('Error (Pixels)')
#plt.xlabel('Emission Line (micrometers)')
#plt.title('Fitting Error')
#plt.scatter(waves,d)
#plt.show()
#********************************************************************************************************

#Solar SPectrum CCDMap without labels***********************************************************************************
#identify fraunhofer
#SED fraunhofer
#(SEDMode(0=Max, 1=Random, 2=Sun, 3=from specFile, 4=from CalibFile), Plot?, specFile/calibfile, 
#Normalize intensity? (0=no, #=range), Distort?, Interpolate, PlotCalibPoints, booPlotLabels, plotbackfile, 
#Gaussian) <-- other options#main(args=['4','1','../solar.txt','1','0','0','0','1','../c_noFlat_sky_0deg_460_median.fits'])
#p = [272.31422902, 90.7157937, 59.6543365, 90.21334551, 89.67646101, 89.82098015, 68.0936684,  65.33694031, 1.19265536, 31.50321471, 199.13548823]
#args=['0','1','../solar.txt','1','0','0','0','0','../c_noFlat_sky_0deg_460_median.fits','0']
#main(p,args)
#********************************************************************************************************


#Hg full spectrum******************************************************************************************************************
#(beam phi, beam theta, prism1 phi, prism1 theta, prism2 phi, prism2 theta, grating phi, grating theta, grating alpha,
#         blaze period (microns), focal length(mm), distortion term)
#(SEDMode(0=Max, 1=Random, 2=Sun, 3=from specFile, 4=from CalibFile), Plot?, specFile/calibfile, 
#Normalize intensity? (0=no, #=range), Distort?, Interpolate, PlotCalibPoints, booPlotLabels, plotbackfile,
# Gaussian) <-- other options
#p = [272.31422902, 90.7157937, 59.6543365, 90.21334551, 89.67646101, 89.82098015, 68.0936684,  65.33694031, 1.19265536, 31.50321471, 199.13548823]
#args=['0','0','../c_noFlat_Hg_0deg_10s.txt','1','0','1','0','0','../c_noFlat_Hg_0deg_10s.fits','0']
#main(p,args)
#*********************************************************************************************************************

