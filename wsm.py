#Imports
import matplotlib.pyplot as plt     #python/matlab
import random                       #random generator package
import pyfits
import os
import numpy as np
import matplotlib.cm as cm
import wsmtools as wt
from constants import *
import defaults

#least square package
from scipy.optimize.minpack import leastsq

#interpolate function
from scipy import interpolate

#Astro Libraries
from astLib import astSED           

import bisect as bis

import matplotlib.image as mpimg



'''
Initial Parameters-------------------------------------------------------------- 
Wavelength'''
#minLambda=0.390 #min wavelength
#maxLambda=0.795    #max wavelength
#deltaLambda=0.001    #step interval
#maxLambda+=deltaLambda

#Can plot orders from  146 to 73 (about 390 to 795nm). If the wavelength range just above does not cover the orders selected here, this code currently fails!
minOrder=80
maxOrder=100
deltaOrder=1
maxOrder+=deltaOrder
booLog=6 
pixelSize= 5.4 



def doSEDMap(SEDMode=SEDModeFlat, minLambda=0.200, maxLambda=1.000, deltaLambda=0.0001, intNormalize=0):
    '''
    Loads the input Spectrum Energy Density map. It simulates the characteristics of the input beam. 
    
    Parameters
    ----------
    SEDMode : int
        Mode for the creation of the SEDMap
        0=Flat, 1=Random, 2=Sun, 3=from specFile, 4=from Calibration file
        
    minLambda : np.float32
        Lower limit for SEDMap.

    maxLambda : np.float32
        Higher limit for SEDMap.

    deltaLambda : np.float32
        Step between wavelengths.
    
    intNormalize : integer
        If !=0, it normalizes to intNormalize value
        
    Returns
    -------
    SEDMap : np.array
        n x 2 np.array with wavelength, Energy
      
    Notes
    -----
    '''  
    if SEDMode==SEDModeFlat: #Flat
        SEDMap = np.column_stack((np.arange(minLambda, maxLambda, deltaLambda),np.ones(np.arange(minLambda, maxLambda, deltaLambda).size)))

    elif SEDMode==SEDModeRandom: #Random
        SEDMap = np.array([0,0])
        
        for Lambda in range(minLambda, maxLambda + deltaLambda, deltaLambda):
#            Intensity=int()
#            newItem=np.np.array([Lambda,random.random(0.0,1.0)])
            SEDMap = np.vstack((SEDMap,np.array([Lambda,random.random(0.0,1.0)])))
            
        SEDMap = SEDMap[1:,]
                 
    elif SEDMode==SEDModeSolar: #Solar
        
        sol = astSED.SOL        
        
        tempA=sol.wavelength.transpose()*1e-4
        tempB=sol.flux.transpose()            
        SEDMap = np.column_stack((tempA, tempB))    
        
        """Remove rows outside the wavelength range"""
        SEDMap = SEDMap[SEDMap[:,0]>=minLambda]     
        SEDMap = SEDMap[SEDMap[:,0]<=maxLambda]     
                                  
    elif SEDMode==SEDModeFile: #From file
        SEDMap = np.array([0,0])
         
        for line in open(specFile):
            Lambda = float(str(line).split()[0]) #Wavelength
            I = float(str(line).split()[1])       #Intensity
            SEDMap = np.vstack((SEDMap,np.array([Lambda,I])))

        SEDMap=SEDMap[1:,]  
         
    elif SEDMode==SEDModeCalib: #From calibration file
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

def doCCDMap(SEDMap, p = [272.31422902, 90.7157937, 59.6543365, 90.21334551, 89.67646101, 89.82098015, 68.0936684,  65.33694031, 1.19265536, 31.50321471, 199.13548823], SEDMode=0, booPlot=False, specFile='c_noFlat_Hg_0deg_10s.txt', intNormalize=1, booDistort=False, booInterpolate=False, booPlotCalibPoints=False, booPlotLabels=False, plotBackImage='c_noFlat_sky_0deg_460_median.fits',booGaussianFit=False):  
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
    Optics = np.append([[n1,[0,0,0],OpticsPrism,1,0]], [[n2,[0,0,0],OpticsPrism,0,0]], 0)
    
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
    Optics = np.append(Optics, [[s,l,OpticsRGrating,0,d]], 0)
    
    
    #Prism surface 3 (surf #2 on return)
    n4=-n2
    Optics = np.append(Optics,[[n4,[0,0,0],OpticsPrism,1,0]], 0)
    
    #Prism surface 4 (surf #1 on return)
    n5=-n1     
    Optics = np.append(Optics, [[n5,[0,0,0],OpticsPrism,0,0]], 0)

    #Distortion np.array
    K = [] #p[11]
       
    #Launch grid loop. Creates an array of (x,y,lambda)
    CCDX, CCDY, CCDLambda, CCDIntensity, CCDOrder = wp.CCDLoop(SEDMap, Beam , Optics, stheta, fLength) #minLambda ,maxLambda ,deltaLambda ,minOrder ,maxOrder ,deltaOrder ,fLength ,stheta) 
     
    #Distort
    if len(K)!=0: CCDX, CCDY = distort(CCDX, CCDY, K)
        
    return CCDX, CCDY, CCDLambda, CCDIntensity, CCDOrder
    
def doPlot(CCDMap,CalibPoints=False,Labels=False,BackImage='c_noFlat_sky_0deg_460_median.fits'):
        CCDX = CCDMap[CCDMapX] 
        CCDY = CCDMap[CCDMapY] 
        CCDLambda = CCDMap[CCDMapLambda] 
        CCDIntensity= CCDMap[CCDMapIntensity] 
        CCDnOrder= CCDMap[CCDMapOrder] 

        colorTable = np.array((wt.wav2RGB(CCDLambda, CCDIntensity))) 
        
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
        ax1.scatter(CCDX, -CCDY ,s=8, color=colorTable , marker='o', alpha =.5)
#        color=colorTable
#        print random.randrange(-30,-10) random()
#        plt.subplots_adjust(bottom = 0.1)
        if Labels==True:
            for label, x, y in zip(Lambda, CCDX, -CCDZ):
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

def doFindFit(calibDataFileName, p_try=[271.92998622,   91.03999719,   59.48997316,   89.83000496,   89.37002499,   89.79999531,   68.03002976,   64.9399939,     1.15998754,   31.52736851,  200.00000005],factor_try=1,diag_try=[1,1,1,1,1,1,1,1,1,.1,1]):
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

    fit = leastsq(wt.fit_errors,p_try, args=[4,False,calibDataFileName,0,False,False,False,True,'c_noFlat_sky_0deg_460_median.fits',False], full_output=True, factor=factor_try, diag=diag_try)


    return fit
