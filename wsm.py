#General Imports
import random       #random generator package
import pyfits       #fits package
import os           #os access

#Aliases
import numpy as np
import matplotlib.cm as cm
import bisect as bis
import matplotlib.image as mpimg
import matplotlib.pyplot as plt     #matlab
from scipy.optimize.minpack import leastsq #least square package
from scipy import interpolate #interpolate function
from astLib import astSED #Astro Libraries           

#Custom Packages
from constants import *
import xml_parser
import wsmtools as wt
import image_calibration as ic



def do_SED_map(SEDMode=SED_MODE_FLAT, minLambda=0.200, maxLambda=1., deltaLambda=0.01, intNormalize=0, specFile=''): 
    #todo change format to SEDMap to be equal SEDFlat
    '''
    Creates/Loads the input Spectrum Energy Density map. It simulates the characteristics of the input beam. 
    
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
    if SEDMode==SED_MODE_FLAT: #Flat
        SEDMap = np.array((np.arange(minLambda, maxLambda, deltaLambda),np.ones(np.arange(minLambda, maxLambda, deltaLambda).size)))

    elif SEDMode==SED_MODE_RANDOM: #Random
        
        np.hstack((range(minLambda, maxLambda + deltaLambda, deltaLambda),[random.random() for _ in range(10)]))
        
        
        SEDMap = np.array([0,0])
        for Lambda in range(minLambda, maxLambda + deltaLambda, deltaLambda):
#            Intensity=int()
#            newItem=np.np.array([Lambda,random.random(0.0,1.0)])
            SEDMap = np.vstack((SEDMap,np.array([Lambda,random.random(0.0,1.0)])))
            
        SEDMap = SEDMap[1:,]
                 
    elif SEDMode==SED_MODE_SOLAR: #Solar
        
        sol = astSED.SOL        
        
        tempA=sol.wavelength.transpose()*1e-4
        tempB=sol.flux.transpose()            
        SEDMap = np.column_stack((tempA, tempB))    
        
        #Remove rows outside the wavelength range
        SEDMap = SEDMap[SEDMap[:,0]>=minLambda]     
        SEDMap = SEDMap[SEDMap[:,0]<=maxLambda]     
                                  
    elif SEDMode==SED_MODE_FILE: #From file
        SEDMap = np.array([])
         
        for line in open(specFile):
            Lambda = float(str(line).split()[0]) #Wavelength
            I = float(str(line).split()[1])       #Intensity
            SEDMap = np.vstack((SEDMap,np.array([Lambda,I])))


    elif SEDMode==SED_MODE_CALIB: #From calibration file
        
        a=np.loadtxt(specFile).transpose()
        SEDMap=np.array((a[2],np.ones(len(a[2]))))
        SEDMap.transpose()

                
    #Normalize the intensity"""   
    if intNormalize!=0:    
        fluxRange=(max(SEDMap[SEDMapLambda])-min(SEDMap[SEDMapLambda]))
        
        if fluxRange==0:
            SEDMap = np.array((SEDMap[SEDMapLambda], np.ones(SEDMap[SEDMapLambda].size)))  
        else:
            SEDMap = np.array((SEDMap[SEDMapLambda], (SEDMap[SEDMapIntensity]-min(SEDMap[SEDMapIntensity]))/(fluxRange+1) ))
       
    return SEDMap
    
def doCCDMap(SEDMap, p = [272.31422902, 90.7157937, 59.6543365, 90.21334551, 89.67646101, 89.82098015, 68.0936684,  65.33694031, 1.19265536, 31.50321471, 199.13548823]):
    '''
    Computes the projection of n beams of monochromatic light passing through an optical system. 

    Parameters
    ----------
    SEDMap : np.array
        n x 2 np.array with wavelength, Energy
    
    p : np np.array (optional)
        (beam phi, beam theta, prism1 phi, prism1 theta, prism2 phi, prism2 theta, grating phi, grating theta, grating alpha,
         blaze period (microns), focal length(mm), distortion term) <-- optical arrangement
         
    
    Returns
    -------
    CCDX : np.array
        x coordinate of the target point 
    CCDY : np.array
        x coordinate of the target point 
    CCDLambda : np.array
        wavelength at x,y
    CCDIntensity : np.array
        Intensity at x,y
    CCDOrder : np.array
        Order at x,y

    Notes
    -----  
    '''   
    
    #Distortion np.array
    K = [] #p[11]
       
       
    Optics, Beams, fLength, p, stheta = xml_parser.read_all()
    
    
    #Launch grid loop. Creates an array of (x,y,lambda, Intensity, Order)
    CCDX, CCDY, CCDLambda, CCDIntensity, CCDOrder = wt.CCDLoop(SEDMap, Beams , Optics, stheta, fLength) #minLambda ,maxLambda ,deltaLambda ,minOrder ,maxOrder ,deltaOrder ,fLength ,stheta) 
     
    #Distort if any distortion data present
    if len(K)!=0: CCDX, CCDY = distort(CCDX, CCDY, K)
        
    return CCDX, CCDY, CCDLambda, CCDIntensity, CCDOrder

def doPlot(CCDMap, CalibPoints=False, Labels=False, BackImage='c_noFlat_sky_0deg_460_median.fits'):
        
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

def do_read_calib(output_filename, image_filename='test.fits', analyze=True):
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
    image_map_x, image_map_y = wt.load_image_map()
    
    image_map_lambda_match = image_map_lambda = np.zeros(len(image_map_x))
    
    #Create SEDMap from Mercury emission
    SEDMap = do_SED_map(SEDMode=SED_MODE_CALIB, specFile='Hg_5lines_double.txt')
    
    #Create the model based on default parameters
    CCDX, CCDY, CCDLambda, CCDIntensity, CCDOrder = doCCDMap(SEDMap)  
    
    #Create list of 
    #todo turn this into a 2D array calculation
    for i in range(len(image_map_x)):
    
        distance_array = np.sqrt((CCDX-image_map_x[i])**2+(CCDY-image_map_y[i])**2)
        closest_point=np.min(distance_array)
        closest_point_index=np.where(distance_array==closest_point)[0][0]       
        image_map_lambda[i] = CCDLambda[closest_point_index]
        
        lambda_distance= abs(SEDMap[SEDMapLambda]-image_map_lambda[i])
        lambda_closest_point=np.min(lambda_distance)
        lambda_closest_point_index=np.where(lambda_distance==lambda_closest_point)[0][0]
        image_map_lambda_match[i] = SEDMap[SEDMapLambda][lambda_closest_point_index]
    
    #Create output file with calibration data
    #todo, add sigmas    
    f = open(output_filename,'w')
    for i in range(len(image_map_x)):
        out_string=str(image_map_x[i])+' '+str(image_map_y[i])+' '+str(image_map_lambda_match[i])+' 1 1\n'
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
    ax1.scatter(image_map_x-imWidth/2, -(image_map_y-imHeight/2) ,s=40, color="red" , marker='o', alpha = 0.3)
    plt.title(str(len(image_map_x))+' point(s) found')
    plt.axis([-imWidth/2 , imWidth/2 , -imHeight/2 , imHeight/2])
    plt.show()
    
    return

