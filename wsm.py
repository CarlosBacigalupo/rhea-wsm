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



def do_sed_map(SEDMode=SED_MODE_FLAT, minLambda=0.4, maxLambda=0.78, deltaLambda=0.0001, intNormalize=0, specFile=''): 

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
        SEDMap = np.array((np.arange(minLambda, maxLambda + deltaLambda, deltaLambda),np.ones(np.arange(minLambda, maxLambda + deltaLambda, deltaLambda).size)))

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
        
        SEDMap = SEDMap.transpose()
                                  
    elif SEDMode==SED_MODE_FILE: #From flat file
        SEDMap = np.array([])
         
        for line in open(specFile):
            Lambda = float(str(line).split()[0]) #Wavelength
            I = float(str(line).split()[1])       #Intensity
            SEDMap = np.vstack((SEDMap,np.array([Lambda,I])))


    elif SEDMode==SED_MODE_CALIB: #From calibration file
        
        a=np.loadtxt(TEMP_DIR+specFile).transpose()
        SEDMap=np.array((a[2],np.ones(len(a[2]))))
        SEDMap.transpose()

                
    #Normalize the intensity  
    if intNormalize!=0:    
        fluxRange=(max(SEDMap[SEDMapLambda])-min(SEDMap[SEDMapLambda]))
        
        if fluxRange==0:
            SEDMap = np.array((SEDMap[SEDMapLambda], np.ones(SEDMap[SEDMapLambda].size)))  
        else:
            SEDMap = np.array((SEDMap[SEDMapLambda], (SEDMap[SEDMapIntensity]-min(SEDMap[SEDMapIntensity]))/(fluxRange+1) ))
       
    return SEDMap

def do_ccd_map(SEDMap):
    '''
    Computes the projection of n beams of monochromatic light passing through an optical system. 

    Parameters
    ----------
    SEDMap : np.array
        n x 2 np.array [wavelength, Energy]
        
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
       
    #Reads xml file
    Optics, Beams, fLength, p, stheta = xml_parser.read_all()
     
    #hack for RHEA. Needs manual reverse prism on beam return. todo
    Optics[4][0]=-Optics[0][0]
    Optics[3][0]=-Optics[1][0]  
    
    #Launch grid loop. Creates an array of (x,y,lambda, Intensity, Order)
    CCDX, CCDY, CCDLambda, CCDIntensity, CCDOrder = wt.ccd_loop(SEDMap, Beams , Optics, stheta, fLength) #minLambda ,maxLambda ,deltaLambda ,minOrder ,maxOrder ,deltaOrder ,fLength ,stheta) 
     
    #Distort if any distortion data present
    if len(K)!=0: CCDX, CCDY = distort(CCDX, CCDY, K)
        
    return CCDX, CCDY, CCDLambda, CCDIntensity, CCDOrder

def do_find_fit(SEDMap, calibration_data_filename, p_try,factor_try=1 ,diag_try=np.ones(11)):
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

    fit = leastsq(wt.fit_errors, p_try, args=[SEDMap, calibration_data_filename], full_output=True, factor=factor_try, diag=diag_try)

    return fit

def do_read_calibration_file(calibration_image_filename, output_filename,  analyze=True):
    '''
    Extracts peaks with sextractor
    Imports found points into arrays
    Plots found points on image
       
    Parameters
    ----------
    image_filename : string
        Name of the calibration  (fits) 

    analyze : boolean
        Run the analysis (reads the last output otherwise)

    Returns
    -------
    nottin
      
    Notes
    -----
    '''      
    if analyze: ic.analyze_image_sex(calibration_image_filename, output_filename)
    
    #Loads from calibration output file
    image_map_x, image_map_y, image_map_sigx, image_map_sigy = wt.load_image_map_sex(output_filename)
    if image_map_x == []: return
    image_map_x -= 3352/2
    image_map_y -= 2532/2
    
    #Create SEDMap from Mercury emission
    SEDMap = do_sed_map(SEDMode=SED_MODE_CALIB, specFile='Hg_5lines_double.txt')
    
    #Create the model based on default parameters
    CCDMap = do_ccd_map(SEDMap)  
    
    #Create output file with calibration data
    f = open(TEMP_DIR + 'c_' + output_filename,'w')
    for i in range(len(image_map_x)):
        out_string = str(image_map_x[i]) + ' ' + str(image_map_y[i]) + ' 0.0 1 1\n'
        f.write(out_string) 
    f.close()
       
    do_plot_calibration_points(SEDMap, calibration_image_filename, 'c_' + output_filename, CCDMap, labels = False, canvasSize=1)
    
    #Find wavelength of detected points
    CCDX = CCDMap[CCD_MAP_X] 
    CCDY = CCDMap[CCD_MAP_Y] 
    CCDLambda = CCDMap[CCD_MAP_LAMBDA] 
    image_map_lambda = wt.identify_image_map_lambda(SEDMap, CCDX, CCDY, CCDLambda, image_map_x, image_map_y)
    
    #Create output file with calibration data
    f = open(TEMP_DIR + 'c_' + output_filename,'w')
    for i in range(len(image_map_x)):
        out_string = str(image_map_x[i]) + ' ' + str(image_map_y[i]) + ' ' + str(image_map_lambda[i]) + ' ' + str(image_map_sigx[i]) + ' ' + str(image_map_sigy[i]) + '\n'
        f.write(out_string) 
    f.close()
       
    return

def do_full_extract_order(CCDMap, nOrder, image):
    
    CCDX = CCDMap[CCD_MAP_X]
    CCDY = CCDMap[CCD_MAP_Y]
    CCDLambda = CCDMap[CCD_MAP_LAMBDA]
    CCDIntensity = CCDMap[CCD_MAP_INTENSITY]
    CCDOrder=CCDMap[CCD_MAP_ORDER]

    xPlot = CCDX[CCDOrder==nOrder]
    yPlot = CCDY[CCDOrder==nOrder]
    LambdaPlot = CCDLambda[CCDOrder==nOrder]
    
    fLambda = interpolate.interp1d(yPlot, LambdaPlot)
    fX = interpolate.interp1d(yPlot, xPlot, 'quadratic', bounds_error=False)
    
    newX, newY, newLambdas = wt.calculate_from_Y(CCDY, fX, fLambda)
    
    flux = wt.extract_order(newX, newY, image)
    
    return flux, newLambdas

def do_plot_ccd_map(CCDMap, CalibPoints=False, Labels=False, canvasSize=2, backImage='c_noFlat_sky_0deg_460_median.fits'):
        
        CCDX = CCDMap[CCD_MAP_X] 
        CCDY = CCDMap[CCD_MAP_Y] 
        CCDLambda = CCDMap[CCD_MAP_LAMBDA] 
        CCDIntensity = CCDMap[CCD_MAP_INTENSITY] 
        CCDOrder = CCDMap[CCD_MAP_ORDER] 

        colorTable = np.array((wt.wav2RGB(CCDLambda, CCDIntensity))) 
        
        hdulist = pyfits.open(FITS_DIR + backImage)
        imWidth = hdulist[0].header['NAXIS1']
        imHeight = hdulist[0].header['NAXIS2']

        im = pyfits.getdata(FITS_DIR + backImage)
        im[im<0] = 0
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
#        color=colorTable
#        print random.randrange(-30,-10) random()
#        plt.subplots_adjust(bottom = 0.1)
#        plt.plot( x,-z, "o", markersize=7, color=colorTable, markeredgewidth=1,markeredgecolor='g', markerfacecolor='None' )
        
        fig = plt.figure()
        ax1 = fig.add_subplot(111)

        plt.imshow(im,extent=[-imWidth/2 , imWidth/2 , -imHeight/2 , imHeight/2])
        plt.set_cmap(cm.Greys_r)
        ax1.scatter(CCDX, -CCDY ,s=8, color=colorTable , marker='o', alpha =.5)
        plt.axis([-imWidth/2 * canvasSize , imWidth/2 * canvasSize , -imHeight/2 * canvasSize , imHeight/2 * canvasSize])
        plt.title('Order Identification')
        plt.ylabel('pixels')
        plt.xlabel('pixels')
        


        if CalibPoints==True:
            x,y,waveList,xSig,ySig = readCalibrationData(specFile)
            ax1.scatter(x-imWidth/2 , -(y-imHeight/2) ,s=400, color='black', marker='x', alpha=1)
        
        plt.show()

def do_plot_sed_map(SEDMap, point_density=1):
        
        SEDMap = SEDMap[:,::point_density]
        
        colorTable = np.array((wt.wav2RGB(SEDMap[SEDMapLambda], SEDMap[SEDMapIntensity]))) 
        bar_width = (max(SEDMap[SEDMapLambda]) - min(SEDMap[SEDMapLambda])) / SEDMap.size

         
        fig = plt.figure()
        ax1 = fig.add_subplot(111, axisbg='black')
        ax1.bar(SEDMap[SEDMapLambda],SEDMap[SEDMapIntensity], width=bar_width, color=colorTable , edgecolor='none')
        
        plt.ylabel('Intensity')
        plt.xlabel('Wavelength ($\mu$m)')
        plt.axis([min(SEDMap[SEDMapLambda]) ,max(SEDMap[SEDMapLambda])*1.01 ,0 , max(SEDMap[SEDMapIntensity])])

        plt.title('Spectral Energy Distribuition')        
        
        plt.show()

def do_plot_calibration_points(SEDMap, back_image_filename, calibration_data_filename, CCDMap=[], labels = False, canvasSize=2):

        if CCDMap!=[]:
            CCDX = CCDMap[CCD_MAP_X] 
            CCDY = CCDMap[CCD_MAP_Y] 
            CCDLambda = CCDMap[CCD_MAP_LAMBDA] 
            CCDIntensity = CCDMap[CCD_MAP_INTENSITY] 
            CCDOrder = CCDMap[CCD_MAP_ORDER] 

        #Loads from calibration output file
        image_map_x, image_map_y, image_map_lambda , image_map_xsig , image_map_ysig  = wt.read_full_calibration_data(calibration_data_filename)
        if image_map_x==[]: return
        
        #Plot
        hdulist = pyfits.open(FITS_DIR + back_image_filename)
        imWidth = hdulist[0].header['NAXIS1']
        imHeight = hdulist[0].header['NAXIS2']
        im = pyfits.getdata(FITS_DIR + back_image_filename)  
#        plt.ion()
        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        plt.imshow(im,extent=[-imWidth/2 , imWidth/2 , imHeight/2 , -imHeight/2])
        plt.set_cmap(cm.Greys_r)
        
        ax1.scatter(image_map_x, image_map_y ,s=40, color="red" , marker='o', alpha = 0.5, label='Calibration Data')
        if CCDMap!=[]: ax1.scatter(CCDX, CCDY ,s=40, color="blue" , marker='o', alpha = 0.5, label='Model Data')
        
        plt.legend()
        plt.title(str(len(image_map_x))+' point(s) found')
        plt.axis([-imWidth/2 * canvasSize, imWidth/2 * canvasSize, -imHeight/2 * canvasSize, imHeight/2 * canvasSize])

        if (labels==True and CCDMap!=[]):
            full_label = ['('+str(CCDX[x])+', '+str(CCDY[x])+')'+str(CCDLambda[x]) for x in np.arange(len(CCDX))]
            
            for x, y, label in zip(CCDX, CCDY, full_label ):
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
#        plt.draw()
        plt.show()
