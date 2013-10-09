#external Packages
import random        
import pyfits       
import os           
import numpy as np
import bisect as bis
import matplotlib.cm as cm
import matplotlib.image as mpimg
import matplotlib.pyplot as plt    
from scipy.optimize.minpack import leastsq
from scipy.optimize import * #least square package
from scipy import interpolate #interpolate function
from astLib import astSED #Astro Libraries           

#internal modules
from wsmtools.constants import *
import wsmtools as wt

def do_sed_map(SEDMode=SED_MODE_FLAT, minLambda=0.4, maxLambda=0.78, deltaLambda=0.0001, intNormalize=0, spectrumFileName=''): 
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
        
    specFile : str
        File with the SED of the source (used when SEDmode = (SED_MODE_CALIB or SED_MODE_FILE))
        
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
        SEDMap = np.array([tempA, tempB])    
        SEDMap = SEDMap.transpose()
                                  
    elif SEDMode==SED_MODE_FILE: #From flat file
        SEDMap = np.loadtxt(APP_PATH + SPECTRUM_PATH + spectrumFileName, unpack=True)
        
    elif SEDMode==SED_MODE_CALIB: #From calibration file
        a=np.loadtxt(TEMP_PATH + spectrumFileName, unpack=True)
        SEDMap=np.array((a[2],np.ones(len(a[2]))))
#        SEDMap.transpose()
                
                
    #Remove rows outside the wavelength range
    SEDMap = np.array([SEDMap[0][SEDMap[0]>=minLambda],SEDMap[1][SEDMap[0]>=minLambda]])     
    SEDMap = np.array([SEDMap[0][SEDMap[0]<=maxLambda],SEDMap[1][SEDMap[0]<=maxLambda]])     
                
                
    #Normalize the intensity  
    if intNormalize!=0:    
        fluxRange=(max(SEDMap[SEDMapLambda])-min(SEDMap[SEDMapLambda]))
        
        if fluxRange==0:
            SEDMap = np.array((SEDMap[SEDMapLambda], np.ones(SEDMap[SEDMapLambda].size)))  
        else:
            SEDMap = np.array((SEDMap[SEDMapLambda], (SEDMap[SEDMapIntensity]-min(SEDMap[SEDMapIntensity]))/(fluxRange+1) ))
       
    return SEDMap

def do_ccd_map(SEDMap ,specXMLFileName, activeCamera=0, p_try = []):
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
            
    #Reads xml file
    Beams, Optics, Cameras, p, stheta = wt.xml.read_all(specXMLFileName, p_try)   
    
    #hack for RHEA. Needs manual reverse of prism on beam return. todo
    if specXMLFileName[-8:]=='rhea.xml':
        Optics[4][0]=-Optics[0][0]
        Optics[3][0]=-Optics[1][0]  
    
    #Launch grid loop. Creates an array of (x,y,lambda, Intensity, Order)
    CCDX, CCDY, CCDLambda, CCDIntensity, CCDOrder, CCDBeamID = wt.optics.ccd_loop(SEDMap, Beams , Optics, Cameras[activeCamera], stheta)
     
    #Distort if any distortion data present
    K = [p[11], p[12], p[13]]
    Xc = p[14]
    Yc = p[15]
    CCDX, CCDY = wt.distort(CCDX, CCDY, K, Xc, Yc)
    
    return CCDX, CCDY, CCDLambda, CCDIntensity, CCDOrder, CCDBeamID

def do_find_fit(SEDMap, specXMLFileName, calibrationDataFileName, activeCameraIndex, p_try = [], factorTry=1 ,diagTry = [], showStats = False, maxfev = 1000, booWriteP = True):
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
    
    if p_try==[]: p_try = wt.xml.read_p(specXMLFileName)
    if diagTry==[]: diagTry = np.ones(len(p_try))
    
    while True:
        fit = leastsq(wt.fit_errors, p_try, args=[SEDMap, specXMLFileName, calibrationDataFileName, activeCameraIndex], full_output=True, factor = factorTry, diag = diagTry, maxfev = maxfev)
        if fit[-1] != 0: #workaround to avoid inconsistent message 'wrong input parameters(error code 0)'
            break
        
    if showStats: wt.fitting_stats(fit)
    
    if booWriteP: wt.xml.write_p(fit[0], specXMLFileName)
    
    return fit

def do_read_calibration_file(calibrationImageFileName, specXMLFileName, outputFileName,  booAnalyse=True, booAvgAdjust = True, booPlotInitialPoints = False, booPlotFinalPoints = False):
    '''
    Extracts peaks with sextractor
    Imports found points into arrays
    Plots found points on image
    Assigns initial wavelength based on proximity
    Asks for user input to refine assignment
       
    Parameters
    ----------
    calibrationImageFileName : string
        Name of the calibration file (fits) 

    specFileName : string
        Run the analysis (reads the last output otherwise)
        
    booAnalyse : boolean
    
    booPlotInitialPoints : boolean



    Returns
    -------
    nottin
      
    Notes
    -----
    '''      
    global Beams, Optics, Cameras, a, b 
    
    if booAnalyse: wt.ia.analyse_image_sex(calibrationImageFileName, outputFileName)
    
    #Loads coordinate points from calibration output file
    imageMapX, imageMapY, image_map_sigx, image_map_sigy = wt.ia.load_image_map_sex(outputFileName)
    if imageMapX == []: return
    imageMapX -= int(Cameras[0][CamerasWidth])/2
    imageMapY -= int(Cameras[0][CamerasHeight])/2
    
    #Create SEDMap from Mercury emission
    SEDMap = do_sed_map(SEDMode=SED_MODE_FILE, spectrumFileName='hg_spectrum.txt')
    
    #Create the model based on default parameters
    CCDMap = do_ccd_map(SEDMap, specXMLFileName)  
    
    #Create initial output file from calibration data
    f = open(TEMP_PATH + 'c_' + outputFileName,'w')
    for i in range(len(imageMapX)):
        out_string = str(imageMapX[i]) + ' ' + str(imageMapY[i]) + ' 0.0 1 1\n'
        f.write(out_string) 
    f.close()
    
    if booPlotInitialPoints:
        do_plot_calibration_points(calibrationImageFileName, 'c_' + outputFileName, CCDMap, booLabels = False, canvasSize=1.3, title = 'Calibration vs Model before wavelength matching ')
    
    #Find wavelength of detected points (first approximation)
    CCDX = CCDMap[CCD_MAP_X] 
    CCDY = CCDMap[CCD_MAP_Y] 
    CCDLambda = CCDMap[CCD_MAP_LAMBDA] 
    imageMapLambda = wt.ia.identify_imageMapLambda_auto(SEDMap, CCDX, CCDY, CCDLambda, imageMapX, imageMapY, booAvgAdjust = booAvgAdjust)
    
    #Create temporary output file with all calibration data
    f = open(TEMP_PATH + 'c_temp_' + outputFileName,'w')
    for i in range(len(imageMapX)):
        out_string = str(imageMapX[i]) + ' ' + str(imageMapY[i]) + ' ' + str(imageMapLambda[i]) + ' ' + str(image_map_sigx[i]) + ' ' + str(image_map_sigy[i]) + '\n'
        f.write(out_string) 
    f.close()

    #Correct/confirm found wavelengths
#    imageMapLambda = wt.ia.identify_imageMapLambda_manual(SEDMap, CCDX, CCDY, CCDLambda, imageMapX, imageMapY, imageMapLambda, calibrationImageFileName, Cameras)

    #Create final output file with all calibration data
    f = open(TEMP_PATH + 'c_' + outputFileName,'w')
    for i in range(len(imageMapX)):
        out_string = str(imageMapX[i]) + ' ' + str(imageMapY[i]) + ' ' + str(imageMapLambda[i]) + ' ' + str(image_map_sigx[i]) + ' ' + str(image_map_sigy[i]) + '\n'
        f.write(out_string) 
    f.close()   
    
    #Plot detected points with assigned wavelengths
    if booPlotFinalPoints:
        do_plot_calibration_points(calibrationImageFileName, 'c_' + outputFileName, CCDMap, booLabels = True, canvasSize=1, title = 'Calibration vs Model (wavelength assigned)')
       
    return

def do_extract_order(CCDMap, nOrder, image):
    
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
    
    flux = wt.ia.extract_order(newX, newY, image)
    
    return flux, newLambdas

def do_plot_ccd_map(CCDMap, canvasSize=1, backImage=''):
        
        CCDX = CCDMap[CCD_MAP_X] 
        CCDY = CCDMap[CCD_MAP_Y] 
        CCDLambda = CCDMap[CCD_MAP_LAMBDA] 
        CCDIntensity = CCDMap[CCD_MAP_INTENSITY] 
        CCDOrder = CCDMap[CCD_MAP_ORDER] 
        CCDBeamID = CCDMap[CCD_MAP_BEAMID] 

        colorTable = np.array((wt.wav2RGB(CCDLambda, CCDIntensity))) 
        
        imWidth = int(Cameras[0][CamerasWidth])
        imHeight = int(Cameras[0][CamerasHeight])
        
        fig = plt.figure()
        ax1 = fig.add_subplot(111, axisbg = 'black')

        if backImage!='':
            im = pyfits.getdata(FITS_PATH + backImage)
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
            plt.imshow(im,extent=[-imWidth/2 , imWidth/2 , -imHeight/2 , imHeight/2])
            plt.set_cmap(cm.Greys_r)
        


        ax1.scatter(CCDX, -CCDY ,s=8, color=colorTable , marker='o', alpha =.5)
        plt.axis([-imWidth/2 * canvasSize , imWidth/2 * canvasSize , -imHeight/2 * canvasSize , imHeight/2 * canvasSize])
        plt.title('Order Identification')
        plt.ylabel('pixels')
        plt.xlabel('pixels')
        
#        if CalibPoints==True:
#            x,y,waveList,xSig,ySig = readCalibrationData(specFile)
#            ax1.scatter(x-imWidth/2 , -(y-imHeight/2) ,s=400, color='black', marker='x', alpha=1)
        
        plt.show()

def do_plot_sed_map(SEDMap, title='', point_density=1):
        
        SEDMap = SEDMap[:,::point_density]
        
        colorTable = np.array((wt.wav2RGB(SEDMap[SEDMapLambda], SEDMap[SEDMapIntensity]))) 
        bar_width = (max(SEDMap[SEDMapLambda]) - min(SEDMap[SEDMapLambda])) / 100

         
        fig = plt.figure()
        ax1 = fig.add_subplot(111, axisbg='black')
        ax1.bar(SEDMap[SEDMapLambda],SEDMap[SEDMapIntensity], width=bar_width, color=colorTable , edgecolor='none')
        
        plt.ylabel('Intensity')
        plt.xlabel('Wavelength ($\mu$m)')
        plt.axis([min(SEDMap[SEDMapLambda]) ,max(SEDMap[SEDMapLambda])*1.01 ,0 , max(SEDMap[SEDMapIntensity])])

        if title=='':title = 'Spectral Energy Distribuition'
        plt.title(title)        
        
        plt.show()

def do_plot_calibration_points(backImageFileName, calibrationDataFileName, CCDMap=[], booLabels = False, canvasSize=2, title=''):

        if CCDMap!=[]:
            CCDX = CCDMap[CCD_MAP_X] 
            CCDY = CCDMap[CCD_MAP_Y] 
            CCDLambda = CCDMap[CCD_MAP_LAMBDA] 
            CCDIntensity = CCDMap[CCD_MAP_INTENSITY] 
            CCDOrder = CCDMap[CCD_MAP_ORDER] 

        #Loads from calibration output file
        imageMapX, imageMapY, imageMapLambda , imageMapXsig , imageMapYsig  = wt.ia.read_full_calibration_data(calibrationDataFileName)
        if imageMapX==[]: return
        
        #Plot
        hdulist = pyfits.open(FITS_PATH + backImageFileName)
        imWidth = hdulist[0].header['NAXIS1']
        imHeight = hdulist[0].header['NAXIS2']
        im = pyfits.getdata(FITS_PATH + backImageFileName) 
        imNorm = wt.ic.normalise_image(im) 
#        plt.ion()
        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        plt.imshow(imNorm,extent=[-imWidth/2 , imWidth/2 , imHeight/2 , -imHeight/2])
        plt.set_cmap(cm.Greys_r)
        
        ax1.scatter(imageMapX, imageMapY ,s=50, color="red" , marker='o', alpha = 0.5, label='Calibration Data')
        if booLabels:
            fullLabel = [str(imageMapLambda[x]) for x in np.arange(len(imageMapLambda))]
            
            for x, y, label in zip(imageMapX, imageMapY, fullLabel ):
                plt.annotate(
                    label, 
                    xy = (x, y), xytext = (0,-25),
                    textcoords = 'offset points', ha = 'right', va = 'bottom',
                    bbox = dict(boxstyle = 'round,pad=0.5', fc = 'white', alpha = 1))
#                    arrowprops = dict(arrowstyle="wedge,tail_width=1.",
#                                fc=(1, 0, 0), ec=(1., 0, 0),
#                                patchA=None,
#                                relpos=(0.2, 0.8),
#                                connectionstyle="arc3,rad=-0.1"), size=10)

        if CCDMap!=[]: 
            ax1.scatter(CCDX, CCDY ,s=50, color="blue" , marker='o', alpha = 0.5, label='Model Data')

            if booLabels:
#                fullLabel = ['('+str(CCDX[x])+', '+str(CCDY[x])+')'+str(CCDLambda[x]) for x in np.arange(len(CCDX))]
                fullLabel = [str(CCDLambda[x]) for x in np.arange(len(CCDLambda))]
                
                for x, y, label in zip(CCDX, CCDY, fullLabel ):
                    plt.annotate(
                        label, 
                        xy = (x, y), xytext = (0,-25),
                        textcoords = 'offset points', ha = 'right', va = 'bottom',
                        bbox = dict(boxstyle = 'round,pad=0.5', fc = 'white', alpha = 0.9))
#                        arrowprops = dict(arrowstyle="wedge,tail_width=1.",
#                                    fc=(0, 0, 1), ec=(1., 1, 1),
#                                    patchA=None,
#                                    relpos=(0.2, 0.8),
#                                    connectionstyle="arc3,rad=-0.1"), size=7)                
        
        plt.legend()
        if title=='':title = str(len(imageMapX))+' point(s) found in the calibration image'
        plt.title(title) 
        plt.axis([-imWidth/2 * canvasSize, imWidth/2 * canvasSize, -imHeight/2 * canvasSize, imHeight/2 * canvasSize])

#        plt.draw()
        plt.show()
        
def do_load_spec(specXMLFileName):        
    global Beams, Optics, Cameras, a, b
    
    Beams, Optics, Cameras, a, b = wt.xml.read_all(specXMLFileName)
    
#do_load_spec('hermes.xml')