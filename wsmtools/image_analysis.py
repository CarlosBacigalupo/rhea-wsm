import numpy as np
import subprocess
import pyfits as pf
import pylab as plt
# import matplotlib.cm as cm

# from constants import *
import constants as c
# import image_calibration as ic

def analyse_image_sex(arcFile, sexParamFile):
    #finds peaks using sextractor, returns program output and generated data filename
    
    
    #adds sex_ to the arc filename for sextractor output
    outputFileName =  'sex_' + ''.join(arcFile).split('.')[0]+ '.txt'
      
    os_command = c.SEXTRACTOR_PATH +'sex ' + arcFile + ' -c ' + sexParamFile
    os_command += ' -CATALOG_NAME ' + outputFileName
#     os_command = '/usr/local/bin/sex'
#     os.system(os_command)
    proc = subprocess.Popen([os_command,c.SEXTRACTOR_PATH], stdout=subprocess.PIPE, shell=True)
    out, err_result = proc.communicate()

    # todo check results of os call and pass err_result
    return out, err_result, outputFileName

def identify_imageMapWavelength_avg(CCDMap, calibrationMap, booAvgAdjust = False):
    #todo turn this into a full array calculation
    
    if booAvgAdjust:
        avgImageX = np.average(calibrationMap[:,c.calibrationMap.x])
        avgImageY = np.average(calibrationMap[:,c.calibrationMap.y])
        avgCCDX = np.average(CCDMap[:,c.CCDMap.x])
        avgCCDY = np.average(CCDMap[:,c.CCDMap.y])
        CCDXShifted = CCDMap[:,c.CCDMap.x] + (avgImageX - avgCCDX)
        CCDYShifted = CCDMap[:,c.CCDMap.y] + (avgImageY - avgCCDY)
        CCDMap[:,c.CCDMap.x] = CCDXShifted 
        CCDMap[:,c.CCDMap.y] = CCDYShifted
        
        
    for i in range(len(calibrationMap[:,c.calibrationMap.x])):
        
        distance_array = np.sqrt((CCDMap[:,c.CCDMap.x]-calibrationMap[i,c.calibrationMap.x])**2+(CCDMap[:,c.CCDMap.y]-calibrationMap[i,c.calibrationMap.y])**2)
        closest_point = np.min(distance_array)
        closest_point_index = np.where(distance_array==closest_point)[0][0]       
        calibrationMap[i,c.calibrationMap.wavelength] = CCDMap[closest_point_index,c.CCDMap.wavelength]
            
    return calibrationMap

def identify_imageMapWavelength_manual(SEDMap, CCDX, CCDY, CCDLambda, imageMapX, imageMapY, imageMapWavelength, backImageFileName, Cameras, canvasSize = 1):    
    
    #Plot
    imWidth = int(Cameras[0][CamerasWidth])
    imHeight = int(Cameras[0][CamerasHeight])
    im = pf.getdata(FITS_PATH + backImageFileName) 
    imNorm = ic.normalise_image(im) 
    
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    plt.imshow(imNorm, extent=[-imWidth/2 , imWidth/2 , imHeight/2 , -imHeight/2])
    plt.set_cmap(cm.Greys_r)
    
    #Model Points
    ax1.scatter(CCDX, CCDY ,s=50, color="blue" , marker='o', alpha = 0.5, label='Model Data')
    fullLabel = [str(CCDLambda[x]) for x in np.arange(len(CCDLambda))]    
    for x, y, label in zip(CCDX, CCDY, fullLabel ):
        plt.annotate(
            label, 
            xy = (x, y), xytext = (0,-25),
            textcoords = 'offset points', ha = 'right', va = 'bottom',
            bbox = dict(boxstyle = 'round,pad=0.5', fc = 'blue', alpha = 0.9))
    
    plt.ion()
    
    #Calibration Points
    ax1.scatter(imageMapX, imageMapY ,s=50, color="red" , marker='o', alpha = 0.5, label='Calibration Data')       
    for i in range(len(imageMapX)):
        x = imageMapX[i]
        y = imageMapY[i]
        label = imageMapWavelength[i]
        plt.annotate(
            label, 
            xy = (x, y), xytext = (0,-25),
            textcoords = 'offset points', ha = 'right', va = 'bottom',
            bbox = dict(boxstyle = 'round,pad=0.5', fc = 'red', alpha = 0.9))
        plt.draw()
               
        default = imageMapWavelength[i]
        userInput = raw_input("Enter a string (default: %s):\n" % default)
        if userInput:
            imageMapWavelength[i] = userInput
        else:
            imageMapWavelength[i] = default

#    plt.legend()
#    title = str(len(imageMapX))+' point(s) found in the calibration image'
#    plt.title(title) 
#    plt.axis([-imWidth/2 * canvasSize, imWidth/2 * canvasSize, -imHeight/2 * canvasSize, imHeight/2 * canvasSize])
    return imageMapWavelength

def read_image_map_sex(imageFilename):
    #reads sextractor output and returns as a 2d array x,y,wavelength,sigx,sigy
    arcFileMap = np.loadtxt(imageFilename, usecols = [0,1,3])
    arcFileMap = np.insert(arcFileMap, 2, 0, axis=1)
    arcFileMap = np.insert(arcFileMap, 3, arcFileMap[:,3], axis=1)
    return arcFileMap 

def read_full_calibration_data(calibrationDataFileName):
    #Reads the calibration data from a calibration format txt file 
    #Returns: Nx5 np array with x, y, wavelength, xSig and ySig
    calibrationMap =  np.loadtxt(calibrationDataFileName)

    if calibrationMap.shape[0]>0:
        calibrationMap = np.array(calibrationMap)
    else:
        calibrationMap = np.array([])
   
    return calibrationMap

 
def extract_order(x,y,image, booShowImage = False):
    
    flux=np.zeros(len(y)) 
    
    if type(image)==np.ndarray:
        im = image
        imWidth = image.shape[1]
        imHeight = image.shape[0]
    elif type(image)==str:
        #Grab image data
        hdulist = pf.open(image)
        im = pf.getdata(image)     
        imWidth = hdulist[0].header['NAXIS1']
        imHeight = hdulist[0].header['NAXIS2']
    
    
#     x += imWidth/2
#     y += imHeight/2
    
    for k in range(len(y)):
        x_int = round(x[k])
        #if e.g. x = 10.4, this goes from 5 to 15, and xv at the center is 0.4
        in_image_temp = im[y[k],x_int-10:x_int+11]
#        in_image_temp=im[y[k],x_int-5]        
#        for i in range(-4,6):
#            in_image_temp = np.hstack((in_image_temp,im[y[k+i-4],x_int+i]))          
        in_image_temp[in_image_temp < 0] = 0
        xv = np.arange(-10,11)  + x_int - x[k]
        flux[k] =  np.sum(in_image_temp * np.exp(-(xv/3.5)**4))

        
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

    if booShowImage:
        
        fig = plt.figure()
        ax1 = fig.add_subplot(111, axisbg = 'black')
        im[im<0] = 0
        im[np.isnan(im)] = 0
        im /= im.max()
#         im = np.sqrt(im)
        im = np.log10(im) 
        plt.imshow(im, origin = 'lower' )
        plt.set_cmap(plt.cm.Greys_r)
        
        ax1.scatter(x, y ,s=8 , color = 'green', marker='o', alpha =.5)
        plt.title('Extraction Mask')
        plt.ylabel('pixels')
        plt.xlabel('pixels')
        plt.show()
        
    return flux