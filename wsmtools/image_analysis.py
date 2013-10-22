import numpy as np
import subprocess
import pyfits
import pylab as plt
import matplotlib.cm as cm

from constants import *
import image_calibration as ic

def analyse_image_sex(calibrationImageFileName, sexParamFile, outputFileName):
      
    os_command = SEXTRACTOR_PATH +'sex ' + calibrationImageFileName + ' -c ' + sexParamFile
    os_command += ' -CATALOG_NAME ' + outputFileName
#     os_command = '/usr/local/bin/sex'
#     os.system(os_command)
    proc = subprocess.Popen([os_command,SEXTRACTOR_PATH], stdout=subprocess.PIPE, shell=True)
    out, err_result = proc.communicate()

    # todo check results of os call and pass err_result
    return err_result, out

def identify_imageMapLambda_avg(SEDMap, CCDX, CCDY, CCDLambda, imageMapX, imageMapY, booAvgAdjust = False):
    #todo turn this into a full array calculation
    imageMapLambda = np.zeros(len(imageMapX))
    
    if booAvgAdjust:
        avgImageX = np.average(imageMapX)
        avgImageY = np.average(imageMapY)
        avgCCDX = np.average(CCDX)
        avgCCDY = np.average(CCDY)
        CCDXShifted = CCDX + (avgImageX- avgCCDX)
        CCDYShifted = CCDY + (avgImageY - avgCCDY)

        #These lines plot the model, calibration and model avg corrected points
        import pylab as ppp
        ppp.scatter(CCDX, CCDY, s=50, color="blue", label='Model Data')
        ppp.scatter(imageMapX, imageMapY, s=50, color="red", label='Calibration Data')
        ppp.scatter(CCDXShifted, CCDYShifted, s=50, color="green", label='Average Corrected Model Data')
        ppp.title('Main Reference Points')
        ppp.legend()
        ppp.show()
        
        CCDX = CCDXShifted 
        CCDY = CCDYShifted
        
        
    for i in range(len(imageMapX)):
        
        distance_array = np.sqrt((CCDX-imageMapX[i])**2+(CCDY-imageMapY[i])**2)
        closest_point = np.min(distance_array)
        closest_point_index = np.where(distance_array==closest_point)[0][0]       
        imageMapLambda[i] = CCDLambda[closest_point_index]
        
#        lambda_distance = abs(SEDMap[SEDMapLambda]-imageMapLambda[i])
#        lambda_closest_point = np.min(lambda_distance)
#        lambda_closest_point_index = np.where(lambda_distance==lambda_closest_point)[0][0]
#        imageMapLambda_match[i] = SEDMap[SEDMapLambda][lambda_closest_point_index]
    
    return imageMapLambda

def identify_imageMapLambda_manual(SEDMap, CCDX, CCDY, CCDLambda, imageMapX, imageMapY, imageMapLambda, backImageFileName, Cameras, canvasSize = 1):    
    
    #Plot
    imWidth = int(Cameras[0][CamerasWidth])
    imHeight = int(Cameras[0][CamerasHeight])
    im = pyfits.getdata(FITS_PATH + backImageFileName) 
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
        label = imageMapLambda[i]
        plt.annotate(
            label, 
            xy = (x, y), xytext = (0,-25),
            textcoords = 'offset points', ha = 'right', va = 'bottom',
            bbox = dict(boxstyle = 'round,pad=0.5', fc = 'red', alpha = 0.9))
        plt.draw()
               
        default = imageMapLambda[i]
        userInput = raw_input("Enter a string (default: %s):\n" % default)
        if userInput:
            imageMapLambda[i] = userInput
        else:
            imageMapLambda[i] = default

#    plt.legend()
#    title = str(len(imageMapX))+' point(s) found in the calibration image'
#    plt.title(title) 
#    plt.axis([-imWidth/2 * canvasSize, imWidth/2 * canvasSize, -imHeight/2 * canvasSize, imHeight/2 * canvasSize])
    print imageMapLambda
    return imageMapLambda

def load_image_map_sex(image_map_filename='calib_out.txt'):
    imageMapX=imageMapY=[]

    try: image_map_file = open(image_map_filename)
    except Exception: return [],[],[],[]
    image_map_file_temp = image_map_file.readlines()

    for lines in image_map_file_temp:
        if lines[0][0] != '#': 
            linetemp = str.split(lines)

            if len(imageMapX)==0:
                imageMapX = np.array([float(linetemp[0])])
                imageMapY = np.array([float(linetemp[1])])
                image_map_sigx = np.array([float(linetemp[3])])
                image_map_sigy = np.array([float(linetemp[3])])

            else:
                imageMapX = np.append(imageMapX,[float(linetemp[0])],0)
                imageMapY = np.append(imageMapY,[float(linetemp[1])],0)
                image_map_sigx = np.append(image_map_sigx,[float(linetemp[3])],0)
                image_map_sigy = np.append(image_map_sigy,[float(linetemp[3])],0)
                
    return imageMapX, imageMapY, image_map_sigx, image_map_sigy 

def read_full_calibration_data(calibrationDataFileName):
    '''
    Reads the calibration data from a txt file and separates the information into 5 separate variables: x, y, wavelength, xSig and ySig.
    '''

    CalibrationMap=np.loadtxt(calibrationDataFileName)
    return CalibrationMap[:,0], CalibrationMap[:,1], CalibrationMap[:,2], CalibrationMap[:,3], CalibrationMap[:,4]
 
def extract_order(x,y,image):
    
    #Don't know what these fluxes are yet!
    flux=np.zeros(len(y))
    flux2=np.zeros(len(y))
    flux3=np.zeros(len(y))

    #Grab image and sizes
    hdulist = pyfits.open(FITS_PATH + image)
    imWidth = hdulist[0].header['NAXIS1']
    imHeight = hdulist[0].header['NAXIS2']
    im = pyfits.getdata(FITS_PATH + image)  
#    x += imWidth/2  #correction for pixel number
#    y += imHeight/2 

    
    for k in range(0,len(y)):
        ###mikes'
        x_int = round(x[k])
        #if e.g. x = 10.4, this goes from 5 to 15, and xv at the center is 0.4
        in_image_temp = im[y[k],x_int-5:x_int+6]
#        in_image_temp=im[y[k],x_int-5]        
#        for i in range(-4,6):
#            in_image_temp = np.hstack((in_image_temp,im[y[k+i-4],x_int+i]))          
        in_image_temp[in_image_temp < 0] = 0
        xv = np.arange(-5,6)  + x_int - x[k]
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
    return flux