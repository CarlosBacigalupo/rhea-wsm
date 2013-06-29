"""
Python module with useful functions to be used by the Wavelength Scale Model code. 

Functions:

do_SED_map     ===> Loads an SED Map of different kinds
nkzfs8         ===> Calculates refractive index of prism for a given wavelength
n              ===> Calculates refractive index of air for a given wavelength
Snell3D        ===> Computes new direction of vector as it goes across a surface
Grating        ===> Computes new direction of vector as it reflects off a grating
RayTrace       ===> Computes location on chip of vector as it progresses through spectrograph for a wavelenth
Intensity      ===> Retrieves or calculates the expected relative intensity based on distance from the central lambda
extractOrder   ===> Extracts an order from the spectrum
find_nearest   ===> finds the nearest element of array to a given value
fftshift1D     ===> shifts an image by sub-pixel amounts in one direction
fftshift       ===> shifts an image by sub-pixel amounts in 2D
   
"""
import numpy as np
import os, pyfits, random
import pylab as pl
from astLib import astSED
from constants import *
from matplotlib import *
import image_calibration as ic
import time

def ccd_loop(SEDMap, Beam, Optics, stheta, fLength): #, intNormalize,Interpolate=False,BackImage='',GaussFit=False):
    ''' todo
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
       
    dataOut=np.zeros(5)
    CCDX = CCDY = CCDLambda = CCDIntensity = CCDOrder = np.array([])
    
    blaze_angle = stheta #Approximately np.arctan(2)

    #Retrieves max and min lambdas for intensity calculation
    minLambda=min(SEDMap[SEDMapLambda])
    maxLambda=max(SEDMap[SEDMapLambda])
#    allFlux=np.array([0])
#    allLambdas=np.array([0])
    
    #Navigates orders within the range given   
    for nOrder in range(minOrder, maxOrder, deltaOrder):

        #loop lambda for current order
        for i in np.arange(len(SEDMap[SEDMapLambda])): 
            
            Lambda=SEDMap[SEDMapLambda][i]
            inI=SEDMap[SEDMapIntensity][i]
                        
            #rhea hack GPeriod from optics, todo
            GPeriod=31.50321471
            
            #the wavelength range is from the blaze wavelength of the next order and the blaze wavelength of the previous order
            if (Lambda >= abs(2*GPeriod*np.sin(blaze_angle)/(nOrder+1)) and Lambda <= abs(2*GPeriod*np.sin(blaze_angle)/(nOrder-1))):

                #Computes the unit vector that results from the optical system for a given wavelength and order
                #This is the actual tracing of the ray for each wavelength             
                start_time = time.time()
                v, isValid = ray_trace_flex(Beam, Lambda, nOrder, Optics, blaze_angle)
                elapsed_time = time.time() - start_time
                print 'Elepased time: ' + str(elapsed_time)
                
                if isValid: #no errors in calculation, within 1 order of blaze wavelength and beam makes it back through the prism
                    x=v[0]*fLength*1000/pixelSize # x-coord in focal plane in pixels
                    z=v[2]*fLength*1000/pixelSize # z-coord in focal plane in pixels
        
                    outI=Intensity(Lambda, minLambda, maxLambda)    
                    
                    #Add results to output arrays
                    if len(CCDX)==0:
                        CCDX = np.array([x])
                        CCDY = np.array([z])
                        CCDLambda = np.array([Lambda])
                        CCDIntensity = np.array([inI*outI])
                        CCDOrder = np.array([nOrder])
                    else:
                        CCDX = np.append(CCDX,[x],0)
                        CCDY = np.append(CCDY,[z],0)
                        CCDLambda = np.append(CCDLambda,[Lambda],0)
                        CCDIntensity = np.append(CCDIntensity,[inI*outI],0)
                        CCDOrder = np.append(CCDOrder,[nOrder],0)
        
    return CCDX, CCDY, CCDLambda, CCDIntensity, CCDOrder

def ray_trace_flex(Beam, Lambda, nOrder, Optics, blaze_angle):
    ''' Traces a beam through the spectrograph. 
    Spectrograph frame of reference, from the opposite end of the camera looking at the camera
    x=to the right, y=to camera, z=up
    u*=beam, n*=surface normals
    s=grating, perp to the grooves.
    l=grating, parallel to the grooves.
    d=blaze period'''
    
    v = Beam
    
    #loops through optics array
    for i in range(len(Optics)):

        if i==0: n_i=n(Lambda,'air')
        optType=Optics[i][OpticsType]
        n_r=n(Lambda,Optics[i][OpticsN])
        
        if optType==OpticsBoundary:
            surfNormal=Optics[i][OpticsCoords1]           
            v_out = Snell3D(n_i, n_r, v, surfNormal)
        elif optType==OpticsRGrating:
            isValid=False
            GPeriod=Optics[i][OpticsGPeriod]
            s=Optics[i][OpticsCoords1]
            l=Optics[i][OpticsCoords2]
            v_out, isValid = Grating(v, s, l, nOrder, Lambda, GPeriod)
#            else:
#                return v, isValid
            
        n_i = n_r
        v = v_out

#                           
#        """Vector transform due to first surface"""
#        u = Snell3D(nAir, nPrism, u, n1)
#                         
#        """Vector transform due to second surface"""
#        u = Snell3D(nPrism, nAir, u, n2)
#        
#        """grating dispersion"""              
#        u, isValid = Grating(u, s, l, nOrder, Lambda, d)
#        
#        if isValid:
#            """Vector transform due to third surface"""
#            u = Snell3D(nAir, nPrism, u, n4)
#                       
#            """Vector transform due to fourth surface"""
#            u = Snell3D(nPrism, nAir, u, n5)
            
    return v, isValid

def identify_image_map_lambda(image_map_x, image_map_y, back_image):
    
    plot(image_map_x,image_map_y)

def load_image_map(image_filename='test.fits', image_map_filename='image_map.txt'):
    

    
    image_map_x=[]

    try: image_map_file = open(image_map_filename)
    except Exception: return 'ERROR' # <-- change this to returning a number
    image_map_file_temp = image_map_file.readlines()

    for lines in image_map_file_temp:
        if lines[0][0] != '#': 
            linetemp = str.split(lines)

            if len(image_map_x)==0:
                image_map_x = np.array([float(linetemp[0])])
                image_map_y = np.array([float(linetemp[1])])
            else:
                image_map_x = np.append(image_map_x,[float(linetemp[0])],0)
                image_map_y = np.append(image_map_y,[float(linetemp[1])],0)
                
    return image_map_x, image_map_y

def load_image_map_sex(image_filename='test.fits', image_map_filename='image_map.txt'):
    

    
    image_map_x=[]

    try: image_map_file = open(image_map_filename)
    except Exception: return 'ERROR' # <-- change this to returning a number
    image_map_file_temp = image_map_file.readlines()

    for lines in image_map_file_temp:
        if lines[0][0] != '#': 
            linetemp = str.split(lines)

            if len(image_map_x)==0:
                image_map_x = np.array([float(linetemp[0])])
                image_map_y = np.array([float(linetemp[1])])
            else:
                image_map_x = np.append(image_map_x,[float(linetemp[0])],0)
                image_map_y = np.append(image_map_y,[float(linetemp[1])],0)
                
    return image_map_x, image_map_y

def nkzfs8(Lambda):
    #### absorved by n()
    '''Function that calculates the refractive index of the prism for a given wavelength'''
    x = float(Lambda)
    n = np.sqrt(1 + (1.62693651*x**2)/(x**2-0.010880863) + 0.24369876*x**2/(x**2-0.0494207753) + 1.62007141*x**2/(x**2-131.009163) )
    
    return n
    
def n(Lambda, medium='air', t=18, p=101325):
    
    if medium=='air':
        '''Function that calculates the refractive index of air for a given wavelength'''
        n = (0.0472326 * (173.3 - (1/Lambda)**2)**(-1))+1
    elif medium=='nkzfs8':  #nkzfs8 only for the moment. Should change into a sellmeier eq with material as input parameter
        '''Function that calculates the refractive index of the prism for a given wavelength'''
        x = float(Lambda)
        n = np.sqrt(1 + (1.62693651*x**2)/(x**2-0.010880863) + 0.24369876*x**2/(x**2-0.0494207753) + 1.62007141*x**2/(x**2-131.009163) )

    return n

def Snell3D(n_i,n_r,u,n):
    """Computes the new direction of a vector when changing medium. 
    n_i, n_r = incident and refractive indices"""
 
    u_p = u - np.dot(u,n)*n
    u_p /= np.linalg.norm(u_p)
    
    theta_i = np.arccos(np.dot(u,n))
    
    if n_i*np.sin(theta_i)/n_r<=1:
        theta_f = np.arcsin(n_i*np.sin(theta_i)/n_r)
        u = u_p*np.sin(np.pi-theta_f) + n*np.cos(np.pi-theta_f)    

    return u       

def Grating(u, s, l, nOrder, Lambda, d):
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
#todo fix
#    x = 0.5*(float(Lambda) - 0.5*(float(maxLambda) - float(minLambda)))/(float(maxLambda) - float(minLambda))
#    x = (((float(Lambda) - float(minLambda))/(float(maxLambda) - float(minLambda)))-0.5)*2
#    if x!=0:
#        z=np.sin(x*np.pi)/(x*np.pi)
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

    x += imWidth/2
    y += imHeight/2

    
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

def find_nearest(array,value):
    idx=(np.abs(array-value)).argmin()
    return array[idx], idx

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

def create_inter_func(CCDX, CCDY, CCDLambda, CCDIntensity, CCDOrder, nOrder):
    
    xPlot = CCDX[CCDOrder==nOrder]
    yPlot = CCDY[CCDOrder==nOrder]
    LambdaPlot = CCDLambda[CCDOrder==nOrder]
    
    fLambda = interpolate.interp1d(yPlot, LambdaPlot)
    fX = interpolate.interp1d(yPlot, xPlot, 'quadratic', bounds_error=False)

    return fX, fLambda

def calculate_from_Y(Y, fX, fLambda):
    
    newY=np.arange(np.min(Y),np.max(Y))           
    newX=fX(newY)            
    nanMap= np.isnan(newX)
    newX=newX[-nanMap]
    newY=newY[-nanMap]            
    Lambdas=fLambda(newY)
    
    return newX, newY, Lambdas
            
def gauss(x, p): 
    #Returs a gaussian probability distribution function based on a mean and standard deviation for a range x
    # p[0]==mean, p[1]==stdev
    return 1.0/(p[1]*np.sqrt(2*np.pi))*np.exp(-(x-p[0])**2/(2*p[1]**2))

def read_calibration_data(calibrationFile):
    '''
    Reads the calibration data from a txt file and separates the information into 5 separate variables: x, y, wavelength, xSig and ySig.

    OLD CODE:
    CalibrationMap = np.array([0,0,0,0,0])
     
    for line in open(calibrationFile):
        x = float(str(line).split()[0]) #x
        y = float(str(line).split()[1]) #y
        Lambda = float(str(line).split()[2]) #Wavelength
        xSig = float(str(line).split()[3]) #xSig
        ySig = float(str(line).split()[4]) #ySig
        CalibrationMap = np.vstack((CalibrationMap,np.array([x,y,Lambda,xSig,ySig])))
    
    CalibrationMap=CalibrationMap[1:,]  
    '''

    CalibrationMap=np.loadtxt(calibrationFile)
    return CalibrationMap[:,0], CalibrationMap[:,1], CalibrationMap[:,2], CalibrationMap[:,3], CalibrationMap[:,4]
 
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
    #old mainArgs=['4','0',calibrationFile,'0','0']
    #x,y, wavelist are the positions of the peaks in calibrationFile.
    #x,y,waveList,xsig,ysig = readCalibrationData(calibrationFile)
    #fit is the output, which is the ideal p vector.

    fit = leastsq(main_errors,p_try, args=[4,False,calibrationFile,0,False,False,False,True,'c_noFlat_sky_0deg_460_median.fits',False], full_output=True, factor=factor_try, diag=diag_try)


    return fit

def fit_errors(p, args):
    #main will return these vectors in a random order. 
    #We assume that there are no rounding errors (probably a hack?)
    #and will use a floating point == to identify the (x,y) corresponding
    #to each wavelength input.
    #NB A much better idea would be to re-write main without vstack but
    #instead with a neatly created output structure that is defined at the start.
    #i.e. when going from SEDMap to SEDMapLoop the information on which line in
    #the input each wavelength came from was lost.
    #print p
    
    SEDMap = args[0]
    calibDataFileName = args[1]
    x,y,waveList,xSig,ySig = read_calibration_data(calibDataFileName)

    hdulist = pyfits.open('test.fits')
    imWidth = hdulist[0].header['NAXIS1']
    imHeight = hdulist[0].header['NAXIS2']
    
    x=x-imWidth/2
    y=y-imHeight/2
    
    x_model, y_model, Lambda, Intensity, Order = wsm.do_ccd_map(SEDMap, p)
    
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