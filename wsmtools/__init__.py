"""
Python module with useful functions to be used by the Wavelength Scale Model code. 

Functions:

doSEDMap       ===> Loads an SED Map of different kinds
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


def nkzfs8(Lambda):
    '''Function that calculates the refractive index of the prism for a given wavelength'''
    x = float(Lambda)
    n = np.sqrt(1 + (1.62693651*x**2)/(x**2-0.010880863) + 0.24369876*x**2/(x**2-0.0494207753) + 1.62007141*x**2/(x**2-131.009163) )
    
    return n
    
def n(Lambda, t=18, p=101325):
    '''Function that calculates the refractive index of air for a given wavelength'''
    n = (0.0472326 * (173.3 - (1/Lambda)**2)**(-1))+1
    
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
        z=np.sin(x*np.pi)/(x*np.pi)
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
