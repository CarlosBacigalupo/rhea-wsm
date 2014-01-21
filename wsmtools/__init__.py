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
extract_order   ===> Extracts an order from the spectrum
find_nearest   ===> finds the nearest element of array to a given value
fftshift1D     ===> shifts an image by sub-pixel amounts in one direction
fftshift       ===> shifts an image by sub-pixel amounts in 2D
   
"""
import numpy as np
import os, pyfits, random
import pylab as pl
from astLib import astSED
from matplotlib import *
import time
import pylab as plt

#Custom Packages
import xml_parser as xml
import image_calibration as ic
import image_analysis as ia
import optics
from constants import *
import wsm


def calculate_from_Y(yRange, fX, fLambda):
    
    newY=np.arange(min(yRange),max(yRange)) #Create continuous array from imHeight       
    
    newX=fX(newY) #Creates newX from continuous newY         
    
    #clean results
#     nanMap= np.isnan(newX)
#     newX=newX[-nanMap]
#     newY=newY[-nanMap]            
    
    Lambdas=fLambda(newY)   #Creates lambdas from continuous newY
    
    return newX, newY, Lambdas
            
def gauss(x, p): 
    #Returs a gaussian probability distribution function based on a mean and standard deviation for a range x
    # p[0]==mean, p[1]==stdev
    return 1.0/(p[1]*np.sqrt(2*np.pi))*np.exp(-(x-p[0])**2/(2*p[1]**2))

def fit_errors(p, args):
    '''
    do_ccd_map will return these vectors in a random order. 
    We assume that there are no rounding errors (probably a hack?)
    and will use a floating point == to identify the (x,y) corresponding
    to each wavelength input.
    NB A much better idea would be to re-write main without vstack but
    instead with a neatly created output structure that is defined at the start.
    i.e. when going from SEDMap to SEDMapLoop the information on which line in
    the input each wavelength came from was lost.
    '''
    
    SEDMap = args[0]
    modelXMLFile = args[1]
    calibrationDataFileName = args[2]
    activeCameraIndex = args[3]    
    
    CCDX_c, CCDY_c, lambda_c, xSig_c, ySig_c = ia.read_full_calibration_data(calibrationDataFileName) #reads calibration points coordinates

    CCDMap = wsm.do_ccd_map(SEDMap, modelXMLFile, activeCameraIndex, p_try = p)
    CCDX_m = CCDMap[CCD_MAP_X]
    CCDY_m = CCDMap[CCD_MAP_Y]
    CCDLambda_m = CCDMap[CCD_MAP_LAMBDA]

##    convert model output to CCD coordinates
#    imWidth = Cameras[activeCameraIndex][CamerasWidth]
#    imHeight = Cameras[activeCameraIndex][CamerasHeight]
#    CCDX_c -= imWidth/2
#    CCDY_c -= imHeight/2
        
    #Loop through the input wavelengths, and find the closest output.
    #The best model fits in x and y (out of 2 options) is called x_best and y_best
    x_best = CCDX_c.copy()
    y_best = CCDY_c.copy()
    for k in range(0,len(lambda_c)):
        ix, = np.where(lambda_c[k] == CCDLambda_m)
        if (ix.size == 0):
            x_best[k]=0
            y_best[k]=0
        else:
            best = ix[np.argmin(np.abs(CCDY_m[ix] - CCDY_c[k]))]
            #The following lines have a potential bug because they assume there is only one "best"
            x_best[k] = CCDX_m[best]
            y_best[k] = CCDY_m[best]
         
#    print x_best + 1663
#    print y_best + 1252
#    print x_best - x
#    print y_best - y
#    print np.hstack([(x_best-x)/xSig,(y_best - y)/ySig])
    xSig_c[xSig_c==0]=1
    ySig_c[ySig_c==0]=1
    
    
    diff_model_calib = np.hstack([(x_best-CCDX_c)/xSig_c,(y_best - CCDY_c)/ySig_c]) #, lambda_c
    print np.sum(diff_model_calib)
    
    return diff_model_calib

def fitting_stats(fit):


    print 'Fitting Stats'
    print ''
    
    print 'Solution Array'
    print fit[0]
    print ''
    
    dict = fit[2]
    
    
    for item in dict:
        print item 
        if item=='nfev':
            desc = 'The number of function calls'
        elif item=='fvec':
            desc = 'The function evaluated at the output'
        elif item=='fjac':
           desc = 'A permutation of the R matrix of a QR factorization of the final approximate Jacobian matrix, stored column wise. Together with ipvt, the covariance of the estimate can be approximated.'
        elif item=='ipvt':
            desc = 'An integer array of length N which defines a permutation matrix, p, such that fjac*p = q*r, where r is upper triangular with diagonal elements of nonincreasing magnitude. Column j of p is column ipvt(j) of the identity matrix.'
        elif item=='qtf':
            desc = 'The vector (transpose(q) * fvec)'

        print desc 
        print dict[item]
        print


