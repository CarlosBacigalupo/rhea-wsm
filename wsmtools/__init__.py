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
# import os, pyfits, random
# from matplotlib import *
# import time
# import pylab as plt

#Custom Packages
import xml_parser as xml
import image_calibration as ic
import image_analysis as ia
import optics
import constants as c

# from constants import *

# import wsm


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


