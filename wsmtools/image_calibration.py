from pyraf import iraf
import os
import numpy as np
from constants import *
import subprocess


def check_if_file_exists(filename):
    #i = 0 # counter to stop this going on forever
    if os.path.isfile(filename): os.remove(filename)
    return filename

def medianCombine(inFiles, outFilename):
    
    for k in range(0,len(inFiles)):
        im = pyfits.getdata(inFiles[k])
        if (k == 0):
            cube=np.empty([im.shape[0],im.shape[1],len(inFiles)])
        cube[:,:,k] = im
    med_im=np.median(cube,axis=2)
    pyfits.writeto(outFilename,med_im)
    
def averageCombine(inFiles, outFilename):
    
    for k in range(0,len(inFiles)):
        im = pyfits.getdata(inFiles[k])
        if (k == 0):
            cube=np.empty([im.shape[0],im.shape[1],len(inFiles)])
        cube[:,:,k] = im
    med_im=np.average(cube,axis=2)
    pyfits.writeto(outFilename,med_im)
    
def batchMedian():
    
    start = range(1,721,5)
    stop = range(5,721,5)
    os.chdir('/media/win-desktop/6hrs/Darks')
    
#    mink=10
#    maxk=15
#    files = 720/5
    for j in range(720/5-1):
        inFiles = ['Hg_20s_%03d.FIT' % k for k in range(start[j],start[j+1])]
        outFilename = 'med' + str(start[j]) + '_' + str(start[j+1]) + '.fits' 
        medianCombine(inFiles, outFilename)

def batchDark():
    
    files = range(1,721)

    os.chdir('/media/win-desktop/6hrs/Darks')

#    mink=10
#    maxk=15
#    files = 720/5
    for j in len(files):
        
        inFile = pyfits.getdata(files[j],0)
        inFiles = ['Hg_20s_%03d.FIT' % k for k in range(start[j],start[j+1])]
        outFilename = 'med' + str(start[j]) + '_' + str(start[j+1]) + '.fits' 
        medianCombine(inFile, darkFile, outFilename)

def calibrateImage(inFileName, biasFileName, darkFileName, flatFileName, outFileName, pixelmaskFileName=''):
    
    hdulist = pyfits.open(inFileName)
    inExpTime=hdulist[0].header['EXPOSURE']

    
    inFile = pyfits.getdata(inFileName)  

    
    
    #Bias subtract
    biasFile = pyfits.getdata(biasFileName) 
    outFile=inFile-biasFile
    
    #Darks
    hdulist = pyfits.open(darkFileName)
    darkExpTime = hdulist[0].header['EXPOSURE']
    
    darkFile=pyfits.getdata(darkFileName)
    darkFile=darkFile/darkExpTime*inExpTime
    
    outFile=outFile-darkFile
    
    #Flats
    flatFile=pyfits.getdata(flatFileName)
    
    #outFile=outFile/flatFile
    
    #Bad Pixels

    
    #write result
    pyfits.writeto(outFileName, outFile)

def subtractDark(inFileName, darkFileName, outFilename):
    
    inFile = pyfits.getdata(inFileName) 
    darkFile = pyfits.getdata(darkFileName) 
    
    outFile = inFile - darkFile
    
    pyfits.writeto(outFilename, outFile)
    
def normalise_image(inFile):
    
    outFile = inFile/np.max(inFile)
    
    return outFile
    
                              