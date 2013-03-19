from pyraf import iraf
import os
import numpy as np

def analyze_image(image_filename='test.fits', image_map_filename='image_map.txt'):
    iraf.noao(_doprint=0)     # load noao
    iraf.digiphot(_doprint=0) # load digiphot
    iraf.apphot(_doprint=0)   # load apphot
    
    
    
    #iraf.daofind.setParam('image',FitsFileName)        #Set ImageName
    iraf.daofind.setParam('verify','no')            #Don't verify
    iraf.daofind.setParam('interactive','no')        #Interactive
    iraf.daofind.setParam('datamin','10000')        #Min Good Value
    iraf.daofind.setParam('fwhmpsf','2.5')        #Interactive   
#    iraf.daofind.setParam('interactive','no')        #Interactive

    
    check_if_file_exists(image_map_filename)
    try: iraf.daofind(image = image_filename, output = image_map_filename)
    except Exception: return 0
#    brightest_star_info = find_brightest_star(outfile)
    return

def find_brightest_star(readinfile):
    try: starfile = open(readinfile)
    except Exception: return 'ERROR' # <-- change this to returning a number
    startemp = starfile.readlines()
    brighteststar = 50
    xpixel = 0
    ypixel = 0
    for lines in startemp:
        if lines[0][0] != '#': #don't want the comments
            linetemp = str.split(lines)
            #print linetemp
            if 1: #float(linetemp[2]) < brighteststar:
                starmag = 1 #float(linetemp[2])
                xpixel = float(linetemp[0])
                ypixel = float(linetemp[1])
                starsharp = float(linetemp[3])
    return [starmag, xpixel, ypixel, starsharp]

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
    
    
                              