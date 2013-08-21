import wsm
import numpy as np
from constants import *
from scipy.constants import codata

def fq_shift(inputFile, deltaV, outputName = '', outputPath = ''):
    '''
    inputFile dat file from dataSimulator
    deltaV radial velocity change. (positive towards the observer)
    '''
    
    speedOfLight = codata.physical_constants['speed of light in vacuum'][0]
    
    shift = deltaV/speedOfLight
    
    SEDMap = np.loadtxt(DATA_SIMULATOR_PATH + DATA_SIMULATOR_INPUT + inputFile)
    print 'file read'
    SEDMap = SEDMap.transpose()    
    SEDMap[0] = SEDMap[0]*(1-shift)
    
    
    if outputPath=='':
        outputPath = DATA_SIMULATOR_PATH + DATA_SIMULATOR_INPUT
         
    if outputName=='':
        outputName = 'S' + str(deltaV) + '_' + inputFile
         
    theFile = open(outputPath + outputName , 'w')
    SEDMap = SEDMap.transpose()
    
    for item in SEDMap:  
        theLine = fmt_dataSim(item[0]) + ' ' + fmt_dataSim(item[1])    
        theFile.write(theLine+'\n')

        
def fmt_dataSim(valueIn):     
    a = ('%.7e' % valueIn)
    b = a.split('e')     
    valueOut = str(b[0]) + 'e' + str('%+04d' % int(b[1]))
    
    return valueOut
    
fq_shift('4750_15_m25p00.ms.spectra', 10000)