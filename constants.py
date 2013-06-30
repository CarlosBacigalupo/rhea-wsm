#Can plot orders from  146 to 73 (about 390 to 795nm). If the wavelength range just above does not cover the orders selected here, this code currently fails!
minOrder=73
maxOrder=146
deltaOrder=1
pixelSize= 5.4 





#Constants start here------------------------------------------------#todo separate non constants
#SED_modes
SED_MODE_FLAT=0
SED_MODE_RANDOM=1
SED_MODE_SOLAR=2
SED_MODE_FILE=3
SED_MODE_CALIB=4

CCD_MAP_X = 0
CCD_MAP_Y = 1 
CCD_MAP_LAMBDA = 2
CCD_MAP_INTENSITY = 3
CCD_MAP_ORDER = 4

#SEDMap
SEDMapLambda=0
SEDMapIntensity=1

#Optics types
OpticsBoundary=0
OpticsRGrating=1

OpticsCoords1=0
OpticsCoords2=1
OpticsType=2
OpticsN=3
OpticsGPeriod=4
OpticsGBlAngle=5

#File info
SPEC_DIR='spectrographs/'
FITS_DIR='fits/'
TEMP_DIR='out_files/'
SEXTRACTOR_DIR='/usr/bin/'