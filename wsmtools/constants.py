global Beams, Optics, Cameras



#Can plot orders from  146 to 73 (about 390 to 795nm). 
minOrder=67
maxOrder=146
deltaOrder=1


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
CCD_MAP_BEAMID = 5

#SEDMap
SEDMapLambda=0
SEDMapIntensity=1

#Optics types
OpticsBoundary=0
OpticsRGrating=1
OpticsVPHGrating=2

OpticsCoords1=0
OpticsCoords2=1
OpticsType=2
OpticsN=3
OpticsGPeriod=4
OpticsGBlAngle=5

CamerasName=0
CamerasFLength=1
CamerasWidth=2
CamerasHeight=3
CamerasPSize=4
CamerasMinLambda=5
CamerasMaxLambda=6
CamerasDistortion1=7
CamerasDistortion2=8
CamerasDistortion3=9
CamerasDistortionCenterX=10
CamerasDistortionCenterY=11

#File info
SPEC_PATH = 'spectrographs/'
SPEC_BKP_PATH = 'spectrographs/bkp/'
SPECTRUM_PATH = 'spectra/'
FITS_PATH = 'fits/'
TEMP_PATH = 'out_files/'
SEXTRACTOR_PATH = '/usr/local/bin/'
APP_PATH = '/Users/Carlos/Documents/workspace/wsm/'
DATA_SIMULATOR_PATH = '/Users/Carlos/Documents/HERMES/dataSimulator-1.64/dataSimulator/'
DATA_SIMULATOR_INPUT = 'ems/'
WORKING_PATH = APP_PATH

