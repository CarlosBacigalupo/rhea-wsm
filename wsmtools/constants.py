# global Beams, Optics, Cameras


class sedMode():
    flat=0 #flat (intensity =1)
    random=1 #random intensity
    solar=2 #solar spectrum
    file=3 #flat_file wavelength, intensity
    calib=4 # from calibration file

class sedMap():
    wavelength=0
    intensity=1
    
class CCDMap():
    x = 0
    y = 1 
    wavelength = 2
    intensity = 3
    order = 4
    beamID = 5    

class optics():
    coords1x=0
    coords1y=1
    coords1z=2
    coords2x=3
    coords2y=4
    coords2z=5
    type=6
    n=7
    gPeriod=8
    gBlAngle=9

class media():
    air = 0
    nkzfs8 = 1
    NSF11 = 2
    
class beams():
    pass

class cameras():
    name=0
    fLength=1
    width=2
    height=3
    pSize=4
    minLambda=5
    maxLambda=6
    distortion1=7
    distortion2=8
    distortion3=9
    distortionCenterX=10
    distortionCenterY=11

class opticTypes():
    boundary=0
    eGrating=1
    VPHGrating=2
    
class calibrationMap():
    x = 0
    y = 1
    wavelength = 2
    sigX = 3
    sigY = 4

    

#Constants start here------------------------------------------------#todo separate non constants


#SED_modes
# SED_MODE_FLAT=0
# SED_MODE_RANDOM=1
# SED_MODE_SOLAR=2
# SED_MODE_FILE=3
# SED_MODE_CALIB=4
# 
# CCD_MAP_X = 0
# CCD_MAP_Y = 1 
# CCD_MAP_LAMBDA = 2
# CCD_MAP_INTENSITY = 3
# CCD_MAP_ORDER = 4
# CCD_MAP_BEAMID = 5
# 
# #SEDMap
# SEDMapLambda=0
# SEDMapIntensity=1
# 
# #Optics types
# OpticsBoundary=0
# OpticsRGrating=1
# OpticsVPHGrating=2
# 
# OpticsCoords1=0
# OpticsCoords2=1
# OpticsType=2
# OpticsN=3
# OpticsGPeriod=4
# OpticsGBlAngle=5
# 
# CamerasName=0
# CamerasFLength=1
# CamerasWidth=2
# CamerasHeight=3
# CamerasPSize=4
# CamerasMinLambda=5
# CamerasMaxLambda=6
# CamerasDistortion1=7
# CamerasDistortion2=8
# CamerasDistortion3=9
# CamerasDistortionCenterX=10
# CamerasDistortionCenterY=11

#File location info
# SPEC_PATH = 'spectrographs/'
# SPEC_BKP_PATH = 'spectrographs/bkp/'
# SPECTRUM_PATH = 'spectra/'
# FITS_PATH = 'fits/'
# TEMP_PATH = 'out_files/'
# SEXTRACTOR_PATH = '/usr/local/bin/'
# APP_PATH = '/Users/Carlos/Documents/workspace/wsm/'
# DATA_SIMULATOR_PATH = '/Users/Carlos/Documents/HERMES/dataSimulator-1.64/dataSimulator/'
# DATA_SIMULATOR_INPUT = 'ems/'
# WORKING_PATH = APP_PATH

