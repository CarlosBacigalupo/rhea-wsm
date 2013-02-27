#Can plot orders from  146 to 73 (about 390 to 795nm). If the wavelength range just above does not cover the orders selected here, this code currently fails!
minOrder=80
maxOrder=100
deltaOrder=1
maxOrder+=deltaOrder
pixelSize= 5.4 




#Options
SEDModeFlat=0
SEDModeRandom=1
SEDModeSolar=2
SEDModeFile=3
SEDModeCalib=4

CCDMapX=0
CCDMapY=1
CCDMapLambda=2
CCDMapIntensity=3
CCDMapOrder=4

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