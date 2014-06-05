import numpy as np
import pylab as plt
import pyfits as pf
import os
from mpl_toolkits.mplot3d import Axes3D

import wsm 

baseDir = '/Users/Carlos/Documents/IZA/reductions/3/'
os.chdir(baseDir)

xml = 'dirhea_v4.xml'
calibFileName='calibMap.txt'
dirhea = wsm.spectrograph()
dirhea.read_xml(xml)

SEDMap = dirhea.create_sed_map(SEDMode=wsm.c.sedMode.flat, minLambda=1.0, maxLambda=2.0, deltaLambda =0.001 )

dirhea.orders = range(40,29,-1)

CCDMap,CCDMapOffset = dirhea.create_ccd_map(SEDMap, xml)

arc1 = np.load('arc1.npy')
arc2 = np.load('arc2.npy')
arc3 = np.load('arc3.npy')
lambda1 = np.load('lambda1.npy')
lambda2 = np.load('lambda2.npy')
lambda3 = np.load('lambda3.npy')

sun1 = np.load('sun1_0.npy')
sun1_n = sun1/np.max(sun1)
sun2 = np.load('sun2_0.npy')
sun2_n = sun2/np.max(sun2)
sun3 = np.load('sun3_0.npy')
sun3_n = sun3/np.max(sun3)

darkSun1 = np.load('darkSun1.npy')
darkSun2 = np.load('darkSun2.npy')
darkSun3 = np.load('darkSun3.npy')

sun1_d = sun1 - darkSun1
sun2_d = sun2 - darkSun2
sun3_d = sun3 - darkSun3
offset, removeCam2, removeCam1 = dirhea.find_offset(arc1, arc2, arc3, lambda1, lambda2, lambda3)
fullScience_d = dirhea.stitch(sun1_d, sun2_d, sun3_d, offset)

for i in range(30, 41):
    x2,y2 = dirhea.extract_order(CCDMapOffset, i, orderFileName='20140502_'+str(i)+'.txt')
    plt.plot(x2,y2)
    plt.show()

# 
# 
# arc1 = np.load('arc1.npy')
# arc2 = np.load('arc2.npy')
# arc3 = np.load('arc3.npy')
# lambda1 = np.load('lambda1.npy')
# lambda2 = np.load('lambda2.npy')
# lambda3 = np.load('lambda3.npy')
# 
# posx = np.array([0, 0, 0]) #20140501
# posy = np.array([150.6, 139.4, 161.2]) #20140501
# posx = np.array([.1, 0, .2]) #20140502
# posy = np.array([149.55, 138.05, 160.7]) #20140502
# posy = np.array([149.7, 138.05, 160.7]) #20140502
# fullImage, x_offset, y_offset = dirhea.stitch(np.sum(arc1,2), np.sum(arc2,2),np.sum(arc3,2), posx, posy)
# 
# SEDMap = dirhea.create_sed_map(SEDMode=wsm.c.sedMode.flat, minLambda=1.4, maxLambda=1.6, deltaLambda =0.001 )
# 
# CCDMap,CCDMapOffset = dirhea.create_ccd_map(SEDMap, xml)
# 
# #need to run offset because ot the triple camera thing
# dirhea.camera.offsetX = x_offset
# dirhea.camera.offsetY = y_offset
# dirhea.offset_CCDMap()
# CCDMapOffset = dirhea.CCDMapOffset
# 
# # wsm.plot().ccd_map(CCDMapOffset, canvasSize=1.2, backImage=fullImage)
# 
# x2,y2 = dirhea.extract_order(CCDMapOffset, 35, fullImage, True)
# # plt.plot(x1,y1)
# plt.plot(x2,y2)
# plt.show()

# import numpy as np
# import pylab as plt
# import pyfits as pf
# import os
# 
# import wsm 
# 
# baseDir = '/Users/Carlos/Documents/IZA/reductions/2/'
# os.chdir(baseDir)
# 
# 
# #instantiate the spectrograph class and load the model xml file into corresponding subclasses
# xml = 'dirhea_v2.xml'
# dirhea = wsm.spectrograph()
# dirhea.read_xml(xml)
# 
# im, wavelengths = dirhea.read_multiple_csv('LaserLines_cam_pos1/')
# peaks = dirhea.find_peaks_im(im, wavelengths)
# np.save('im1',im)
# np.save('wavelengths1', wavelengths)
# print peaks
# 
# im, wavelengths = dirhea.read_multiple_csv('LaserLines_cam_pos2/')
# peaks = dirhea.find_peaks_im(im, wavelengths)
# np.save('im2',im)
# np.save('wavelengths2', wavelengths)
# print peaks
# 
# im, wavelengths = dirhea.read_multiple_csv('LaserLines_cam_pos3/')
# peaks = dirhea.find_peaks_im(im, wavelengths)
# np.save('im3',im)
# np.save('wavelengths3', wavelengths)
# print peaks
# # im2, lambda2 = dirhea.read_multiple_csv('LaserLines_cam_pos2/')
# # peaks2 = dirhea.find_peaks_im(im2, lambda2)
# # print peaks2





# import time
# import numpy as np
# import pylab as plt
# import pyfits as pf
# import os
# import glob
# import wsm 
# import constants as c
# import Image
# from mpl_toolkits.mplot3d import Axes3D
# 
# 
# 
# baseDir = '/Users/Carlos/Documents/IZA/reductions/1/'
# os.chdir(baseDir)
# 
# # xml = 'dirhea_v1.xml'
# xml = 'dirhea_v2.xml'
# # xml = 'rhea_v7.xml'
# im = Image.open('arc_all.png')
# dirhea = wsm.spectrograph()
# dirhea.read_xml(xml)
# 
# # SEDMap = dirhea.create_sed_map(SEDMode=wsm.c.sedMode.flat, minLambda=1.4, maxLambda=1.6, deltaLambda=0.001)
# SEDMap = np.array([ 1.5034,  1.51  ,  1.52  ,  1.5228,  1.555 ,  1.56  ,  1.566 , 1.5213,  1.5228,  1.531 ,  1.5364,  1.54  ,  1.5667,  1.5713, 1.578 ,  1.58 ,1.5   ,  1.503 ,  1.506 ,  1.5299,  1.5372,  1.5404,  1.5437,  1.55])
# SEDMap = np.vstack((SEDMap, np.ones(len(SEDMap)))).transpose()
# # SEDMap =np.array([[1.5,1]])
# 
# #beam
# # dirhea.p[0] = 272 #phi
# # dirhea.p[1] = 90 #theta
# 
# # phi_p = 90
# # theta_p = 30
# # dirhea.p[2] = phi_p
# # dirhea.p[3] = theta_p
# # dirhea.p[4] = phi_p              
# # dirhea.p[5] = theta_p+60.              
# 
# # dirhea.p[6] = 27 #phi  
# # dirhea.p[7] = 124 #theta
# # dirhea.p[8] = 180 #alpha
#             
# CCDMap,CCDMapOffset = dirhea.create_ccd_map(SEDMap, xml,p_try = dirhea.p)            
# 
# if CCDMapOffset.shape[0]>0:
#     wsm.plot().ccd_map(CCDMapOffset,  title = str(dirhea.p), backImage = im)
#     fig = plt.figure()
#     ax = fig.add_subplot(111, projection='3d')
#     for i in range(0,len(dirhea.CCDMap3D[:,0,0]),4):
# #     for i in [0,len(dirhea.CCDMap3D[:,0,0])/2,len(dirhea.CCDMap3D[:,0,0])-1]:
#         x = dirhea.CCDMap3D[i,0,:]
#         y = dirhea.CCDMap3D[i,1,:]
#         z = dirhea.CCDMap3D[i,2,:]
#         ax.plot(x, y, z)
#         

#     x = y = np.arange(-3.0, 3.0, 0.05)
#     X, Y = np.meshgrid(x, y)
#     zs = np.array([fun(x,y) for x,y in zip(np.ravel(X), np.ravel(Y))])
#     Z = zs.reshape(X.shape)
#     ax.plot_surface(X, Y, Z)

#     ax.set_xlabel('X (into breadboard')
#     ax.set_ylabel('Y (into Camera)')
#     ax.set_zlabel('Z (from grating to prism)')
#     ax.axis([ -1000, 5000 , -1000, 5000])
#     plt.show()
#     ax = fig.add_subplot(121)
#     ax.scatter(CCDMapOffset[:,c.CCDMap.x],CCDMapOffset[:,c.CCDMap.y])
#     plt.axis((0,640,0,512))
#     plt.show()
# import time
# import numpy as np
# import pylab as plt
# import pyfits as pf
# import wsm 
# import constants as c
# import sys
# 
# rhea = wsm.spectrograph()
# 
# rhea.baseDir = '/Users/Carlos/Documents/RHEA/reductions/1/'
# sexParamFile = rhea.baseDir + 'rhea.sex'
# arcFile = rhea.baseDir + '20140208220659.fits'
# scienceFile = rhea.baseDir + '20140208220925.fits'
# 
# xml = rhea.baseDir + 'rhea_v1.xml'
# calibFileNoWl = rhea.baseDir + 'calNoWl_20140208220659.txt'
# rhea.read_xml(xml, arcFile)
# 
# # calibFileNoWl = rhea.read_arc_file(arcFile, sexParamFile)
# 
# # calibrationMap = wsm.wt.ia.read_full_calibration_data(calibFileNoWl)
# 
# SEDMap1 = rhea.create_sed_map(SEDMode=wsm.c.sedMode.file, spectrumFileName='/Users/Carlos/Documents/workspace/wsm/spectra/hg_spectrum.txt', minI = 0.2)
# SEDMap2 = rhea.create_sed_map(SEDMode=wsm.c.sedMode.file, spectrumFileName='/Users/Carlos/Documents/workspace/wsm/spectra/ar_spectrum.txt', intNormalize = 1, minI=0.3) #create SEDMap from ar emission file
# SEDMap = np.vstack((SEDMap1, SEDMap2))
# 
# CCDMap,CCDMapOffset = rhea.create_ccd_map(SEDMap, xml)
# 
# # wsm.plot().ccd_map(CCDMapOffset,backImage = arcFile)
# 
# 
# rhea.arcMapMask = np.array([False, False, False, False, False, False,  True,  True, False,
#                             False, False, False, False,  True,  True, False, False,  True,
#                             True,  True, False,  True, False, False, False, False,  True,
#                             False,  True, False, False, False, False,  True])
# 
# # wsm.plot().ccd_map(CCDMapOffset[rhea.arcMapMask],backImage = arcFile)
# 
# calibFile = rhea.assign_initial_wavelength(calibFileNoWl, xml, SEDMap, booAvgAdjust = True)
# calibrationMap = wsm.wt.ia.read_full_calibration_data(calibFile)
# 
# 
# print 'CCD points ',rhea.CCDMap[rhea.arcMapMask].shape[0]
# print 'Calib points',calibrationMap.shape[0]
# # a,b = wsm.plot().calib_points_CCDMap(arcFile, calibrationMap, CCDMapOffset, booLabels = True, canvasSize=1, title = 'Calibration vs Model (wavelength assigned)', booInteractive=0, arcMapMask = rhea.arcMapMask)
# 
# rhea.find_fit(SEDMap, xml, calibrationMap, p_try=[], factorTry=1, diagTry=[], showStats=True, maxfev=1000, booWriteP=True)
# 
# 
# print 'aaa'
# 
# # import wsm_old as ws
# # 
# # import pylab as plt
# # import numpy as np
# # import pyfits
# # import sys
# # 
# # 
# # 
# # import wsm
# # 
# # WORKING_PATH = 'set_11/'
# # sexParamFile = WORKING_PATH + 'rhea.sex'
# # arcFile = WORKING_PATH + '20140130231031.fits'
# # xml = 'rhea_v1.xml'
# # 
# # wsm.read_xml(xml)
# # 
# # calibFileNoWl = wsm.read_arc_file(arcFile, sexParamFile, True)
# # sedMap1 = wsm.create_sed_map(SEDMode=wsm.constants.sedMode.file, spectrumFileName='spectra/hg_spectrum.txt', minI=0.6)
# # sedMap2 = wsm.create_sed_map(SEDMode=wsm.constants.sedMode.file, spectrumFileName='spectra/ar_spectrum.txt', intNormalize = 1, minI=0.9) #create SEDMap from ar emission file
# # sedMap = np.hstack((sedMap1, sedMap2))
# # 
# # wsm.plot().calib_points2(arcFile, calibFileNoWl)
# 
# # a = wsm.create_ccd_map(sedMap, xml)
# # 
# # wsm.plot().ccd_map(a, backImage = arcFile)
# 
# 
# # wsm.assign_initial_wavelength(arcFile, calibFileNoWl, xml, sedMap,  booPlotPoints = True)
# 
# 
# 
# 
# sys.exit()
# 
# def launch_task(task):  
#     global wsm, arcFile, modelXMLFile, modelXMLFile2, outputFileName, sexParamFile, finalOutputFileName, scienceFile, WORKING_PATH, flatFile
#     
#     if task==1: # Plot solar SED
#         a=wsm.do_sed_map(SEDMode = SED_MODE_SOLAR)
#         wsm.do_plot_sed_map(a)
#         
#     elif task==2: #Create and plot CCD map from continuous source 
#         a = wsm.do_sed_map(minLambda=0.3, maxLambda=0.8, deltaLambda=0.0001) 
#         b = wsm.do_ccd_map(a, modelXMLFile)
#         wsm.do_load_spec(modelXMLFile)
#         wsm.do_plot_ccd_map(b, backImage = scienceFile)
#     
#     elif task==3: #Read a calibration fits file
#         '''Extract centroids.
#         Assign wavelengths and uncertainties.
#         User confirm assignment
#         Write output file.'''
#         wsm.do_load_spec(modelXMLFile)
#         wsm.do_read_calibration_file(arcFile, modelXMLFile, outputFileName, booAvgAdjust = True, booPlotInitialPoints = True, booPlotFinalPoints = True)
#        
#     elif task==4: #Find fit
# 
#         activeCamera = 0
#         diagTry = [1.5, 0.5, 0.9 , 0.2, 1.2, 1.2, 1.6, 1, 1.2, 6.7, 50, 50] 
#         showStats = True
#         wsm.do_load_spec(WORKING_PATH + modelXMLFile)
#         
#         SEDMap1 = wsm.do_sed_map(SEDMode=wsm.SED_MODE_FILE, spectrumFileName='spectra/hg_spectrum.txt', minI=0.6)
#         SEDMap2 = wsm.do_sed_map(SEDMode=wsm.SED_MODE_FILE, spectrumFileName='spectra/ar_spectrum.txt', intNormalize = 1, minI=0.9) #create SEDMap from ar emission file
#         SEDMap = np.hstack((SEDMap1, SEDMap2))
#                 
#         p_out = wsm.do_find_fit(SEDMap, WORKING_PATH + modelXMLFile, WORKING_PATH +  finalOutputFileName, 0, showStats = True, diagTry = diagTry)      
#         p_out = p_out[0]
#                 
#         #plot fit to see new model
# #         SEDMap = wsm.do_sed_map()
#         CCDMap = wsm.do_ccd_map(SEDMap, WORKING_PATH + modelXMLFile, p_try = p_out)
#         wsm.do_plot_ccd_map(CCDMap, backImage = WORKING_PATH + arcFile)
#         wsm.do_plot_calibration_points(WORKING_PATH + arcFile, WORKING_PATH + arcFile, CCDMap, booLabels = True, canvasSize=1, title = 'Calibration points vs Model comparison ')
#     
#     elif task==5: #Full loop from calibration file to model
# 
# 
# 
#         #step 1 - Read fits, extract calibration points to finalOutputFileName (n x 5 np.array)        
# 
#         
#         #Create SEDMap from hgar emission
#         SEDMap1 = wsm.do_sed_map(SEDMode=wsm.SED_MODE_FILE, spectrumFileName='hg_spectrum.txt', minI=0.6)
#         SEDMap2 = wsm.do_sed_map(SEDMode=wsm.SED_MODE_FILE, spectrumFileName='ar_spectrum.txt', intNormalize = 1, minI=1) #create SEDMap from ar emission file
#         SEDMap = np.hstack((SEDMap1, SEDMap2))
#         
#         wsm.do_load_spec(modelXMLFile)
#         
#         wsm.do_read_calibration_file(arcFile,  modelXMLFile, 
#                                     outputFileName,  sexParamFile, 
#                                     finalOutputFileName, SEDMap,
#                                     booPlotInitialPoints = True, booPlotFinalPoints = True, 
#                                     booAvgAdjust = False)
#         
#         #step 2 - find best fit to results extracted
# #         diagTry = [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]
#         diagTry = [1,1,1,1,1,1,1,1,1,1,1,0.001,0.1,0.1,1,0.0001]
#         p_out = wsm.do_find_fit(SEDMap, modelXMLFile, finalOutputFileName, 0, diagTry = diagTry, showStats = True)      
#         p_out = p_out[0]
#         
#     
#         #step 3 - plot the optimised solution
# #         CCDMap = wsm.do_ccd_map(SEDMap, modelXMLFile, p_try = p_out)
# #         wsm.do_plot_calibration_points(arcFile, finalOutputFileName, CCDMap, booLabels = False, canvasSize=1, title = 'Calibration points vs Model after fit')
#         
#     
#         #step 4 - extract using model
# #         SEDMap = wsm.do_sed_map()
# #         SEDMap = wsm.do_sed_map(SEDMode=wsm.SED_MODE_FILE, spectrumFileName='thar_spectrum.txt', 
# #                                 intNormalize=1, maxLambda=1) #create SEDMap from ar emission file
# 
#         CCDMap = wsm.do_ccd_map(SEDMap, modelXMLFile, p_try = p_out)  
#         wsm.do_plot_ccd_map(CCDMap, backImage = scienceFile)
# #         flux = wsm.do_extract_order(CCDMap , 87, WORKING_PATH + scienceFile)
# #         wsm.do_plot_flux(flux)
# #         outputFile = WORKING_PATH + 'thar_map.txt'
# #         wsm.do_export_CCDMap(CCDMap, WORKING_PATH + scienceFile, outputFile)  
#         
#     elif task==6: #extract spectrum
#         a = wsm.do_sed_map(minLambda=0.3, maxLambda=0.8, deltaLambda=0.0001) 
#         b = wsm.do_ccd_map(a, modelXMLFile)
#         wsm.do_load_spec(modelXMLFile)
#         c, d = wsm.do_extract_order(b , 90, scienceFile, True)
#         wsm.do_plot_flux(c, d)
#         
#     elif task==7: #plot points from calibration file
#         SEDMap1 = wsm.do_sed_map(SEDMode=wsm.SED_MODE_FILE, spectrumFileName='hg_spectrum.txt', minI=0.6)
#         SEDMap2 = wsm.do_sed_map(SEDMode=wsm.SED_MODE_FILE, spectrumFileName='ar_spectrum.txt', intNormalize = 1, minI=1) #create SEDMap from ar emission file
#         SEDMap = np.hstack((SEDMap1, SEDMap2))
#         CCDMap = wsm.do_ccd_map(SEDMap, modelXMLFile)
#         wsm.do_plot_calibration_points(arcFile, outputFileName, CCDMap=CCDMap, booLabels = True, canvasSize=1, title = 'Calibration points vs Model comparison ')
#     
#     elif task==8: #find fit, play with parameter space 
#         calibrationDataFileName = 'c_hg_rhea_sample1.txt'
#         modelXMLFile = 'rhea_initial.xml'
#         arcFile = 'hg_rhea_sample1.fits'
#         activeCamera = 0
#         showStats = True
#         
#         iter = np.arange(1)
#         dist = np.zeros(len(iter))
#         avgDist = np.zeros(len(iter))
#         
#         SEDMap = wsm.do_sed_map(SEDMode=SED_MODE_FILE, spectrumFileName='hg_spectrum.txt') #create SEDMap from flat Hg emission file
# 
#         plt.ion()
#         for i in np.arange(len(iter)):
#             diagTry = [1.5, 0.5, 0.9 , 0.2, 1.2, 1.2, 1.6, 1, 1.2, 67, 1, 1]  
#             factorTry = 1
#             p_out = wsm.do_find_fit(SEDMap, modelXMLFile, calibrationDataFileName, activeCamera, diagTry = diagTry, factorTry = factorTry, showStats = showStats)      
#             fvec = p_out[2]['fvec']
#             x = fvec[:len(fvec)/2]
#             y = fvec[len(fvec)/2:]
#             dist[i] = sum(np.sqrt(x**2+y**2))
#             avgDist[i] = dist[i]/len(x)
#             
#             plt.clf()
# #            plt.scatter(iter[:i], dist[:i])
#             plt.scatter(iter[:i+1], avgDist[:i+1], color = 'red')
#             plt.grid()
#             plt.draw()           
#             
#         #plot fit to see new model
# #        CCDMap = wsm.do_ccd_map(SEDMap, modelXMLFile, p_try = p_out)
# #        wsm.do_plot_calibration_points(arcFile, calibrationDataFileName, CCDMap, booLabels = True, canvasSize=1, title = 'Calibration points vs Model comparison ')
#         plt.ioff()
#         plt.show()
#         
#     elif task==9: #Create and plot CCD map from continuous source 
#         WORKING_PATH = wsm.APP_PATH + 'set_1/'
#         a = wsm.do_sed_map(minLambda=0.3, maxLambda=0.9, deltaLambda=0.001) 
#         b = wsm.do_ccd_map(a, WORKING_PATH + 'rhea.xml')
#         wsm.do_load_spec(WORKING_PATH + 'rhea.xml')
#         wsm.do_plot_ccd_map(b, backImage='sky_rhea_sample.fits')        
# 
#     elif task==10: #Distortion tests
#         inX = np.arange(-20,21)
#         inX = np.tile(inX,len(inX))
#         inY = np.sort(inX)
#         K=[-0.8,0,0.5]
#         
#         outX, outY = wsm.wt.distort(inX, inY, K, 4, -7)
#         
#         plt.scatter(outX, outY)
#         plt.show()
#         
#     elif task==11: #New spectra test
#        
#        SEDMap = wsm.do_sed_map(SEDMode=wsm.SED_MODE_FILE, spectrumFileName='thar_spectrum.txt') #create SEDMap from ar emission file
#        wsm.do_plot_sed_map(SEDMap)
#         
#     elif task==12: #cross correlation test
#         from mpl_toolkits.mplot3d.axes3d import Axes3D
#         from scipy.signal import correlate2d
#         
#         arcFile = 'hgar.fit'
#         modelXMLFile = 'rhea.xml'
#         outputFileName = 'hg_rhea2.txt'
#         sexParamFile = 'rhea.sex'
#         finalOutputFileName = 'hg_rhea2.txt'
#         scienceFile = 'thar.fit'
#         WORKING_PATH = wsm.APP_PATH + 'set_6/'
# 
#         wsm.do_load_spec(WORKING_PATH + modelXMLFile)
#         
#         #Create SEDMap from hgar emission
#         SEDMap1 = wsm.do_sed_map(SEDMode=wsm.SED_MODE_FILE, spectrumFileName='hg_spectrum.txt', minI=0.6)
#         SEDMap2 = wsm.do_sed_map(SEDMode=wsm.SED_MODE_FILE, spectrumFileName='ar_spectrum.txt', intNormalize = 1, minI=0.9) #create SEDMap from ar emission file
#         SEDMap = np.hstack((SEDMap1, SEDMap2))
# #        SEDMap = wsm.do_sed_map(SEDMode=wsm.SED_MODE_FILE, spectrumFileName='thar_spectrum.txt', 
# #                                 intNormalize=1, minI=0.2, maxLambda=1) #create SEDMap from ar emission file
#         
#         
#         CCDMap = wsm.do_ccd_map(SEDMap, WORKING_PATH + modelXMLFile)#, p_try = p_out)
#         
#         hdulist = pyfits.open(WORKING_PATH + arcFile)
#         im1Width = hdulist[0].header['NAXIS1']
#         im1Height = hdulist[0].header['NAXIS2']
#         im1Width_range = range(im1Width)
#         im1Height_range = range(im1Height)
#         X1,Y1  = np.meshgrid(im1Width_range,im1Height_range)
#         im1 = pyfits.getdata(WORKING_PATH + arcFile)  
#         
# #         wsm.do_read_calibration_file(WORKING_PATH + arcFile, WORKING_PATH + modelXMLFile, 
# #                                      WORKING_PATH +  outputFileName, WORKING_PATH +  sexParamFile, 
# #                                      WORKING_PATH +  finalOutputFileName, SEDMap,
# #                                      booPlotInitialPoints = True, booPlotFinalPoints = True)
#         #Create empty grid
# #         im2Width = max(CCDMap[wsm.CCD_MAP_X]) - min(CCDMap[wsm.CCD_MAP_X]) 
# #         im2Height = max(CCDMap[wsm.CCD_MAP_Y]) - min(CCDMap[wsm.CCD_MAP_Y])  
#         im2Width_range = range( int(min(CCDMap[wsm.CCD_MAP_X])), int(max(CCDMap[wsm.CCD_MAP_X])))
#         im2Height_range = range( int(min(CCDMap[wsm.CCD_MAP_Y])), int(max(CCDMap[wsm.CCD_MAP_Y])))        
#         X2,Y2  = np.meshgrid(im2Width_range,im2Height_range)
#         im2 = np.zeros()
#         
#         
# #         for i in range(len(CCDMap[CCD_MAP_X])):
#         sigma = 2.
#         b = 2*sigma**2       
#         xCent=50.
#         yCent=50.
#         im2 = np.exp(-((X2-xCent)**2/b + (Y2-yCent)**2/b))
#         
# #         xCent=100.
# #         yCent=100.
# #         Z2 = np.exp(-((X2-xCent)**2/b + (Y2-yCent)**2/b))
#         a =correlate2d(im1,im2, mode='same')
# #         print a
# #         print a.shape
# #         print np.argmax(np.max(a,axis=0))+1
# #         print np.argmax(np.max(a,axis=1))+1
# #         phi_m = linspace(0, 2*pi, 100)
# #         phi_p = linspace(0, 2*pi, 100)
# #         X,Y = meshgrid(phi_p, phi_m)
# #         Z = flux_qubit_potential(X, Y)/flux_qubit_potential(X, Y)
# #         
# #         fig, ax = plt.subplots()
# #         
# #         p = ax.pcolor(X/(2*pi), Y/(2*pi), Z, cmap=cm.RdBu, vmin=abs(Z).min(), vmax=abs(Z).max())
# #         cb = fig.colorbar(p, ax=ax)
#     
#         fig = plt.figure(figsize=(12,6))
# 
#         ax = fig.add_subplot(1,3,1, projection='3d')
#         ax.plot_surface(X1, Y1, im1, rstride=4, cstride=4, alpha=0.25)
#         
#         ax = fig.add_subplot(1,3,2, projection='3d')
#         ax.plot_surface(X2, Y2, im2, rstride=4, cstride=4, alpha=0.25)
#         
# #         ax = fig.add_subplot(1,3,3, projection='3d')
# #         ax.contour(range(a.shape[0]),range(a.shape[1]), a)
# 
#         plt.show()
#         
#     elif task==13: #New calibration bit
#         
#         arcFile = 'hgar.fit'
#         modelXMLFile = 'rhea.xml'
#         spectrumFileName='hg_spectrum.txt'
#         calibFile = 'calib_hg_rhea.txt'
#         WORKING_PATH = wsm.APP_PATH + 'set_7/'
# 
#         #step 1 - Read fits, extract calibration points to finalOutputFileName (n x 5 np.array) 
#         wsm.do_load_spec(WORKING_PATH + modelXMLFile)
#         
#         wsm.do_create_calibFile( WORKING_PATH + arcFile, 
#                                  WORKING_PATH +  modelXMLFile, 
#                                  WORKING_PATH + calibFile, 
#                                  WORKING_PATH + spectrumFileName)
#         
#         
#         #Create SEDMap from hgar emission
#         SEDMap = wsm.do_sed_map(SEDMode=wsm.SED_MODE_FILE, spectrumFileName = WORKING_PATH + 'hg_spectrum.txt', minI=0.6)
#         
#         wsm.do_read_calibration_file(WORKING_PATH + arcFile, WORKING_PATH + modelXMLFile, 
#                                      WORKING_PATH +  outputFileName, WORKING_PATH +  sexParamFile, 
#                                      WORKING_PATH +  finalOutputFileName, SEDMap,
#                                      booPlotInitialPoints = True, booPlotFinalPoints = True)
#         
#         #step 2 - find best fit to results extracted
#         diagTry = [1,1,1,1,1,1,1,1,1,1,1,0.001,0.1,0.1,1,0.0001]
#         p_out = wsm.do_find_fit(SEDMap, WORKING_PATH + modelXMLFile, WORKING_PATH +  finalOutputFileName, 0, diagTry = diagTry, showStats = True)      
#         p_out = p_out[0]
#         
#     
#         #step 3 - plot the optimised solution
#         CCDMap = wsm.do_ccd_map(SEDMap, WORKING_PATH + modelXMLFile, p_try = p_out)
#         wsm.do_plot_calibration_points(WORKING_PATH +  arcFile, WORKING_PATH +  finalOutputFileName, CCDMap, booLabels = False, canvasSize=1, title = 'Calibration points vs Model after fit')
#     
#     elif task==14: #extract full spectrum + 
#         import MyAstroLib as MAL
#         import math
#         
#         c = MAL.const.c
#         
#         SEDMap = wsm.do_sed_map(minLambda=0.5, maxLambda=0.6, deltaLambda=0.001) 
#         
#         CCDMapStar = wsm.do_ccd_map(SEDMap, modelXMLFile)
#         wsm.do_load_spec(modelXMLFile)
#         specOrderStar, specLambdaStar, specFluxStar = wsm.do_extract_full_spectrum(CCDMapStar, scienceFile)
#         
#         CCDMapFlat = wsm.do_ccd_map(SEDMap, modelXMLFile2)
#         wsm.do_load_spec(modelXMLFile2)
#         specOrderFlat, specLambdaFlat, specFluxFlat = wsm.do_extract_full_spectrum(CCDMapFlat, flatFile)
# 
#         flatLambda, flatFlux = wsm.wt.ic.divide_flat_full(specOrderStar, specLambdaStar, specFluxStar, specOrderFlat, specLambdaFlat, specFluxFlat, booPlot = False)
# #         for i in range(len(flatLambda)):
# #             Q, dRV = MAL.RVS.QdRV(flatLambda[i], flatFlux[i])
# #             print Q, dRV
#             
#         fullLambda, fullFlux = wsm.wt.ic.stitch(flatLambda, flatFlux)
#         
#         np.save(WORKING_PATH + 'all_lambda.npy',fullLambda)
#         np.save(WORKING_PATH + 'all_flux.npy',fullFlux)
#         
#         wsm.do_plot_flux(fullLambda, fullFlux) 
#                
#  
# taskDict = {'plotCCD':2,
#             'findModel':5,
#             'plotCalib':7,
#             'extractOrder':6,
#             'extractAll':14}
#  
#  
#  
# WORKING_PATH = wsm.APP_PATH + 'set_11/'
# sexParamFile = WORKING_PATH + 'rhea.sex'
#  
# outputFileName = WORKING_PATH + 'hg_rhea2.txt'
# finalOutputFileName = WORKING_PATH + 'hg_rhea2.txt'
#  
# modelXMLFile = WORKING_PATH + 'rhea_v1.xml'
# modelXMLFile2 = WORKING_PATH + 'rhea_v3.xml'
#  
# arcFile = WORKING_PATH + '20140130231031.fits'
# scienceFile = WORKING_PATH + '20140130225508.fits'
# # scienceFile = WORKING_PATH + '20131203_Masterflat.fit'
# # scienceFile = WORKING_PATH + '20131203_HgAr.fit'
# # flatFile = WORKING_PATH + '20131203_Masterflat.fit'
#  
# # wsm.do_load_spec(modelXMLFile)
#  
# launch_task(taskDict['findModel'])
# 
# 
# #a=wsm.do_sed_map(minLambda=0.3, maxLambda=0.9)
# #b=wsm.do_ccd_map(a)
# #c=wsm.do_extract_order(b , 87, 'c_noFlat_sky_0deg_460_median.fits')
# #
# #
# #
# #
# #
# #print 'do stuff now'
# #wsm.do_plot_ccd_map(b)
# 
# 
# #SEDMap=wsm.do_sed_map(SEDMode=wsm.SED_MODE_CALIB, specFile='calib_out.txt')
# #p=wsm.do_find_fit(SEDMap, 'calib_out.txt')
# #print p
# #import image_calibration as ic
# 
# #
# #im = wt.load_image_map()
# #print im
# 
# #ic.analyseImage('test.fits', 'a.txt')
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# #
# #a=wsm.do_sed_map(SEDMode=wsm.SED_MODE_CALIB, specFile='Hg_5lines_double.txt')
# # 
# #b=wsm.do_ccd_map(a)
# #print b
# #wsm.do_plot_ccd_map(b)
# 
# 
# 
# 
# #from scipy import stats
# #import numpy as np
# #import ds9
# # 
# ## Make a 2D gaussian image that is stored in a 2D numpy array
# #x = np.arange(-3, 3, 0.1)
# #xx, yy = np.meshgrid(x, x)
# #gauss2d = stats.norm.pdf(xx) * stats.norm.pdf(yy)
# # 
# ## Now open ds9 (this assumes no ds9 instance is yet running)
# #d = ds9.ds9()
# # 
# ## Load up our 2D gaussian
# #d.set_np2arr(gauss2d)
# ## ~/Documents/workspace/rhea-wsm/fits
# ## Zoom to fit
# #d.set('zoom to fit')
# # 
# ## Change the colormap and scaling
# #d.set('cmap bb')
# #d.set('scale log')
# # 
# ## Add a label
# #d.set('regions command {text 30 20 #text="Fun with pyds9" font="times 18 bold"}')
# # 
# ## Now you can play in ds9 to your heart's content.
# ## Check back to see what the current color scale is.
# #print d.get('scale')
# # 
# ## Finally, save your completed image (including regions or labels)
# #d.set('saveimage png my_pyds9_img.png')
# 
# 
# 
# #print "3"
# #import xml_parser
# #
# #
# #p=xml_parser.read_All()
# #
# #print p
# 
# #os.chdir('/media/win-desktop/6hrs/Darks')
# #inFiles = ['Hg_6hrs_20s_D%01d.FIT' % k for k in range(1,6)]
# #outFilename = 'masterDark.fits' 
# #medianCombine(inFiles, outFilename)
# 
# #
# #fit,success=findFit('c_noFlat_Hg_0deg_10s.txt')
# #print fit
# #main(p=fit)
# 
# #main()
# 
# #compareTemp('6hrs/Medians/med1_5.fits','6hrs/Medians/med711_715.fits','6hrs/Results/out1.fits')
# 
# #batchAverage2()
# 
# #batchMedian()
# 
# #inFiles = ['sky_0deg_460_%01d.fits' % k for k in range(1,4)]
# #medianCombine(inFiles, 'sky_0deg_460_median.fits')
# #inFiles = ['sky_0deg_1550_%01d.fits' % k for k in range(1,4)]
# #medianCombine(inFiles, 'sky_0deg_1550_median.fits')
# 
# #inFiles = ['sky_0deg_460_%01d.fits' % k for k in range(1,4)]
# #averageCombine(inFiles, 'sky_0deg_460_average.fits')
# #inFiles = ['sky_0deg_1550_%01d.fits' % k for k in range(1,4)]
# #averageCombine(inFiles, 'sky_0deg_1550_average.fits')
# 
# #calibrateImage('sky_0deg_460_median.fits','mbias.fits','mdark.fits','mflat.fits','c_noFlat_sky_0deg_460_median.fits')
# #calibrateImage('sky_0deg_1550_median.fits','mbias.fits','mdark.fits','mflat.fits','c_noFlat_sky_0deg_1550_median.fits')
# 
# #calibrateImage('Hg_0deg_10s.fits','mbias.fits','mdark.fits','mflat.fits','c_Hg_0deg_10s.fits')
# 
# #calibrateImage('sky_0deg_460_1.fits','mbias.fits','mdark.fits','mflat.fits','c_noFlat_sky_0deg_460_1.fits')
# #calibrateImage('sky_0deg_460_2.fits','mbias.fits','mdark.fits','mflat.fits','c_noFlat_sky_0deg_460_2.fits')
# #calibrateImage('sky_0deg_460_3.fits','mbias.fits','mdark.fits','mflat.fits','c_noFlat_sky_0deg_460_3.fits')
# 
# #calibrateImage('sky_0deg_1550_median.fits','mbias.fits','mdark.fits','mflat.fits','c_noFlat_sky_0deg_1550_median.fits')
# 
# 
# 
# 
# #inFiles = ['dark_a%01d.fits' % k for k in range(1,4)]
# #medianCombine(inFiles, 'dark_a_median.fits')
# #
# #inFiles = ['dark_a%01d.fits' % k for k in range(1,4)]
# #averageCombine(inFiles,  'dark_a_average.fits')
# #
# #
# #subtractDark('arcturus_1min_1.fits','dark_a_average.fits','c_arcturus_1min_1.fits')
# #subtractDark('arcturus_1min_2.fits','dark_a_average.fits','c_arcturus_1min_2.fits')
# #subtractDark('arcturus_1min_3.fits','dark_a_average.fits','c_arcturus_1min_3.fits')
# #subtractDark('arcturus_1min_4.fits','dark_a_average.f#inFiles = ['dark_a%01d.fits' % k for k in range(1,4)]
# #medianCombine(inFiles, 'dark_a_median.fits')
# #
# #inFiles = ['dark_a%01d.fits' % k for k in range(1,4)]
# #averageCombine(inFiles,  'dark_a_average.fits')
# #
# #
# #subtractDark('arcturus_1min_1.fits','dark_a_average.fits','c_arcturus_1min_1.fits')
# #subtractDark('arcturus_1min_2.fits','dark_a_average.fits','c_arcturus_1min_2.fits')
# #subtractDark('arcturus_1min_3.fits','dark_a_average.fits','c_arcturus_1min_3.fits')
# #subtractDark('arcturus_1min_4.fits','dark_a_average.fits','c_arcturus_1min_4.fits')
# #
# #inFiles = ['c_arcturus_1min_%01d.fits' % k for k in range(1,5)]
# #medianCombine(inFiles, 'c_arcturus_1min_median.fits')
# #
# #inFiles = ['c_arcturus_1min_%01d.fits' % k for k in range(1,5)]
# #averageCombine(inFiles, 'c_arcturus_1min_average.fits')its','c_arcturus_1min_4.fits')
# #
# #inFiles = ['c_arcturus_1min_%01d.fits' % k for k in range(1,5)]
# #medianCombine(inFiles, 'c_arcturus_1min_median.fits')
# #
# #inFiles = ['c_arcturus_1min_%01d.fits' % k for k in range(1,5)]
# #averageCombine(inFiles, 'c_arcturus_1min_average.fits')
# 
# 
# #print x(0.502)
# 
# #p=main_errors(p = [ 271.92998622,   91.03999719,   59.48997316,   89.83000496,   89.37002499,   89.79999531,   68.03002976,   64.9399939,     1.15998754,   31.52736851,  200.00000005], SEDMode=4,booPlot=True,specFile='c_noFlat_Hg_0deg_10s.txt')
# 
# #main(p = [279, 90.6, 59, 90, 89, 90, 70.75, 65, 0, 31.64, 210])
# #main(p = [ 271.92998622,   91.03999719,   59.48997316,   89.83000496,   89.37002499,   89.79999531,   68.03002976,   64.9399939,     1.15998754,   31.52736851,  200.00000005]) ####best fit so far....
# ##
# ##print arrayOut
# ##
# #xi = np.arange(0,9)
# #A = np.array([ xi, np.ones(9)])
# ## linearly generated sequence
# #y = [19, 20, 20.5, 21.5, 22, 23, 23, 25.5, 24]
# #w = np.linalg.lstsq(A.T,y)[0] # obtaining the parameters
# #s=1
# 
# 
# #best p for previous fitting:
# #p = [ 271.92998622,   91.03999719,   59.48997316,   89.83000496,   89.37002499,   89.79999531,   68.03002976,   64.9399939,     1.15998754,   31.52736851,  200.00000005]
# #'c_noFlat_Hg_0deg_10s.fits')#,'solar.png')#,'mdark.fits')#,'simple_flat.fits')#
# 
# 
# #Fraunhofer Lines-CCDMap with labels***********************************************************************************
# #identify fraunhofer
# #SED fraunhofer
# #(SEDMode(0=Max, 1=Random, 2=Sun, 3=from specFile, 4=from CalibFile), Plot?, specFile/calibfile, 
# #Normalize intensity? (0=no, #=range), Distort?, Interpolate, PlotCalibPoints, booPlotLabels, plotbackfile, 
# #Gaussian) <-- other options#main(args=['4','1','solar.txt','1','0','0','0','1','c_noFlat_sky_0deg_460_median.fits'])
# #p = [272.31422902, 90.7157937, 59.6543365, 90.21334551, 89.67646101, 89.82098015, 68.0936684,  65.33694031, 1.19265536, 31.50321471, 199.13548823]
# #main(p,SEDMode=4,booPlot=True,specFile='solar.txt')
# #********************************************************************************************************
# 
# 
# #Arcturus callibration**********************************************************************************
# #p_try=[ 272.45928478  , 91.17527895 ,  59.25339245 ,  89.61147631 ,  89.63791054,   89.81178189 ,  68.1669475 ,   63.3555271 ,   1.10221798 ,  31.8848023,  199.70165672] #best
# ##p_try=[250.6, 64.2, 58.3, 77.9, 89.5, 85.5, 58.5, 61.7, 1.6, 33.7, 193.5] #doeesn't work
# ## ok p_try=[279, 90.6, 59, 90, 89, 90, 70.75, 65, 0, 31.64, 210]
# ##p_try = [278.6, 90.6, 61, 90, 91, 90, 70.75, 63.95, 0, 31.9, 210]
# ##p_try = [278.2,90.8,59.0,90.4,88.9,89.7,70.8,65.1,0.758,31.67,203.1]
# #only h-alpha
# #p_try = [ 272.46778143,   91.14782003 ,  59.2663218 ,   89.65582564 ,  89.62544383,   89.76937348  , 68.19750843 ,  63.29873297  ,  1.10998689  , 31.92112407,  199.70153725]
# #
# ##next
# #p_try = [ 272.45928478,   91.17527895  , 59.25339245 ,  89.61147631 ,  89.63791054,  89.81178189 ,  68.1669475 ,   63.3555271,     1.10221798 ,  31.8848023,  199.70165672]
# #
# #diag_try=[1,1,1,1,1,1,1,1,1,0.1,1]
# #factor_try=0.3
# #fit, success=findFit('arcturus.txt',p_try, factor_try, diag_try)
# #print fit
# #main(fit,SEDMode=4,booPlot=True,specFile='arcturus.txt',booPlotCalibPoints=True,plotBackImage='c_arcturus_1min_4.fits')
# #********************************************************************************************************
# 
# 
# #Arcturus plot********************************************************************************************************
# #(beam phi, beam theta, prism1 phi, prism1 theta, prism2 phi, prism2 theta, grating phi, grating theta, grating alpha,
# #         blaze period (microns), focal length(mm), distortion term)
# #(SEDMode(0=Max, 1=Random, 2=Sun, 3=from specFile, 4=from CalibFile), Plot?, specFile/calibfile, 
# #Normalize intensity? (0=no, #=range), Distort?, Interpolate, PlotCalibPoints, booPlotLabels, plotbackfile, 
# #Gaussian) <-- other options
# #p=[ 272.45928478  , 91.5 ,  59.25339245 ,  89.61147631 ,  89.63791054,   89.81178189 ,  68.1669475 ,   63.3555271 ,   1.10221798 ,  31.8848023,  200] 
# #p=[ 272.35233865 ,  91.17544429  , 59.24907075 ,  90.19912466  , 89.63148374,   89.2621423  ,  68.1372743   , 62.95098288   , 1.0992334  ,  32.24603515,  199.65078345]
# #p = [ 271.9473,   91.137 ,  59.2663218 ,   89.65582564 ,  89.62544383,   89.76937348  , 68.19750843 ,  63.29873297  ,  1.65  , 31.92112407,  199.70153725]
# #main(p,SEDMode=0,booPlot=True,specFile='arcturus.txt',booInterpolate=True,plotBackImage='c_arcturus_1min_4.fits')
# #*********************************************************************************************************************
# 
# #Spectral Resolution-CCDMap mercury with labels******************************************************************************************************************
# #(beam phi, beam theta, prism1 phi, prism1 theta, prism2 phi, prism2 theta, grating phi, grating theta, grating alpha,
# #         blaze period (microns), focal length(mm))
# #(SEDMode(0=Max, 1=Random, 2=Sun, 3=from specFile, 4=from CalibFile), Plot?, specFile/calibfile, 
# #Normalize intensity? (0=no, #=range), Distort?, Interpolate, PlotCalibPoints, booPlotLabels, plotbackfile,
# # Gaussian) <-- other options
# #p = [272.31422902, 90.7157937, 59.6543365, 90.21334551, 89.67646101, 89.82098015, 68.0936684,  65.33694031, 1.19265536, 31.50321471, 199.13548823]
# #main(p,SEDMode=4,booPlot=True,specFile='c_noFlat_Hg_0deg_10s.txt',booPlotLabels=True,plotBackImage='c_noFlat_sky_0deg_460_median.fits')
# #*********************************************************************************************************************
# 
# 
# #Spectral Resolution-Gaussian Fit in Orders******************************************************************************************************************
# #(beam phi, beam theta, prism1 phi, prism1 theta, prism2 phi, prism2 theta, grating phi, grating theta, grating alpha,
# #         blaze period (microns), focal length(mm), distortion term)
# #(SEDMode(0=Max, 1=Random, 2=Sun, 3=from specFile, 4=from CalibFile), Plot?, specFile/calibfile, 
# #Normalize intensity? (0=no, #=range), Distort?, Interpolate, PlotCalibPoints, booPlotLabels, plotbackfile,
# # Gaussian) <-- other options
# #p = [272.31422902, 90.7157937, 59.6543365, 90.21334551, 89.67646101, 89.82098015, 68.0936684,  65.33694031, 1.19265536, 31.50321471, 199.13548823]
# #args=['0','0','c_noFlat_Hg_0deg_10s.txt','1','0','1','0','0','c_noFlat_Hg_0deg_10s.fits','1']
# #main(p,SEDMode=0,specFile='c_noFlat_Hg_0deg_10s.txt',booInterpolate=True,plotBackImage='c_noFlat_Hg_0deg_10s.fits',booGaussianFit=True)
# #*********************************************************************************************************************
# 
# 
# #Complete Solar Spectrum***********************************************************************************
# #(SEDMode(0=Max, 1=Random, 2=Sun, 3=from specFile, 4=from CalibFile), Plot?, specFile/calibfile, 
# #Normalize intensity? (0=no, #=range), Distort?, Interpolate, PlotCalibPoints, booPlotLabels, plotbackfile, 
# #Gaussian) <-- other options#main(args=['4','1','solar.txt','1','0','0','0','1','c_noFlat_sky_0deg_460_median.fits'])
# #p = [272.31422902, 90.7157937, 59.6543365, 90.21334551, 89.67646101, 89.82098015, 68.0936684,  65.33694031, 1.19265536, 31.50321471, 199.13548823]
# #main(p,SEDMode=0,specFile='solar.txt',booInterpolate=True,plotBackImage='c_noFlat_sky_0deg_460_median.fits')
# #********************************************************************************************************
# 
# 
# #Spectral Resolution-CCDMap mercury without labels and with calibpoints******************************************************************************************************************
# #(beam phi, beam theta, prism1 phi, prism1 theta, prism2 phi, prism2 theta, grating phi, grating theta, grating alpha,
# #         blaze period (microns), focal length(mm))
# #(SEDMode(0=Max, 1=Random, 2=Sun, 3=from specFile, 4=from CalibFile), Plot?, specFile/calibfile, 
# #Normalize intensity? (0=no, #=range), Distort?, Interpolate, PlotCalibPoints, booPlotLabels, plotbackfile,
# # Gaussian) <-- other options
# #p = [272.31422902, 90.7157937, 59.6543365, 90.21334551, 89.67646101, 89.82098015, 68.0936684,  65.33694031, 1.19265536, 31.50321471, 199.13548823]
# #main(p,SEDMode=4,booPlot=True,specFile='c_noFlat_Hg_0deg_10s.txt',booPlotCalibPoints=True,plotBackImage='c_noFlat_Hg_0deg_10s.fits')
# #*********************************************************************************************************************
# 
# #Plot Errors Example**********************************************************************************
# #p_try=[ 272.45928478  , 91.17527895 ,  59.25339245 ,  89.61147631 ,  89.63791054,   89.81178189 ,  68.1669475 ,   63.3555271 ,   1.10221798 ,  31.8848023,  199.70165672] #best
# ##p_try=[250.6, 64.2, 58.3, 77.9, 89.5, 85.5, 58.5, 61.7, 1.6, 33.7, 193.5] #doeesn't work
# ## ok p_try=[279, 90.6, 59, 90, 89, 90, 70.75, 65, 0, 31.64, 210]
# ##p_try = [278.6, 90.6, 61, 90, 91, 90, 70.75, 63.95, 0, 31.9, 210]
# ##p_try = [278.2,90.8,59.0,90.4,88.9,89.7,70.8,65.1,0.758,31.67,203.1]
# #only h-alpha
# #p_try = [ 272.46778143,   91.14782003 ,  59.2663218 ,   89.65582564 ,  89.62544383,   89.76937348  , 68.19750843 ,  63.29873297  ,  1.10998689  , 31.92112407,  199.70153725]
# #
# ###next
# #p = [272.31422902, 90.7157937, 59.6543365, 90.21334551, 89.67646101, 89.82098015, 68.0936684,  65.33694031, 1.19265536, 31.50321471, 199.13548823]
# #temp , waves = main_errors(p, SEDMode=4,specFile='c_noFlat_Hg_0deg_10s.txt',booPlotCalibPoints=True,PlotBackImage='c_noFlat_Hg_0deg_10s.fits')
# #x2=temp[0]
# #y2=temp[1]
# #
# #d=np.sqrt(x2**2+y2**2)
# #plt.ylabel('Error (Pixels)')
# #plt.xlabel('Emission Line (micrometers)')
# #plt.title('Fitting Error')
# #plt.scatter(waves,d)
# #plt.show()
# #********************************************************************************************************
# 
# #Solar SPectrum CCDMap without labels***********************************************************************************
# #identify fraunhofer
# #SED fraunhofer
# #(SEDMode(0=Max, 1=Random, 2=Sun, 3=from specFile, 4=from CalibFile), Plot?, specFile/calibfile, 
# #Normalize intensity? (0=no, #=range), Distort?, Interpolate, PlotCalibPoints, booPlotLabels, plotbackfile, 
# #Gaussian) <-- other options#main(args=['4','1','solar.txt','1','0','0','0','1','c_noFlat_sky_0deg_460_median.fits'])
# #p = [272.31422902, 90.7157937, 59.6543365, 90.21334551, 89.67646101, 89.82098015, 68.0936684,  65.33694031, 1.19265536, 31.50321471, 199.13548823]
# #main(p,SEDMode=0,booPlot=True,specFile='solar.txt',plotBackImage='c_noFlat_sky_0deg_460_median.fits')
# #********************************************************************************************************
# 
# 
# #Hg full spectrum******************************************************************************************************************
# #(beam phi, beam theta, prism1 phi, prism1 theta, prism2 phi, prism2 theta, grating phi, grating theta, grating alpha,
# #         blaze period (microns), focal length(mm), distortion term)
# #(SEDMode(0=Max, 1=Random, 2=Sun, 3=from specFile, 4=from CalibFile), Plot?, specFile/calibfile, 
# #Normalize intensity? (0=no, #=range), Distort?, Interpolate, PlotCalibPoints, booPlotLabels, plotbackfile,
# # Gaussian) <-- other options
# #p = [272.31422902, 90.7157937, 59.6543365, 90.21334551, 89.67646101, 89.82098015, 68.0936684,  65.33694031, 1.19265536, 31.50321471, 199.13548823]
# #main(p,SEDMode=0,specFile='c_noFlat_Hg_0deg_10s.txt',booInterpolate=True,plotBackImage='c_noFlat_Hg_0deg_10s.fits')
# #*********************************************************************************************************************
