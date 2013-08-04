import wsm
from constants import *
import numpy as np
import xml_parser as xml
import pylab as plt
import optics as o


def launch_task(task):  
    
    if task==1: # Plot solar SED
        a=wsm.do_sed_map(SEDMode = SED_MODE_SOLAR)
        wsm.do_plot_sed_map(a)
        
    elif task==2: #Create and plot CCD map from continuous source 
        a = wsm.do_sed_map(minLambda=0.3, maxLambda=0.9) 
        b = wsm.do_ccd_map(a)
        wsm.do_plot_ccd_map(b)
    
    elif task==3: #Read a calibration fits file
        '''Extract centroids.
        Assign wavelengths and uncertainties.
        Write output file.'''
        calibrationImageFileName = 'hg_rhea_sample1.fits'
        specXMLFileNamefactorTry = 'rhea.xml'
        outputFileName = 'hg_rhea_sample1.txt'
        wsm.do_read_calibration_file(calibrationImageFileName, specXMLFileName, outputFileName, booPlotInitialPoints = False, booPlotFinalPoints = True)
       
    elif task==4: #Find fit
        calibrationDataFileName = 'c_hg_rhea_sample1.txt'
        specXMLFileName = 'rhea_initial.xml'
        calibrationImageFileName = 'hg_rhea_sample1.fits'
        activeCamera = 0
        diagTry = [1.5, 0.5, 0.9 , 0.2, 1.2, 1.2, 1.6, 1, 1.2, 6.7, 50, 50] 
        showStats = False

        SEDMap = wsm.do_sed_map(SEDMode=SED_MODE_FILE, spectrumFileName='hg_spectrum.txt') #create SEDMap from flat Hg emission file
        p_out = wsm.do_find_fit(SEDMap, specXMLFileName, calibrationDataFileName, activeCamera, diagTry = diagTry, showStats = showStats)      
        print p_out
        
        #plot fit to see new model
        CCDMap = wsm.do_ccd_map(SEDMap, specXMLFileName, p_try = p_out[0])
        wsm.do_plot_calibration_points(calibrationImageFileName, calibrationDataFileName, CCDMap, booLabels = True, canvasSize=1, title = 'Calibration points vs Model comparison ')
    
    elif task==5: #Full loop from calibration file to model
        #wsm.do_plot_sed_map(SEDMap, title='Hg SED raw input')
        
        #step 1 - Read fits, extract calibration points to c_+outputFileName (n x 5 np.array) 
        calibrationImageFileName = 'hg_rhea_sample2.fits'
        outputFileName = 'hg_rhea_sample2.txt'
        specXMLFileName =  'rhea_v2.xml'
        wsm.do_read_calibration_file(calibrationImageFileName, specXMLFileName, outputFileName)
        
        #step 2 - find best fit to results extracted
        SEDMap = wsm.do_sed_map(SEDMode=SED_MODE_FILE, specFile='hg_spectrum.txt') #create SEDMap from flat Hg emission file
        p_out = wsm.do_find_fit(SEDMap, specXMLFileName, calibrationDataFileName =  outputFileName)      
        print p_out
        p_out = p_out[0]
    
        #step 3 - plot the optimised solution
        #p_out=np.array([272. , 90.7157937, 59.6543365, 90.21334551, 89.67646101, 89.82098015, 68.0936684, 65.33694031, 1.19265536, 31.50321471, 199.13548823])
        CCDMap = wsm.do_ccd_map(SEDMap, specXMLFileName, p_try = p_out)
        wsm.do_plot_calibration_points(calibrationImageFileName, 'c_' + outputFileName, CCDMap, labels = False, canvasSize=1, title = 'Calibration points vs Model comparison ')
    
    elif task==6: #extract spectrum
        a = wsm.do_sed_map(SEDMode=SED_MODE_CALIB, specFile='Hg_5lines_double.txt')
        wsm.do_plot_sed_map(a)
        b = wsm.do_ccd_map(a)  
        wsm.do_plot_ccd_map(b)
        c = wsm.do_full_extract_order(b , 87, 'c_noFlat_sky_0deg_460_median.fits')
    
    elif task==7: #plot points from calibration file
        calibrationImageFileName = 'hg_rhea_sample1.fits'
        wsm.do_plot_calibration_points(calibrationImageFileName, 'c_hg_rhea_sample1.txt', booLabels = True, canvasSize=1, title = 'Calibration points vs Model comparison ')
    
    elif task==8: #find fit, play with parameter space 
        calibrationDataFileName = 'c_hg_rhea_sample1.txt'
        specXMLFileName = 'rhea_initial.xml'
        calibrationImageFileName = 'hg_rhea_sample1.fits'
        activeCamera = 0
        showStats = True
        
        iter = np.arange(1)
        dist = np.zeros(len(iter))
        avgDist = np.zeros(len(iter))
        
        SEDMap = wsm.do_sed_map(SEDMode=SED_MODE_FILE, spectrumFileName='hg_spectrum.txt') #create SEDMap from flat Hg emission file

        plt.ion()
        for i in np.arange(len(iter)):
            diagTry = [1.5, 0.5, 0.9 , 0.2, 1.2, 1.2, 1.6, 1, 1.2, 67, 1, 1]  
            factorTry = 1
            p_out = wsm.do_find_fit(SEDMap, specXMLFileName, calibrationDataFileName, activeCamera, diagTry = diagTry, factorTry = factorTry, showStats = showStats)      
            fvec = p_out[2]['fvec']
            x = fvec[:len(fvec)/2]
            y = fvec[len(fvec)/2:]
            dist[i] = sum(np.sqrt(x**2+y**2))
            avgDist[i] = dist[i]/len(x)
            
            plt.clf()
#            plt.scatter(iter[:i], dist[:i])
            plt.scatter(iter[:i+1], avgDist[:i+1], color = 'red')
            plt.grid()
            plt.draw()           
            
        #plot fit to see new model
#        CCDMap = wsm.do_ccd_map(SEDMap, specXMLFileName, p_try = p_out)
#        wsm.do_plot_calibration_points(calibrationImageFileName, calibrationDataFileName, CCDMap, booLabels = True, canvasSize=1, title = 'Calibration points vs Model comparison ')
        plt.ioff()
        plt.show()
        
    elif task==100: #random stuff
        a=wsm.do_sed_map(SEDMode=SED_MODE_FILE, specFile = 'hg_spectrum.txt')
    #    diag_try = [10,1,1,1,1,1,1,1,1,1,1]
        diag_try = []
        b = wsm.do_find_fit(a, 'rhea_initial.xml', 'c_hg_rhea_sample1.txt', factor_try = 1.5, diag_try = diag_try)
        print b
        print 'Finished'

    elif task==101: #Distortion tests
        inX = np.arange(-50,51)
        inX = np.tile(inX,len(inX))
        inY = np.sort(inX)
        K=-0.3
        
        outX, outY = o.distort(inX, inY, K, 20, 20)
        
        plt.scatter(outX, outY)
        plt.show()

    elif task==102: #test filename maker
        import wsmtools as wt
        a = wt.find_specXMLFileName('rhea.xml')
        
launch_task(102)


#a=wsm.do_sed_map(minLambda=0.3, maxLambda=0.9)
#b=wsm.do_ccd_map(a)
#c=wsm.do_full_extract_order(b , 87, 'c_noFlat_sky_0deg_460_median.fits')
#
#
#
#
#
#print 'do stuff now'
#wsm.do_plot_ccd_map(b)


#SEDMap=wsm.do_sed_map(SEDMode=wsm.SED_MODE_CALIB, specFile='calib_out.txt')
#p=wsm.do_find_fit(SEDMap, 'calib_out.txt')
#print p
#import image_calibration as ic

#
#im = wt.load_image_map()
#print im

#ic.analyseImage('test.fits', 'a.txt')













#
#a=wsm.do_sed_map(SEDMode=wsm.SED_MODE_CALIB, specFile='Hg_5lines_double.txt')
# 
#b=wsm.do_ccd_map(a)
#print b
#wsm.do_plot_ccd_map(b)




#from scipy import stats
#import numpy as np
#import ds9
# 
## Make a 2D gaussian image that is stored in a 2D numpy array
#x = np.arange(-3, 3, 0.1)
#xx, yy = np.meshgrid(x, x)
#gauss2d = stats.norm.pdf(xx) * stats.norm.pdf(yy)
# 
## Now open ds9 (this assumes no ds9 instance is yet running)
#d = ds9.ds9()
# 
## Load up our 2D gaussian
#d.set_np2arr(gauss2d)
## ~/Documents/workspace/rhea-wsm/fits
## Zoom to fit
#d.set('zoom to fit')
# 
## Change the colormap and scaling
#d.set('cmap bb')
#d.set('scale log')
# 
## Add a label
#d.set('regions command {text 30 20 #text="Fun with pyds9" font="times 18 bold"}')
# 
## Now you can play in ds9 to your heart's content.
## Check back to see what the current color scale is.
#print d.get('scale')
# 
## Finally, save your completed image (including regions or labels)
#d.set('saveimage png my_pyds9_img.png')



#print "3"
#import xml_parser
#
#
#p=xml_parser.read_All()
#
#print p

#os.chdir('/media/win-desktop/6hrs/Darks')
#inFiles = ['Hg_6hrs_20s_D%01d.FIT' % k for k in range(1,6)]
#outFilename = 'masterDark.fits' 
#medianCombine(inFiles, outFilename)

#
#fit,success=findFit('c_noFlat_Hg_0deg_10s.txt')
#print fit
#main(p=fit)

#main()

#compareTemp('6hrs/Medians/med1_5.fits','6hrs/Medians/med711_715.fits','6hrs/Results/out1.fits')

#batchAverage2()

#batchMedian()

#inFiles = ['sky_0deg_460_%01d.fits' % k for k in range(1,4)]
#medianCombine(inFiles, 'sky_0deg_460_median.fits')
#inFiles = ['sky_0deg_1550_%01d.fits' % k for k in range(1,4)]
#medianCombine(inFiles, 'sky_0deg_1550_median.fits')

#inFiles = ['sky_0deg_460_%01d.fits' % k for k in range(1,4)]
#averageCombine(inFiles, 'sky_0deg_460_average.fits')
#inFiles = ['sky_0deg_1550_%01d.fits' % k for k in range(1,4)]
#averageCombine(inFiles, 'sky_0deg_1550_average.fits')

#calibrateImage('sky_0deg_460_median.fits','mbias.fits','mdark.fits','mflat.fits','c_noFlat_sky_0deg_460_median.fits')
#calibrateImage('sky_0deg_1550_median.fits','mbias.fits','mdark.fits','mflat.fits','c_noFlat_sky_0deg_1550_median.fits')

#calibrateImage('Hg_0deg_10s.fits','mbias.fits','mdark.fits','mflat.fits','c_Hg_0deg_10s.fits')

#calibrateImage('sky_0deg_460_1.fits','mbias.fits','mdark.fits','mflat.fits','c_noFlat_sky_0deg_460_1.fits')
#calibrateImage('sky_0deg_460_2.fits','mbias.fits','mdark.fits','mflat.fits','c_noFlat_sky_0deg_460_2.fits')
#calibrateImage('sky_0deg_460_3.fits','mbias.fits','mdark.fits','mflat.fits','c_noFlat_sky_0deg_460_3.fits')

#calibrateImage('sky_0deg_1550_median.fits','mbias.fits','mdark.fits','mflat.fits','c_noFlat_sky_0deg_1550_median.fits')




#inFiles = ['dark_a%01d.fits' % k for k in range(1,4)]
#medianCombine(inFiles, 'dark_a_median.fits')
#
#inFiles = ['dark_a%01d.fits' % k for k in range(1,4)]
#averageCombine(inFiles,  'dark_a_average.fits')
#
#
#subtractDark('arcturus_1min_1.fits','dark_a_average.fits','c_arcturus_1min_1.fits')
#subtractDark('arcturus_1min_2.fits','dark_a_average.fits','c_arcturus_1min_2.fits')
#subtractDark('arcturus_1min_3.fits','dark_a_average.fits','c_arcturus_1min_3.fits')
#subtractDark('arcturus_1min_4.fits','dark_a_average.f#inFiles = ['dark_a%01d.fits' % k for k in range(1,4)]
#medianCombine(inFiles, 'dark_a_median.fits')
#
#inFiles = ['dark_a%01d.fits' % k for k in range(1,4)]
#averageCombine(inFiles,  'dark_a_average.fits')
#
#
#subtractDark('arcturus_1min_1.fits','dark_a_average.fits','c_arcturus_1min_1.fits')
#subtractDark('arcturus_1min_2.fits','dark_a_average.fits','c_arcturus_1min_2.fits')
#subtractDark('arcturus_1min_3.fits','dark_a_average.fits','c_arcturus_1min_3.fits')
#subtractDark('arcturus_1min_4.fits','dark_a_average.fits','c_arcturus_1min_4.fits')
#
#inFiles = ['c_arcturus_1min_%01d.fits' % k for k in range(1,5)]
#medianCombine(inFiles, 'c_arcturus_1min_median.fits')
#
#inFiles = ['c_arcturus_1min_%01d.fits' % k for k in range(1,5)]
#averageCombine(inFiles, 'c_arcturus_1min_average.fits')its','c_arcturus_1min_4.fits')
#
#inFiles = ['c_arcturus_1min_%01d.fits' % k for k in range(1,5)]
#medianCombine(inFiles, 'c_arcturus_1min_median.fits')
#
#inFiles = ['c_arcturus_1min_%01d.fits' % k for k in range(1,5)]
#averageCombine(inFiles, 'c_arcturus_1min_average.fits')


#print x(0.502)

#p=main_errors(p = [ 271.92998622,   91.03999719,   59.48997316,   89.83000496,   89.37002499,   89.79999531,   68.03002976,   64.9399939,     1.15998754,   31.52736851,  200.00000005], SEDMode=4,booPlot=True,specFile='c_noFlat_Hg_0deg_10s.txt')

#main(p = [279, 90.6, 59, 90, 89, 90, 70.75, 65, 0, 31.64, 210])
#main(p = [ 271.92998622,   91.03999719,   59.48997316,   89.83000496,   89.37002499,   89.79999531,   68.03002976,   64.9399939,     1.15998754,   31.52736851,  200.00000005]) ####best fit so far....
##
##print arrayOut
##
#xi = np.arange(0,9)
#A = np.array([ xi, np.ones(9)])
## linearly generated sequence
#y = [19, 20, 20.5, 21.5, 22, 23, 23, 25.5, 24]
#w = np.linalg.lstsq(A.T,y)[0] # obtaining the parameters
#s=1


#best p for previous fitting:
#p = [ 271.92998622,   91.03999719,   59.48997316,   89.83000496,   89.37002499,   89.79999531,   68.03002976,   64.9399939,     1.15998754,   31.52736851,  200.00000005]
#'c_noFlat_Hg_0deg_10s.fits')#,'solar.png')#,'mdark.fits')#,'simple_flat.fits')#


#Fraunhofer Lines-CCDMap with labels***********************************************************************************
#identify fraunhofer
#SED fraunhofer
#(SEDMode(0=Max, 1=Random, 2=Sun, 3=from specFile, 4=from CalibFile), Plot?, specFile/calibfile, 
#Normalize intensity? (0=no, #=range), Distort?, Interpolate, PlotCalibPoints, booPlotLabels, plotbackfile, 
#Gaussian) <-- other options#main(args=['4','1','solar.txt','1','0','0','0','1','c_noFlat_sky_0deg_460_median.fits'])
#p = [272.31422902, 90.7157937, 59.6543365, 90.21334551, 89.67646101, 89.82098015, 68.0936684,  65.33694031, 1.19265536, 31.50321471, 199.13548823]
#main(p,SEDMode=4,booPlot=True,specFile='solar.txt')
#********************************************************************************************************


#Arcturus callibration**********************************************************************************
#p_try=[ 272.45928478  , 91.17527895 ,  59.25339245 ,  89.61147631 ,  89.63791054,   89.81178189 ,  68.1669475 ,   63.3555271 ,   1.10221798 ,  31.8848023,  199.70165672] #best
##p_try=[250.6, 64.2, 58.3, 77.9, 89.5, 85.5, 58.5, 61.7, 1.6, 33.7, 193.5] #doeesn't work
## ok p_try=[279, 90.6, 59, 90, 89, 90, 70.75, 65, 0, 31.64, 210]
##p_try = [278.6, 90.6, 61, 90, 91, 90, 70.75, 63.95, 0, 31.9, 210]
##p_try = [278.2,90.8,59.0,90.4,88.9,89.7,70.8,65.1,0.758,31.67,203.1]
#only h-alpha
#p_try = [ 272.46778143,   91.14782003 ,  59.2663218 ,   89.65582564 ,  89.62544383,   89.76937348  , 68.19750843 ,  63.29873297  ,  1.10998689  , 31.92112407,  199.70153725]
#
##next
#p_try = [ 272.45928478,   91.17527895  , 59.25339245 ,  89.61147631 ,  89.63791054,  89.81178189 ,  68.1669475 ,   63.3555271,     1.10221798 ,  31.8848023,  199.70165672]
#
#diag_try=[1,1,1,1,1,1,1,1,1,0.1,1]
#factor_try=0.3
#fit, success=findFit('arcturus.txt',p_try, factor_try, diag_try)
#print fit
#main(fit,SEDMode=4,booPlot=True,specFile='arcturus.txt',booPlotCalibPoints=True,plotBackImage='c_arcturus_1min_4.fits')
#********************************************************************************************************


#Arcturus plot********************************************************************************************************
#(beam phi, beam theta, prism1 phi, prism1 theta, prism2 phi, prism2 theta, grating phi, grating theta, grating alpha,
#         blaze period (microns), focal length(mm), distortion term)
#(SEDMode(0=Max, 1=Random, 2=Sun, 3=from specFile, 4=from CalibFile), Plot?, specFile/calibfile, 
#Normalize intensity? (0=no, #=range), Distort?, Interpolate, PlotCalibPoints, booPlotLabels, plotbackfile, 
#Gaussian) <-- other options
#p=[ 272.45928478  , 91.5 ,  59.25339245 ,  89.61147631 ,  89.63791054,   89.81178189 ,  68.1669475 ,   63.3555271 ,   1.10221798 ,  31.8848023,  200] 
#p=[ 272.35233865 ,  91.17544429  , 59.24907075 ,  90.19912466  , 89.63148374,   89.2621423  ,  68.1372743   , 62.95098288   , 1.0992334  ,  32.24603515,  199.65078345]
#p = [ 271.9473,   91.137 ,  59.2663218 ,   89.65582564 ,  89.62544383,   89.76937348  , 68.19750843 ,  63.29873297  ,  1.65  , 31.92112407,  199.70153725]
#main(p,SEDMode=0,booPlot=True,specFile='arcturus.txt',booInterpolate=True,plotBackImage='c_arcturus_1min_4.fits')
#*********************************************************************************************************************

#Spectral Resolution-CCDMap mercury with labels******************************************************************************************************************
#(beam phi, beam theta, prism1 phi, prism1 theta, prism2 phi, prism2 theta, grating phi, grating theta, grating alpha,
#         blaze period (microns), focal length(mm))
#(SEDMode(0=Max, 1=Random, 2=Sun, 3=from specFile, 4=from CalibFile), Plot?, specFile/calibfile, 
#Normalize intensity? (0=no, #=range), Distort?, Interpolate, PlotCalibPoints, booPlotLabels, plotbackfile,
# Gaussian) <-- other options
#p = [272.31422902, 90.7157937, 59.6543365, 90.21334551, 89.67646101, 89.82098015, 68.0936684,  65.33694031, 1.19265536, 31.50321471, 199.13548823]
#main(p,SEDMode=4,booPlot=True,specFile='c_noFlat_Hg_0deg_10s.txt',booPlotLabels=True,plotBackImage='c_noFlat_sky_0deg_460_median.fits')
#*********************************************************************************************************************


#Spectral Resolution-Gaussian Fit in Orders******************************************************************************************************************
#(beam phi, beam theta, prism1 phi, prism1 theta, prism2 phi, prism2 theta, grating phi, grating theta, grating alpha,
#         blaze period (microns), focal length(mm), distortion term)
#(SEDMode(0=Max, 1=Random, 2=Sun, 3=from specFile, 4=from CalibFile), Plot?, specFile/calibfile, 
#Normalize intensity? (0=no, #=range), Distort?, Interpolate, PlotCalibPoints, booPlotLabels, plotbackfile,
# Gaussian) <-- other options
#p = [272.31422902, 90.7157937, 59.6543365, 90.21334551, 89.67646101, 89.82098015, 68.0936684,  65.33694031, 1.19265536, 31.50321471, 199.13548823]
#args=['0','0','c_noFlat_Hg_0deg_10s.txt','1','0','1','0','0','c_noFlat_Hg_0deg_10s.fits','1']
#main(p,SEDMode=0,specFile='c_noFlat_Hg_0deg_10s.txt',booInterpolate=True,plotBackImage='c_noFlat_Hg_0deg_10s.fits',booGaussianFit=True)
#*********************************************************************************************************************


#Complete Solar Spectrum***********************************************************************************
#(SEDMode(0=Max, 1=Random, 2=Sun, 3=from specFile, 4=from CalibFile), Plot?, specFile/calibfile, 
#Normalize intensity? (0=no, #=range), Distort?, Interpolate, PlotCalibPoints, booPlotLabels, plotbackfile, 
#Gaussian) <-- other options#main(args=['4','1','solar.txt','1','0','0','0','1','c_noFlat_sky_0deg_460_median.fits'])
#p = [272.31422902, 90.7157937, 59.6543365, 90.21334551, 89.67646101, 89.82098015, 68.0936684,  65.33694031, 1.19265536, 31.50321471, 199.13548823]
#main(p,SEDMode=0,specFile='solar.txt',booInterpolate=True,plotBackImage='c_noFlat_sky_0deg_460_median.fits')
#********************************************************************************************************


#Spectral Resolution-CCDMap mercury without labels and with calibpoints******************************************************************************************************************
#(beam phi, beam theta, prism1 phi, prism1 theta, prism2 phi, prism2 theta, grating phi, grating theta, grating alpha,
#         blaze period (microns), focal length(mm))
#(SEDMode(0=Max, 1=Random, 2=Sun, 3=from specFile, 4=from CalibFile), Plot?, specFile/calibfile, 
#Normalize intensity? (0=no, #=range), Distort?, Interpolate, PlotCalibPoints, booPlotLabels, plotbackfile,
# Gaussian) <-- other options
#p = [272.31422902, 90.7157937, 59.6543365, 90.21334551, 89.67646101, 89.82098015, 68.0936684,  65.33694031, 1.19265536, 31.50321471, 199.13548823]
#main(p,SEDMode=4,booPlot=True,specFile='c_noFlat_Hg_0deg_10s.txt',booPlotCalibPoints=True,plotBackImage='c_noFlat_Hg_0deg_10s.fits')
#*********************************************************************************************************************

#Plot Errors Example**********************************************************************************
#p_try=[ 272.45928478  , 91.17527895 ,  59.25339245 ,  89.61147631 ,  89.63791054,   89.81178189 ,  68.1669475 ,   63.3555271 ,   1.10221798 ,  31.8848023,  199.70165672] #best
##p_try=[250.6, 64.2, 58.3, 77.9, 89.5, 85.5, 58.5, 61.7, 1.6, 33.7, 193.5] #doeesn't work
## ok p_try=[279, 90.6, 59, 90, 89, 90, 70.75, 65, 0, 31.64, 210]
##p_try = [278.6, 90.6, 61, 90, 91, 90, 70.75, 63.95, 0, 31.9, 210]
##p_try = [278.2,90.8,59.0,90.4,88.9,89.7,70.8,65.1,0.758,31.67,203.1]
#only h-alpha
#p_try = [ 272.46778143,   91.14782003 ,  59.2663218 ,   89.65582564 ,  89.62544383,   89.76937348  , 68.19750843 ,  63.29873297  ,  1.10998689  , 31.92112407,  199.70153725]
#
###next
#p = [272.31422902, 90.7157937, 59.6543365, 90.21334551, 89.67646101, 89.82098015, 68.0936684,  65.33694031, 1.19265536, 31.50321471, 199.13548823]
#temp , waves = main_errors(p, SEDMode=4,specFile='c_noFlat_Hg_0deg_10s.txt',booPlotCalibPoints=True,PlotBackImage='c_noFlat_Hg_0deg_10s.fits')
#x2=temp[0]
#y2=temp[1]
#
#d=np.sqrt(x2**2+y2**2)
#plt.ylabel('Error (Pixels)')
#plt.xlabel('Emission Line (micrometers)')
#plt.title('Fitting Error')
#plt.scatter(waves,d)
#plt.show()
#********************************************************************************************************

#Solar SPectrum CCDMap without labels***********************************************************************************
#identify fraunhofer
#SED fraunhofer
#(SEDMode(0=Max, 1=Random, 2=Sun, 3=from specFile, 4=from CalibFile), Plot?, specFile/calibfile, 
#Normalize intensity? (0=no, #=range), Distort?, Interpolate, PlotCalibPoints, booPlotLabels, plotbackfile, 
#Gaussian) <-- other options#main(args=['4','1','solar.txt','1','0','0','0','1','c_noFlat_sky_0deg_460_median.fits'])
#p = [272.31422902, 90.7157937, 59.6543365, 90.21334551, 89.67646101, 89.82098015, 68.0936684,  65.33694031, 1.19265536, 31.50321471, 199.13548823]
#main(p,SEDMode=0,booPlot=True,specFile='solar.txt',plotBackImage='c_noFlat_sky_0deg_460_median.fits')
#********************************************************************************************************


#Hg full spectrum******************************************************************************************************************
#(beam phi, beam theta, prism1 phi, prism1 theta, prism2 phi, prism2 theta, grating phi, grating theta, grating alpha,
#         blaze period (microns), focal length(mm), distortion term)
#(SEDMode(0=Max, 1=Random, 2=Sun, 3=from specFile, 4=from CalibFile), Plot?, specFile/calibfile, 
#Normalize intensity? (0=no, #=range), Distort?, Interpolate, PlotCalibPoints, booPlotLabels, plotbackfile,
# Gaussian) <-- other options
#p = [272.31422902, 90.7157937, 59.6543365, 90.21334551, 89.67646101, 89.82098015, 68.0936684,  65.33694031, 1.19265536, 31.50321471, 199.13548823]
#main(p,SEDMode=0,specFile='c_noFlat_Hg_0deg_10s.txt',booInterpolate=True,plotBackImage='c_noFlat_Hg_0deg_10s.fits')
#*********************************************************************************************************************
