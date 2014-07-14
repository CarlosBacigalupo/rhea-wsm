#external libs
import numpy as np
import pyfits as pf
import pylab as plt
from scipy.optimize import leastsq
from scipy import interpolate
import os

#internal libs
import wsmtools as wt


class spectrograph():
    
    name = 'rhea'
    arcMapMask = ''
    
    def __init__(self):
        self.camera = cameras()
        self.optics = optics()
        self.beams = beams()
        self.orders = np.array([34,35,36,37])

        
    def create_sed_map(self, SEDMode=wt.c.sedMode.flat, minLambda=0.4, maxLambda=0.78, deltaLambda=0.0001, intNormalize=0, spectrumFileName='', minI=0, maxI=1e9): 
        '''
        Creates/Loads the input Spectrum Energy Density map. It simulates the characteristics of the input beam. 
        
        Parameters
        ----------
        SEDMode : int
            Mode for the creation of the SEDMap
            0=Flat, 1=Random, 2=Sun, 3=from specFile, 4=from Calibration file
            
        minLambda : np.float32
            Lower limit for SEDMap.
    
        maxLambda : np.float32
            Higher limit for SEDMap.
    
        deltaLambda : np.float32
            Step between wavelengths.
        
        intNormalize : integer
            If !=0, it normalizes to intNormalize value
            
        specFile : str
            File with the SED of the source (used when SEDmode = (SED_MODE_CALIB or SED_MODE_FILE))
            
        Returns
        -------
        SEDMap : np.array
            n x 2 np.array with wavelength, Energy
          
        Notes
        -----
        '''  
        if SEDMode==wt.c.sedMode.flat: #Flat
            SEDMap = np.array((np.arange(minLambda, maxLambda + deltaLambda, deltaLambda),np.ones(np.arange(minLambda, maxLambda + deltaLambda, deltaLambda).size)))
    
        elif SEDMode==wt.c.sedMode.random: #Random        
            np.hstack((range(minLambda, maxLambda + deltaLambda, deltaLambda),[random.random() for _ in range(10)]))    
            SEDMap = np.array([0,0])
            for Lambda in range(minLambda, maxLambda + deltaLambda, deltaLambda):
    #            Intensity=int()
    #            newItem=np.np.array([Lambda,random.random(0.0,1.0)])
                SEDMap = np.vstack((SEDMap,np.array([Lambda,random.random(0.0,1.0)])))     
            SEDMap = SEDMap[1:,]
                     
        elif SEDMode==wt.c.sedMode.solar: #Solar   
            sol = astSED.SOL        
            tempA=sol.wavelength.transpose()*1e-4
            tempB=sol.flux.transpose()            
            SEDMap = np.array([tempA, tempB])    
            SEDMap = SEDMap.transpose()
                                      
        elif SEDMode==wt.c.sedMode.file: #From line atlas
            SEDMap = np.loadtxt(spectrumFileName, unpack=True)
            
        elif SEDMode==wt.c.sedMode.calib: #From calibration file
            a=np.loadtxt(TEMP_PATH + spectrumFileName, unpack=True)
            SEDMap=np.array((a[2],np.ones(len(a[2]))))
    #        SEDMap.transpose()
                    
                    
        #Remove rows outside the wavelength range
        SEDMap = np.array([SEDMap[0][SEDMap[0]>=minLambda],SEDMap[1][SEDMap[0]>=minLambda]])     
        SEDMap = np.array([SEDMap[0][SEDMap[0]<=maxLambda],SEDMap[1][SEDMap[0]<=maxLambda]])     
                    
                    
        #Normalize the intensity  
        if intNormalize!=0:    
            fluxRange=(max(SEDMap[wt.c.sedMap.intensity])-min(SEDMap[wt.c.sedMap.intensity]))
            
            if fluxRange==0:
                SEDMap = np.array((SEDMap[wt.c.sedMap.wavelength], np.ones(SEDMap[wt.c.sedMap.wavelength].size)))  
            else:
                SEDMap = np.array((SEDMap[wt.c.sedMap.wavelength], (SEDMap[wt.c.sedMap.intensity]-min(SEDMap[wt.c.sedMap.intensity]))/(fluxRange+1) ))
    
        #Remove rows outside the intensity range
        SEDMap = np.array([SEDMap[0][SEDMap[wt.c.sedMap.intensity]>=minI],SEDMap[1][SEDMap[wt.c.sedMap.intensity]>=minI]])     
        SEDMap = np.array([SEDMap[0][SEDMap[wt.c.sedMap.intensity]<=maxI],SEDMap[1][SEDMap[wt.c.sedMap.intensity]<=maxI]])     
           
        return SEDMap.transpose()

    
    def clear_calib_from_CCDMap(self, calibrationMap, CCDMap, distMin =10):
        '''Removes points in calibration map that are further than distMin from a CCDMap point.
        '''
        
        newCalibrationMap = np.array([[]])
        
        for i in range(len(CCDMap[:,wt.c.CCDMap.x])):
            distance_array = np.sqrt((CCDMap[i,wt.c.CCDMap.x]-calibrationMap[:,wt.c.calibrationMap.x])**2+(CCDMap[i,wt.c.CCDMap.y]-calibrationMap[:,wt.c.calibrationMap.y])**2)
            closest_point = np.min(distance_array)
            if closest_point<=distMin:
                closest_point_index = np.where(distance_array==closest_point)[0][0]       
                if newCalibrationMap.shape[1]==0:
                    newCalibrationMap = np.array([calibrationMap[closest_point_index,:]])
                else:
                    newCalibrationMap = np.append(newCalibrationMap, [calibrationMap[closest_point_index,:]], axis=0)

        return newCalibrationMap

                
    def read_arc_file(self, arcFile, sexParamFile):
        '''
        Extracts peaks with sextractor
        Imports found points into arrays
           
        Parameters
        ----------
        arcFile : string
            fits file with arc light
            
        sexParamFile : string
            sextractor input parameters
    
    
        Returns
        -------
         
          
        Notes
        -----
        '''
        
        out, err_result, sexOutputFileName = wt.ia.analyse_image_sex(arcFile, sexParamFile)
        
        #Loads coordinate points from calibration sextractor output file
        arcFileMap = wt.ia.read_image_map_sex(sexOutputFileName)
        if arcFileMap.shape[0] == 0: 
            print 'sextractor: Point detection returned empty. Check .sex file in base_dir'
            return
        
        # shifts to center of the chip (COC) coords    
        # sextractor is 1 based array
        arcFileMap[:,wt.c.calibrationMap.x] -= 1 
        arcFileMap[:,wt.c.calibrationMap.y] -= 1
            
        #Create initial output file from calibration data
        calibFileNoWl =  'calNoWl_' + ''.join(arcFile).split('.')[0]+ '.txt'
        np.savetxt(calibFileNoWl, arcFileMap, delimiter=' ', newline = '\n')
        
        return calibFileNoWl


    def read_xml(self, modelXMLFile, fitsFile='', p_try = []):        
        #reads spectrograph xml file into corresponding classes
        #fitsfile is any fits file of the right size to verify xml file size parameters
        global a, b
        
        Beams, Optics, Cameras, p, stheta = wt.xml.read_all(modelXMLFile, p_try = p_try)
        
        self.beams.raw = Beams
#         print 'Beam:',Beams
        self.optics.raw = Optics
        self.optics.coords1x = Optics[:,wt.c.optics.coords1x]
        self.optics.coords1y = Optics[:,wt.c.optics.coords1y]
        self.optics.coords1z = Optics[:,wt.c.optics.coords1z]
        self.optics.coords2x = Optics[:,wt.c.optics.coords2x]
        self.optics.coords2y = Optics[:,wt.c.optics.coords2y]
        self.optics.coords2z = Optics[:,wt.c.optics.coords2z]
        self.optics.type = Optics[:,wt.c.optics.type]
        self.optics.n = Optics[:,wt.c.optics.n]
        self.optics.gPeriod = Optics[:,wt.c.optics.gPeriod]
        self.optics.gBlAngle = Optics[:,wt.c.optics.gBlAngle]
#         print 'Optics:',Optics
        
        self.camera.raw = Cameras
        self.camera.name = Cameras[0][wt.c.cameras.name]
        self.camera.width = int(Cameras[0][wt.c.cameras.width])
        self.camera.height = int(Cameras[0][wt.c.cameras.height])
        self.camera.offsetX = (int(Cameras[0][wt.c.cameras.width])-1)/2.
        self.camera.offsetY = (int(Cameras[0][wt.c.cameras.height])-1)/2.
        self.camera.pSize = float(Cameras[0][wt.c.cameras.pSize])
        self.camera.fLength = float(Cameras[0][wt.c.cameras.fLength])
        self.stheta = stheta
        self.p = p
        
        if fitsFile!='':
            hdulist = pf.open(fitsFile)
            imWidth = hdulist[0].header['NAXIS1']
            imHeight = hdulist[0].header['NAXIS2']
            if imWidth!=self.camera.width: print'File width inconsistent between xml data and fitsfile', self.camera.width, imWidth
            if imHeight!=self.camera.height: print'File height inconsistent between xml data and fitsfile', self.camera.height, imHeight

    
    def offset_CCDMap(self):
        #shifts the coords of spectrograph.CCDMap by camera.offset... amount into spectrograph.CCDMapOffset
        if self.CCDMap.shape[1]>0:
        
            self.CCDMapOffset = self.CCDMap.copy()
            self.CCDMapOffset[:,wt.c.CCDMap.x] += self.camera.offsetX
            self.CCDMapOffset[:,wt.c.CCDMap.y] += self.camera.offsetY
        else:
            print 'Found empty CCDMap when trying to apply offset'

        
    def assign_initial_wavelength(self, calibFileNoWl, modelXMLFile, SEDMap, booAvgAdjust = False):
        '''
        Assigns initial wavelength based on proximity
        Optional correct by average
        
           
        Parameters
        ----------
        calibFileNoWl : string
             
    
        modelXMLFile : string
            
            
        SEDMap : n x 2 np.array
            
            
        booPlotPoints : boolean
    
    
        Returns
        -------
        nottin
          
        Notes
        -----
        '''      
        
        calibFile = 'cal_' + ''.join(''.join(calibFileNoWl).split('_')[-1]).split('.')[0] + '.txt'
        
        #Loads coordinate points from calibration initial (no Wl) output file
        calibrationMap = wt.ia.read_full_calibration_data(calibFileNoWl)
        if calibrationMap.shape[0] == 0:
            print 'Calibration data not found in ',  calibFileNoWl
            return

        #Create the model based on initial parameters        
        CCDMap, CCDMapOffset = self.create_ccd_map(SEDMap, modelXMLFile, resetArcMap = False)  
        
        #Find wavelength of detected points (first approximation)
        calibrationMap = wt.ia.identify_imageMapWavelength_avg(CCDMapOffset[self.arcMapMask], calibrationMap, booAvgAdjust = booAvgAdjust)
        
        #Create calibFile with all calibration data
        np.savetxt(calibFile, calibrationMap, delimiter=' ', newline = '\n')
    
        return calibFile

    
    def create_ccd_map(self, SEDMap ,modelXMLFile, resetArcMap = True, p_try = [], booAddCero=False):
        '''
        Computes the projection of n beams of monochromatic light passing through an optical system. 
    
        Parameters
        ----------
        SEDMap : np.array
            n x 2 np.array [wavelength, Energy]
            
        modelXMLFile : str
            name of model xml file
            
        Returns
        -------
        CCDMap : n x 6 np.array [CCDX, CCDY, CCDWavelength, CCDIntensity, CCDOrder, CCDBeam ]
    
        Notes
        -----  
        '''   
                
        #Reads xml file
        self.read_xml(modelXMLFile, p_try=p_try)   
        
        #hack for RHEA. Needs manual reverse of prism on beam return. todo
    #     if modelXMLFile[-8:]=='rhea.xml':
        self.optics.raw[4,0:3]=-self.optics.raw[0,0:3]
        self.optics.raw[3,0:3]=-self.optics.raw[1,0:3]  
        
        self.camera.fLength = self.p[10]
        #Launch grid loop. Creates an array of (x,y,lambda, Intensity, Order, beamID)
        CCDMap, CCDMap3D = wt.optics.ccd_loop(SEDMap, self.beams.raw , self.optics.raw, self.camera , self.stheta, self.orders)
        
        #Distort if any distortion data present
#         K = [self.p[11], self.p[12], self.p[13]]
#         Xc = self.p[14]
#         Yc = self.p[15]
        if CCDMap.shape[1]>0:
#             CCDMap[:,wt.c.CCDMap.x], CCDMap[:,wt.c.CCDMap.y] = wt.optics.distort(CCDMap[:,wt.c.CCDMap.x], CCDMap[:,wt.c.CCDMap.y], K, Xc, Yc)
        
            self.CCDMap = CCDMap
            if booAddCero: self.CCDMap = np.append(self.CCDMap,[[0,0,0.7,1,0,0]],0)
            self.CCDMap3D = CCDMap3D
            self.offset_CCDMap()
            if resetArcMap: self.arcMapMask = np.ones(CCDMap.shape[0]).astype(bool)
        
        return CCDMap, self.CCDMapOffset

    
    def save_arcMaskMap(self):
        np.savetxt(self.name + '_arcMask.txt',self.arcMapMask)
        
 
    def find_fit(self, SEDMap, modelXMLFile, calibrationMap, p_try = [], factorTry=1 ,diagTry = [], showStats = False, maxfev = 1000, booWriteP = True):
        '''
        Wrapper for reading the calibration file, and launching the fitting function
           
        Parameters
        ----------
        calibrationFile : string
            Name of the file with the data from the spectrograph
    
        Returns
        -------
        fit : np np.array
            1 x 12 np np.array with fitted arguments (p np.array)
          
        Notes
        -----
        '''  
        #old mainArgs=['4','0',calibrationFile,'0','0']
        #x,y, wavelist are the positions of the peaks in calibrationFile.
        #x,y,waveList,xsig,ysig = readCalibrationData(calibrationFile)
        #fit is the output, which is the ideal p vector.
        
        if p_try==[]: p_try = wt.xml.read_p(modelXMLFile)
        if diagTry==[]: diagTry = np.ones(len(p_try))
        
        while True:
            fit = leastsq(self.fit_errors, p_try, args=[SEDMap, modelXMLFile, calibrationMap], full_output=True, factor = factorTry, diag = diagTry, maxfev = maxfev)
            if fit[-1] != 0: #workaround to avoid inconsistent message 'wrong input parameters(error code 0)'
                break
            
        if showStats: self.fitting_stats(fit)
        
        if booWriteP: 
            wt.xml.write_p(fit[0], modelXMLFile)
        return fit


    def fit_errors(self, p, args):
        '''
        do_ccd_map will return these vectors in a random order. 
        We assume that there are no rounding errors (probably a hack?)
        and will use a floating point == to identify the (x,y) corresponding
        to each wavelength input.
        NB A much better idea would be to re-write main without vstack but
        instead with a neatly created output structure that is defined at the start.
        i.e. when going from SEDMap to SEDMapLoop the information on which line in
        the input each wavelength came from was lost.
        '''

        SEDMap = args[0]
        modelXMLFile = args[1]
        calibrationMap = args[2]   
    
        CCDMap, CCDMapOffset = self.create_ccd_map(SEDMap, modelXMLFile, resetArcMap = False, p_try = p)
        
        CCDX_m = CCDMapOffset[:,wt.c.CCDMap.x]#[self.arcMapMask]
        CCDY_m = CCDMapOffset[:,wt.c.CCDMap.y]#[self.arcMapMask]
        CCDLambda_m = CCDMapOffset[:,wt.c.CCDMap.wavelength]#[self.arcMapMask]
            
        CCDX_c = calibrationMap[:,wt.c.calibrationMap.x]
        CCDY_c = calibrationMap[:,wt.c.calibrationMap.y]
        CCDLambda_c = calibrationMap[:,wt.c.calibrationMap.wavelength]
        xSig_c = calibrationMap[:,wt.c.calibrationMap.sigX]  
        ySig_c = calibrationMap[:,wt.c.calibrationMap.sigY]
        
        
        #Loop through the input wavelengths, and find the closest output.
        #The best model fits in x and y (out of 2 options) is called x_best and y_best
        x_best = CCDX_c.copy()
        y_best = CCDY_c.copy()
        for k in range(0,len(CCDX_c)):
            ix, = np.where(CCDLambda_c[k] == CCDLambda_m)
            if (ix.size == 0):
                x_best[k]=0
                y_best[k]=0
            else:
                best = ix[np.argmin(np.abs(CCDY_m[ix] - CCDY_c[k]))]
                #The following lines have a potential bug because they assume there is only one "best"
                x_best[k] = CCDX_m[best]
                y_best[k] = CCDY_m[best]
             
        #eliminate devide by 0     
        xSig_c[xSig_c==0]=1
        ySig_c[ySig_c==0]=1
        
        
        diff_model_calib = np.hstack([(x_best - calibrationMap[:,wt.c.calibrationMap.x])/calibrationMap[:,wt.c.calibrationMap.sigX],
                                      (y_best - calibrationMap[:,wt.c.calibrationMap.y])/calibrationMap[:,wt.c.calibrationMap.sigY] ]) 
        print 'RMS =', np.sqrt(np.mean(diff_model_calib**2))
        
        return diff_model_calib

    
    def fitting_stats(self, fit):
    
    
        print 'Fitting Stats'
        print ''
        
        print 'Solution Array'
        print fit[0]
        print ''
        
        dict = fit[2]
        
        
        for item in dict:
            print item 
            if item=='nfev':
                desc = 'The number of function calls'
            elif item=='fvec':
                desc = 'The function evaluated at the output'
            elif item=='fjac':
               desc = 'A permutation of the R matrix of a QR factorization of the final approximate Jacobian matrix, stored column wise. Together with ipvt, the covariance of the estimate can be approximated.'
            elif item=='ipvt':
                desc = 'An integer array of length N which defines a permutation matrix, p, such that fjac*p = q*r, where r is upper triangular with diagonal elements of nonincreasing magnitude. Column j of p is column ipvt(j) of the identity matrix.'
            elif item=='qtf':
                desc = 'The vector (transpose(q) * fvec)'
    
            print desc 
            print dict[item]
            print


    def extract_order(self, CCDMap, nOrder, image = '', orderFileName = '', booShowImage = False):
        
        im = []
        #Grab image data
        if type(image)==np.ndarray:
            im = image
            imWidth = image.shape[1]
            imHeight = image.shape[0]
        elif ((type(image)==str) and (image!='')):
            hdulist = pf.open(image)
            im = pf.getdata(image)     
            imWidth = hdulist[0].header['NAXIS1']
            imHeight = hdulist[0].header['NAXIS2']     
        else:
            imWidth = self.camera.width
            imHeight = self.camera.height     
            
        CCDMapTemp = CCDMap.copy()
        orderMap = CCDMapTemp[:,wt.c.CCDMap.order]==nOrder
        yPlot = CCDMapTemp[:,wt.c.CCDMap.y][orderMap]
        xPlot = CCDMapTemp[:,wt.c.CCDMap.x][orderMap]
        LambdaPlot = CCDMapTemp[:,wt.c.CCDMap.wavelength][orderMap]
        
        sortIndex = np.argsort(yPlot)
        yPlot = yPlot[sortIndex]
        xPlot = xPlot[sortIndex]
        LambdaPlot = LambdaPlot[sortIndex]
        
#         yRange = [min(yPlot),max(yPlot)]    
#         if yRange[0]<0: yRange[0]=0
#         if yRange[1]>imHeight: yRange[1]=imHeight

        yRangeFilter = ((yPlot>=0) & (yPlot<imHeight))
        yPlot = yPlot[yRangeFilter]
        xPlot = xPlot[yRangeFilter]
        LambdaPlot = LambdaPlot[yRangeFilter]

        fLambda = interpolate.interp1d(yPlot, LambdaPlot)
        fX = interpolate.interp1d(yPlot, xPlot, 'quadratic', bounds_error=False)

        yRange = [min(yPlot),max(yPlot)]
        newX, newY, newLambdas = wt.calculate_from_Y(yRange, fX, fLambda)
        
        xRangeFilter = ((newX>=0) & (newX<imWidth))
        newX = newX[xRangeFilter]
        newY = newY[xRangeFilter]
        newLambdas = newLambdas[xRangeFilter]
        
        if orderFileName!='':
            a = np.zeros((len(newX),3))
            a[:,0] = newX
            a[:,1] = newY
            a[:,2] = newLambdas
            np.savetxt(orderFileName, a)
        
#         if type(im)==np.ndarray: 
        flux = wt.ia.extract_order(newX, newY, image, booShowImage)
#         else:
#             flux = np.zeros(len(newLambdas))
            
        return newLambdas, flux
    
    
    def read_multiple_laser_csv(self, folder, darkScale = 1):
        
        dark = self.create_dark(folder + 'darks/', scale = darkScale)
        
        # read all files into allFileNames
        _,_,allFileNames = os.walk(folder).next()
        allFileNames = np.array(allFileNames)
        
        #retrieves the unique line filenames
        linesList = []
        wavelengths = []
        for i in allFileNames:
            linesList = np.append(linesList, i[:35])
            wavelengths = np.append(wavelengths, (int(i[10:14])*1e2+int(i[15:17]))/1e5)
        linesList = np.unique(linesList)
        wavelengths = np.unique(wavelengths)
        #creates an height, width, # lines array
        im_temp = np.loadtxt(folder + allFileNames[0], delimiter=';')
        im = np.zeros((im_temp.shape[1],im_temp.shape[0],len(wavelengths))) #rotated 90 CCW
        
        #loads the im array
        for i in range(len(linesList)):
            im_temp = np.zeros((im_temp.shape[0],im_temp.shape[1],32))
            for j in range(32):
                currFile = str(linesList[i]) + '_' + str(j) + '.csv'
                im_temp[:,:,j] = np.loadtxt(folder + str(currFile), delimiter=';')
            im_temp = np.sum(im_temp, axis=2)
            im[:,:,i] = np.rot90(im_temp) - dark

        return im, wavelengths
    
    
    def read_multiple_science_csv(self, folder, darkScale = 1):
        
        dark = self.create_dark(folder + 'darks/', scale = darkScale)
        
        # read all files into allFileNames
        _,_,allFileNames = os.walk(folder).next()
        allFileNames = np.array(allFileNames)
        
        #loads the dark array
        im = np.loadtxt(folder + str(allFileNames[0]), delimiter=';')
        
        #loads the im array
        for currFile in allFileNames[1:]:
            im += np.loadtxt(folder + str(currFile), delimiter=';')
        im = np.rot90(im)

        return im, dark


    def find_peaks_im(self, im, wavelengths):
        a = np.where(im[:,:,0]==np.max(im[:,:,0]))
        array_out = np.array([[a[1][0], a[0][0],wavelengths[0], 10,0]])
        
        #hack to fix bad peak in lambda = 1.55
        im[603,287,6] = 0
        im[603,279,6] = 0
        
        
        #Extracts the wavelength list and peaks (max)
        for i in range(1,len(wavelengths)):
            a = np.where(im[:,:,i]==np.max(im[:,:,i]))
            array_out = np.append(array_out,np.array([[a[1][0], a[0][0], wavelengths[i], 10,0]]), axis=0) #lambda, x, y  (lambda, col, row)

        return array_out
        
        
    def create_dark(self, folder, scale = 1):
        # read all files into allFileNames
        _,_,allFileNames = os.walk(folder).next()
        allFileNames = np.array(allFileNames)
        
        #loads the dark array
        dark = np.loadtxt(folder + str(allFileNames[0]), delimiter=';')
        for currFile in allFileNames[1:]:
            dark += np.loadtxt(folder + str(currFile), delimiter=';')
        dark = np.rot90(dark)
        
        dark *= scale
        print 'Dark Created'
        return dark
   
   
    def create_flat(self, folder, booNormalise = True):
        # read all files into allFileNames
        _,_,allFileNames = os.walk(folder).next()
        allFileNames = np.array(allFileNames)
        
        #loads the flat array
        flat_temp = np.loadtxt(folder + str(allFileNames[0]), delimiter=';')
        flat = np.zeros((flat_temp.shape[0],flat_temp.shape[1],len(allFileNames)))
        for i in range(len(allFileNames)):
            flat[:,:,i] = np.loadtxt(folder + str(allFileNames[i]), delimiter=';')
        flat = np.median(flat,2)
        flat = np.rot90(flat)
        
        if booNormalise==True: flat /= np.max(flat)
        print 'Flat Created'
        return flat
        
        
    def stitch(self, cam1, cam2, cam3, offset):
       
        h=self.camera.height
        w=self.camera.width
        
#         fullImage = np.zeros((h*3,w))
        
#         posx_rel_cam2 = posx - posx[1]
#         posx_px_rel_cam2 = posx_rel_cam2*1000./self.camera.pSize 
#         posy_rel_cam2 = posy - posy[1]
#         posy_px_rel_cam2 = posy_rel_cam2*1000./self.camera.pSize 
        
        fullImage = np.zeros((h, w + offset[1,0]))
        
        fullImage[0:640,0:512]=cam2
        fullImage[offset[0,1]:offset[0,1]+640,offset[0,0]:offset[0,0]+512]=cam1
        fullImage[offset[1,1]:offset[1,1]+640,offset[1,0]:offset[1,1]+512]=cam3

        return fullImage
    
    
    def find_offset(self, cam1, cam2, cam3, lambda1, lambda2, lambda3):
        
        peaks1 = self.find_peaks_im(cam1, lambda1)
        peaks2 = self.find_peaks_im(cam2, lambda2)
        peaks3 = self.find_peaks_im(cam3, lambda3)
        
        #get common wavelengths
        commonLambda12 = lambda1[np.in1d(lambda1,lambda2)]
        commonLambda13 = lambda1[np.in1d(lambda1,lambda3)]
        
        #read idx of slices to be removed from arc files
        #repeated points are removed from the lower(<y) camera because it gets written over by the camera above.
        removeCam2 = np.where(np.in1d(lambda2,lambda1)==True)[0]
        removeCam1 = np.where(np.in1d(lambda1,lambda3)==True)[0]
        
        #get coords of common points in cam1
        x1 = np.zeros(len(commonLambda12))
        y1 = np.zeros(len(commonLambda12))
        for i in range(len(commonLambda12)):
            x1[i],y1[i] =  peaks1[np.where(peaks1[:,2]==commonLambda12[i])[0],0:2][0]
        
        #get coords of common points in cam2
        x2 = np.zeros(len(commonLambda12))
        y2 = np.zeros(len(commonLambda12))
        for i in range(len(commonLambda12)):
            x2[i],y2[i] =  peaks2[np.where(peaks2[:,2]==commonLambda12[i])[0],0:2][0]
        
        #get first offset (to be applied to cam1)
        offset = np.zeros((2,2))
        offset[0,:] = np.average(x2-x1), np.average(y2-y1)
        
        print x2, y2
        print x1, y1
        print x2-x1, y2-y1
        print ''
        
        #get coords of common points in cam1
        x1 = np.zeros(len(commonLambda13))
        y1 = np.zeros(len(commonLambda13))
        for i in range(len(commonLambda13)):
            x1[i],y1[i] =  peaks1[np.where(peaks1[:,2]==commonLambda13[i])[0],0:2][0]
        
        #get coords of common points in cam3    
        x3 = np.zeros(len(commonLambda13))
        y3 = np.zeros(len(commonLambda13))
        for i in range(len(commonLambda13)):
            x3[i],y3[i] =  peaks3[np.where(peaks3[:,2]==commonLambda13[i])[0],0:2][0]
        
        # complete offset points and add displace on cam1 (to be used in camera3)
        offset[1,:] = np.average(x1-x3), np.average(y1-y3)
        offset[1,:] +=  offset[0,:]
        offset = offset.astype(int)

        print x1, y1
        print x3, y3
        print x1-x3, y1-y3        
        
        return offset, removeCam2, removeCam1
    

    def write_csv_from_im(self, fileName, im):
        
        np.savetxt(fileName, im, fmt='%05d', delimiter=';')


class cameras():    
    name = 'blue'
    raw = ''
    name = ''
    width = ''
    height = ''
    offsetX = 0
    offsetY = 0
    CCDMap = []
    CCDMapOffset = []
    
class beams():
    pass

class optics():
    pass

    
class plot():
       
    def ccd_map(self, CCDMap, canvasSize=1, backImage='', title = '', backgroundFile = ''):
        #Plot CCD map 
        
        colorTable = np.array((wt.optics.wav2RGB(CCDMap[:,wt.c.CCDMap.wavelength], CCDMap[:,wt.c.CCDMap.intensity]))) 
        
        if backImage!='':
            im = backImage
            plt.imshow(im, origin='lower')
            
        if backgroundFile!='':
            hdulist = pf.open(backgroundFile)
            imWidth = hdulist[0].header['NAXIS1']
            imHeight = hdulist[0].header['NAXIS2']
            im = hdulist[0].data
            im[im<0] = 0
            im[np.isnan(im)] = 0
#             im /= np.max(im)
#             im = np.max(im)
            im = np.log10(im) #Remove this line for Hg
            plt.imshow(im, origin='lower')
            plt.set_cmap(plt.cm.Greys_r)
            plt.axis([ 0, imWidth * canvasSize , 0, imHeight * canvasSize])
            
        plt.scatter(CCDMap[:,wt.c.CCDMap.x], CCDMap[:,wt.c.CCDMap.y] ,s=20, color=colorTable , marker='o', alpha =.5)
        if title=='': plt.title('CCD Map')
        plt.ylabel('pixels')
        plt.xlabel('pixels')
        plt.show()

    
    def calib_points(self, calibrationDataFileName, backImage='', booLabels = False, canvasSize=2, title='', booInteractive = False,  backgroundFile = ''):

        #Loads from calibration output file
        calibrationMap = wt.ia.read_full_calibration_data(calibrationDataFileName)
        if calibrationMap.shape[0]==0: 
            print  'Reading calibration file returned an empty array.'
            return
        
        if backImage!='':
            im = backImage
            plt.imshow(im, origin='lower')

        if backgroundFile!='':
            hdulist = pf.open(backImage)
            imWidth = hdulist[0].header['NAXIS1']
            imHeight = hdulist[0].header['NAXIS2']
            im = hdulist[0].data
            fig = plt.figure()
            plt.imshow(im, origin='lower')#,extent=[0, imWidth , 0, imHeight])
            plt.set_cmap(plt.cm.Greys_r)
        
        
        
        plt.scatter(calibrationMap[:,wt.c.calibrationMap.x], calibrationMap[:,wt.c.calibrationMap.y] ,s=50, color="red" , marker='o', alpha = 0.5, label='Calibration Data')
#         plt.axis([-imWidth/2 * canvasSize, imWidth/2 * canvasSize, -imHeight/2 * canvasSize, imHeight/2 * canvasSize])
#         plt.axis([0, imWidth * canvasSize, 0, imHeight * canvasSize])
        plt.legend()
        if title=='': title = str(calibrationMap.shape[0])+' point(s) found in the calibration file'
        
        if booLabels:
            fullLabel = calibrationMap[:,wt.c.calibrationMap.wavelength] 
            for x, y, label in zip(calibrationMap[:,wt.c.calibrationMap.x], calibrationMap[:,wt.c.calibrationMap.y], fullLabel ):
                plt.annotate(
                    label, 
                    xy = (x, y), xytext = (0,-25),
                    textcoords = 'offset points', ha = 'right', va = 'bottom',
                    bbox = dict(boxstyle = 'round,pad=0.5', fc = 'white', alpha = 1))
#                    arrowprops = dict(arrowstyle="wedge,tail_width=1.",
#                                fc=(1, 0, 0), ec=(1., 0, 0),
#                                patchA=None,
#                                relpos=(0.2, 0.8),
#                                connectionstyle="arc3,rad=-0.1"), size=10)

        plt.title(title) 
        if booInteractive>0:
            plt.ion()
            plt.draw()
            a = plt.ginput(-1)
            plt.close()
            plt.ioff()
        else:
            plt.show()

            a = []
            
        return a
    
    def calib_points_CCDMap(self, calibrationMap, backImage = '', backgroundFile='', CCDMap = [], booLabels = False, canvasSize=2, title='', booInteractive = 0, arcMapMask = []):
                 
        if backImage!='':
            im = backImage
            plt.imshow(im, origin='lower')

        if backgroundFile!='':
            hdulist = pf.open(backgroundFile)
            imWidth = hdulist[0].header['NAXIS1']
            imHeight = hdulist[0].header['NAXIS2']
            im = hdulist[0].data
#             fig = plt.figure()
            plt.imshow(im, origin='lower')#,extent=[0, imWidth , 0, imHeight])
            plt.set_cmap(plt.cm.Greys_r)

        plt.scatter(calibrationMap[:,wt.c.calibrationMap.x], calibrationMap[:,wt.c.calibrationMap.y] ,s=50, color="red" , marker='o', alpha = 0.5, label='Calibration Data')
        if booLabels:
            fullLabel = calibrationMap[:,wt.c.calibrationMap.wavelength]
            for x, y, label in zip(calibrationMap[:,wt.c.calibrationMap.x], calibrationMap[:,wt.c.calibrationMap.y], fullLabel ):
                plt.annotate(
                    label, 
                    xy = (x, y), xytext = (0,-25),
                    textcoords = 'offset points', ha = 'right', va = 'bottom',
                    bbox = dict(boxstyle = 'round,pad=0.5', fc = 'white', alpha = 1))

        if CCDMap!=[]: 
            CCDX = CCDMap[:,wt.c.CCDMap.x][arcMapMask]
            CCDY = CCDMap[:,wt.c.CCDMap.y][arcMapMask] 
            CCDWavelength = CCDMap[:,wt.c.CCDMap.wavelength][arcMapMask] 
            CCDIntensity = CCDMap[:,wt.c.CCDMap.intensity][arcMapMask]
            CCDOrder = CCDMap[:,wt.c.CCDMap.order][arcMapMask]
            plt.scatter(CCDX, CCDY ,s=50, color="blue" , marker='o', alpha = 0.5, label='Model Data')
            if booLabels:
#                fullLabel = ['('+str(CCDX[x])+', '+str(CCDY[x])+')'+str(CCDLambda[x]) for x in np.arange(len(CCDX))]
#                 fullLabel = [str(CCDWavelength[x]) for x in np.arange(len(CCDWavelength))]
                fullLabel = CCDMap[:,wt.c.CCDMap.wavelength][arcMapMask].astype('|S10')
#                 print fullLabel
                for x, y, label in zip(CCDX, CCDY, fullLabel ):
                    plt.annotate(
                        label, 
                        xy = (x, y), xytext = (0,-25),
                        textcoords = 'offset points', ha = 'right', va = 'bottom',
                        bbox = dict(boxstyle = 'round,pad=0.5', fc = 'white', alpha = 0.9))
        
        plt.legend()
        if title=='':title = str(len(calibrationMap[:,wt.c.calibrationMap.x]))+' point(s) found in the calibration image'
        plt.title(title) 
#         plt.axis([-imWidth/2 * canvasSize, imWidth/2 * canvasSize, -imHeight/2 * canvasSize, imHeight/2 * canvasSize])
        
        if booInteractive:
            plt.draw()
            removePoints = plt.ginput(0, timeout = -1)
            plt.close()
            for i in removePoints:
                distCCDX = np.abs(CCDMap[:,wt.c.CCDMap.x] - i[0])
                distCCDY = np.abs(CCDMap[:,wt.c.CCDMap.y] - i[1])
                distCalibX = np.abs(calibrationMap[:,wt.c.calibrationMap.x] - i[0])
                distCalibY = np.abs(calibrationMap[:,wt.c.calibrationMap.y] - i[1])
                closestCCD = np.min(distCCDX+distCCDY)
                closestCalib = np.min(distCalibX+distCalibY)
                if closestCCD<closestCalib:
                    if closestCCD<booInteractive:
                        removeCCDIx = np.where((distCCDX+distCCDY)==closestCCD)[0][0]
                        arcMapMask[removeCCDIx] = 0
                        print 'removed CCDMap: ', CCDMap[removeCCDIx,wt.c.CCDMap.wavelength]
                else:
                    if closestCalib<booInteractive:
                        removeCalibIx = np.where((distCalibX+distCalibY)==closestCalib)[0][0]
                        print 'removed CailbPoint: ', calibrationMap[removeCalibIx,wt.c.calibrationMap.wavelength]
                        calibrationMap = np.delete(calibrationMap,removeCalibIx,0)
            print 'CCDMap points: ',CCDMap[arcMapMask].shape[0]
            print 'Calibration Points: ' ,calibrationMap.shape[0]
        else:
            plt.show()
            
        return calibrationMap, arcMapMask
        

#         
#         
#         
#         #Loads from calibration output file
#         imageMapX, imageMapY, imageMapWavelength , imageMapXsig , imageMapYsig  = wt.ia.read_full_calibration_data(calibrationDataFileName)
#         if imageMapX==[]: return
#         
#         #Plot
#         hdulist = pf.open(backImageFileName)
#         imWidth = hdulist[0].header['NAXIS1']
#         imHeight = hdulist[0].header['NAXIS2']
#         im = pf.getdata(backImageFileName) 
# #         imNorm = wt.ic.normalise_image(im) 
#         imNorm = im
#         fig = plt.figure()
#         ax1 = fig.add_subplot(111)
#         # (left, right, bottom, top) 
# #         plt.imshow(imNorm,extent=[-imWidth/2 , imWidth/2 , imHeight/2 , -imHeight/2])
#         plt.imshow(imNorm, origin='lower')#,extent=[0, imWidth , 0, imHeight])
#         plt.set_cmap(plt.cm.Greys_r)
#         
#         #-0.5 to match centre of the pixel fits convention
#         
#         ax1.scatter(imageMapX, imageMapY ,s=50, color="red" , marker='o', alpha = 0.5, label='Calibration Data')
# #         plt.axis([-imWidth/2 * canvasSize, imWidth/2 * canvasSize, -imHeight/2 * canvasSize, imHeight/2 * canvasSize])
# #         plt.axis([0, imWidth * canvasSize, 0, imHeight * canvasSize])
# 
#         plt.legend()
#         if title=='': title = str(len(imageMapX))+' point(s) found in the calibration file'
#         
#         if booLabels:
#             fullLabel = [str(imageMapX[x])+','+ str(imageMapY[x]) for x in np.arange(len(imageMapWavelength))] 
#             for x, y, label in zip(imageMapX, imageMapY, fullLabel ):
#                 plt.annotate(
#                     label, 
#                     xy = (x, y), xytext = (0,-25),
#                     textcoords = 'offset points', ha = 'right', va = 'bottom',
#                     bbox = dict(boxstyle = 'round,pad=0.5', fc = 'white', alpha = 1))
# #                    arrowprops = dict(arrowstyle="wedge,tail_width=1.",
# #                                fc=(1, 0, 0), ec=(1., 0, 0),
# #                                patchA=None,
# #                                relpos=(0.2, 0.8),
# #                                connectionstyle="arc3,rad=-0.1"), size=10)
# 
#         plt.title(title) 
#         plt.show()

        