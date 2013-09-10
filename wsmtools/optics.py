import numpy as np

def ccd_loop(SEDMap, Beams, Optics, Camera, stheta): #, intNormalize,Interpolate=False,BackImage='',GaussFit=False):
    '''
    Computes the projection of n beams of monochromatic light passing through an optical system. 

    Parameters
    ----------
    SEDMap : np.array
        n x 2 np.array with wavelength, Energy
    
    p : np np.array (optional)
        (beam phi, beam theta, prism1 phi, prism1 theta, prism2 phi, prism2 theta, grating phi, grating theta, grating alpha,
         blaze period (microns), focal length(mm), distortion term) <-- optical arrangement
         
    args : np np.array
        (SEDMode(0=Max, 1=Random, 2=Sun, 3=from specFile, 4=from CalibFile), Plot?, specFile, Normalize intensity? (0=no, #=range), Distort?, Interpolate, PlotCalibPoints) <-- other options
   
    Returns
    -------
    x : np np.array
        x coordinate of the target point 
    y : np np.array
        x coordinate of the target point 
    lambda : np np.array
        wavelength at x,y

    Notes
    -----  
    '''   
       
    dataOut=np.zeros(5)
    CCDX = CCDY = CCDLambda = CCDIntensity = CCDOrder = np.array([])
    
    pixelSize = float(Camera[CamerasPSize])
    fLength = float(Camera[CamerasFLength])
    blaze_angle = stheta #Approximately np.arctan(2)

    #Retrieves max and min lambdas for intensity calculation
    minLambda=min(SEDMap[SEDMapLambda])
    maxLambda=max(SEDMap[SEDMapLambda])
#    allFlux=np.array([0])
#    allLambdas=np.array([0])
    
    #One per beam
    for Beam in Beams:
        
        #Navigates orders within the range given   
        for nOrder in range(minOrder, maxOrder, deltaOrder):
    
            #loop lambda for current order
            for i in np.arange(len(SEDMap[SEDMapLambda])): 
                
                Lambda=SEDMap[SEDMapLambda][i]
                inI=SEDMap[SEDMapIntensity][i]
                            
                #rhea hack GPeriod from optics, todo
                GPeriod=31.50321471
                
                #the wavelength range is from the blaze wavelength of the next order and the blaze wavelength of the previous order
                if (Lambda >= abs(2*GPeriod*np.sin(blaze_angle)/(nOrder+1)) and Lambda <= abs(2*GPeriod*np.sin(blaze_angle)/(nOrder-1))):
    
                    #Computes the unit vector that results from the optical system for a given wavelength and order
                    #This is the actual tracing of the ray for each wavelength             
                    start_time = time.time()
                    v, isValid = ray_trace_flex(Beam, Lambda, nOrder, Optics, blaze_angle)
                    elapsed_time = time.time() - start_time
                    #print 'Elepased time: ' + str(elapsed_time)
                    
                    if isValid: #no errors in calculation, within 1 order of blaze wavelength and beam makes it back through the prism
                        x=v[0]*fLength*1000/pixelSize # x-coord in focal plane in pixels
                        z=v[2]*fLength*1000/pixelSize # z-coord in focal plane in pixels
            
                        outI=Intensity(Lambda, minLambda, maxLambda)    
                        
                        #Add results to output arrays
                        if len(CCDX)==0:
                            CCDX = np.array([x])
                            CCDY = np.array([z])
                            CCDLambda = np.array([Lambda])
                            CCDIntensity = np.array([inI*outI])
                            CCDOrder = np.array([nOrder])
                        else:
                            CCDX = np.append(CCDX,[x],0)
                            CCDY = np.append(CCDY,[z],0)
                            CCDLambda = np.append(CCDLambda,[Lambda],0)
                            CCDIntensity = np.append(CCDIntensity,[inI*outI],0)
                            CCDOrder = np.append(CCDOrder,[nOrder],0)
        
    return CCDX, CCDY, CCDLambda, CCDIntensity, CCDOrder

def ray_trace_flex(Beam, Lambda, nOrder, Optics, blaze_angle):
    ''' Traces a beam through the spectrograph. 
    Spectrograph frame of reference, from the opposite end of the camera looking at the camera
    x=to the right, y=to camera, z=up
    u*=beam, n*=surface normals
    s=grating, perp to the grooves.
    l=grating, parallel to the grooves.
    d=blaze period'''
    
    v = Beam
    
    #loops through optics array
    for i in range(len(Optics)):

        if i==0: n_i=n(Lambda,'air')
        optType=Optics[i][OpticsType]
        n_r=n(Lambda,Optics[i][OpticsN])
        
        if optType==OpticsBoundary:
            surfNormal=Optics[i][OpticsCoords1]           
            v_out = Snell3D(n_i, n_r, v, surfNormal)
        elif optType==OpticsRGrating:
            isValid=False
            GPeriod=Optics[i][OpticsGPeriod]
            s=Optics[i][OpticsCoords1]
            l=Optics[i][OpticsCoords2]
            v_out, isValid = RGrating(v, s, l, nOrder, Lambda, GPeriod)
        elif optType==OpticsVPHGrating:
            isValid=False
            GPeriod=Optics[i][OpticsGPeriod]
            s=Optics[i][OpticsCoords1]
            l=Optics[i][OpticsCoords2]
            v_out, isValid = VPHGrating(v, s, l, nOrder, Lambda, GPeriod)
            isValid = True
#            else:
#                return v, isValid
            
        n_i = n_r
        v = v_out

#                           
#        """Vector transform due to first surface"""
#        u = Snell3D(nAir, nPrism, u, n1)
#                         
#        """Vector transform due to second surface"""
#        u = Snell3D(nPrism, nAir, u, n2)
#        
#        """grating dispersion"""              
#        u, isValid = RGrating(u, s, l, nOrder, Lambda, d)
#        
#        if isValid:
#            """Vector transform due to third surface"""
#            u = Snell3D(nAir, nPrism, u, n4)
#                       
#            """Vector transform due to fourth surface"""
#            u = Snell3D(nPrism, nAir, u, n5)
            
    return v, isValid

def RGrating(u, s, l, nOrder, Lambda, d):
    """Computes the new direction of a vector when hitting a grating."""
    
    isValid=False

    n = np.cross(s, l) 
    u_l = np.dot(u, l)
    u_s = np.dot(u, s) + nOrder*Lambda/d  
    if (1-u_l**2 -u_s**2)>=0: 
        u_n = np.sqrt(1- u_l**2 - u_s**2)
        u = u_l*l + u_s*s + u_n*n
        isValid=True

    return u, isValid     

def VPHGrating(u, s, l, nOrder, Lambda, d):
    """Computes the new direction of a vector when hitting a grating."""
    
    isValid=False

    n = np.cross(s, l) 
    u_l = np.dot(u, l)
    u_s = np.dot(u, s) + nOrder*Lambda/d  
    if (1-u_l**2 -u_s**2)>=0: 
        u_n = np.sqrt(1- u_l**2 - u_s**2)
        u = u_l*l + u_s*s + u_n*n
        isValid=True

    return u, isValid     

def Snell3D(n_i,n_r,u,n):
    """Computes the new direction of a vector when changing medium. 
    n_i, n_r = incident and refractive indices"""
 
    u_p = u - np.dot(u,n)*n
    u_p /= np.linalg.norm(u_p)
    
    theta_i = np.arccos(np.dot(u,n))
    
    if n_i*np.sin(theta_i)/n_r<=1:
        theta_f = np.arcsin(n_i*np.sin(theta_i)/n_r)
        u = u_p*np.sin(np.pi-theta_f) + n*np.cos(np.pi-theta_f)    

    return u       

def n(Lambda, medium='air', t=18, p=101325):
    
    if medium=='air':
        '''Function that calculates the refractive index of air for a given wavelength'''
        n = (0.0472326 * (173.3 - (1/Lambda)**2)**(-1))+1
    elif medium=='nkzfs8':  #nkzfs8 only for the moment. Should change into a sellmeier eq with material as input parameter
        '''Function that calculates the refractive index of the prism for a given wavelength'''
        x = float(Lambda)
        n = np.sqrt(1 + (1.62693651*x**2)/(x**2-0.010880863) + 0.24369876*x**2/(x**2-0.0494207753) + 1.62007141*x**2/(x**2-131.009163) )

    return n

def distort( inX, inY, K=0, Xc=0 , Yc=0):        
    '''
    Transforms coordinate values based on Brown's 1966 model.
       
    Parameters
    ----------
    inX : np np.array
        1 x n np.array of X-coordinates to be transformed.
    inY : np np.array
        1 x n np.array of Y-coordinates to be transformed.
    K : np np.array
        1 x n np.array of distortion factors.
        The length of the np.array determines the order of the polynomial. 
    Xc : float
        X-coord of distortion center.
    Yc : float
        Y-coord of distortion center.

    Returns
    -------
    outX : np np.array
        1 x n np.array of X-coordinates transformed.
    outY : np np.array
        1 x n np.array of Y-coordinates transformed.
      
    Notes
    -----
    
    '''
    
    Delta=1
    r = np.sqrt((inX-np.ones(len(inX))*Xc)**2+(inY-np.ones(len(inY))*Yc)**2)
    r /= max(r)
       
#    for i in np.arange(len(K)):
#        Delta += K[i] * r**((i+1)*2)
    Delta += K * r**2
        
    outX = (inX-np.ones(len(inX))*Xc)*Delta
    outY = (inY-np.ones(len(inY))*Yc)*Delta

    return outX, outY

def wav2RGB(Lambda, Intensity):
    
    """Converts Lambda into RGB"""
    out=np.array([0,0,0])
    Intensity /= np.max(np.abs(Intensity),axis=0)   
    for i in range(Lambda.size):
    
        w = Lambda[i]
        I = Intensity[i]
    
        # colour
        if w >= .380 and w < .440:
            R = -(w - .440) / (.440 - .350)
            G = 0.0
            B = 1.0
        elif w >= .440 and w < .490:
            R = 0.0
            G = (w - .440) / (.490 - .440)
            B = 1.0
        elif w >= .490 and w < .510:
            R = 0.0
            G = 1.0
            B = -(w - .510) / (.510 - .490)
        elif w >= .510 and w < .580:
            R = (w - .510) / (.580 - .510)
            G = 1.0
            B = 0.0
        elif w >= .580 and w < .645:
            R = 1.0
            G = -(w - .645) / (.645 - .580)
            B = 0.0
        elif w >= .645 and w <= .780:
            R = 1.0
            G = 0.0
            B = 0.0
        else:
            R = 1.0
            G = 1.0
            B = 1.0
    
        # intensity correction
        if w >= .3800 and w < .4200:
            SSS = 0.3 + 0.7*(w - .3500) / (.4200 - .3500)
        elif w >= .4200 and w <= .7000:
            SSS = 1.0
        elif w > .7000 and w <= .7800:
            SSS = 0.3 + 0.7*(.7800 - w) / (.7800 - .7000)
        else:
            SSS = 1.0

        SSS *= (I)

        out=np.vstack((out,np.array( [float(SSS*R), float(SSS*G), float(SSS*B)]))) 
        
    return out[1:,]

def fftshift1D(inImage, shift):
    '''
    This program shifts an image by sub-pixel amounts.
       
    Parameters
    ----------
    inImage :  image
        Input image
    shift : array
        (x,y) pixel shift
    Returns
    -------
    outImage : Image
        Shifted Image
      
    Notes
    -----
    '''  
    
## An example - shift a rough Gaussian by 0.5 pixels.
#in_image = np.array([[0,0,0.06,0,0],[0,0.25,0.5,0.25,0],[0.06,0.5,1,0.5,0.06],[0,0.25,0.5,0.25,0],[0,0,0.06,0,0]])
#out_image = fftshift(in_image,(0,0))
#print out_image

    ftin = np.fft.fft(inImage)
    sh = inImage.shape
    
    #The following line makes a meshgrid np.array as floats. Not sure if there is a neater way.
    xy = np.mgrid[0:sh[0]] + 0.0
    xy[:] = (((xy[0,:] + sh[0]/2) % sh[0]) - sh[0]/2)/float(sh[0])
    #xy[1,:,:] = (((xy[1,:,:] + sh[1]/2) % sh[1]) - sh[1]/2)/float(sh[1])

    return np.real(np.fft.ifft(ftin*np.exp(np.complex(0,-2*np.pi)*(xy[0,:,:]*shift[0] + xy[1,:,:]*shift[1])))) 

def fftshift(inImage, shift):

    '''
    This program shifts an image by sub-pixel amounts.
       
    Parameters
    ----------
    inImage :  image
        Input image
    shift : array
        (x,y) pixel shift
    Returns
    -------
    outImage : Image
        Shifted Image
      
    Notes
    -----
    '''  
    
## An example - shift a rough Gaussian by 0.5 pixels.
#in_image = np.array([[0,0,0.06,0,0],[0,0.25,0.5,0.25,0],[0.06,0.5,1,0.5,0.06],[0,0.25,0.5,0.25,0],[0,0,0.06,0,0]])
#out_image = fftshift(in_image,(0,0))
#print out_image

    ftin = np.fft.fft2(inImage)
    sh = inImage.shape
    
    #The following line makes a meshgrid np.array as floats. Not sure if there is a neater way.
    xy = np.mgrid[0:sh[0],0:sh[1]] + 0.0
    xy[0,:,:] = (((xy[0,:,:] + sh[0]/2) % sh[0]) - sh[0]/2)/float(sh[0])
    xy[1,:,:] = (((xy[1,:,:] + sh[1]/2) % sh[1]) - sh[1]/2)/float(sh[1])
    
    return np.real(np.fft.ifft2(ftin*np.exp(np.complex(0,-2*np.pi)*(xy[0,:,:]*shift[0] + xy[1,:,:]*shift[1]))))

def Intensity(Lambda, minLambda, maxLambda):
    '''
    Retrieves or calculates the expected relative intensity based on distance from the central lambda value
       
    Parameters
    ----------
    Lambda : np.float32
        Wavelength.
    minLambda : np.float32
        Lower end of wavelength range.
    maxLambda : np.float32
        Higher end of wavelength range. 
   
    Returns
    -------
    z : np.float32
        result 0 to 1.
      
    Notes
    -----
    
    '''
#todo fix
#    x = 0.5*(float(Lambda) - 0.5*(float(maxLambda) - float(minLambda)))/(float(maxLambda) - float(minLambda))
#    x = (((float(Lambda) - float(minLambda))/(float(maxLambda) - float(minLambda)))-0.5)*2
#    if x!=0:
#        z=np.sin(x*np.pi)/(x*np.pi)
    z=1

#    print x,z
                
#        """Look for wavelength"""
#        EMax = max(spectrum_dict.iteritems(), key=operator.itemgetter(1))[1]   
#        a= spectrum_dict[str(Lambda)]            
#        if (1==2):
#            E=random.random() * (EMax-EMin)+EMin            
#        else:
#            j = spectrum_dict
#            topLambda = max([i for i in j if i >= 5])
#            topLambda = max(spectrum_dict.iteritems(), key=operator.itemgetter(1))[1]
#            bottomLambda=max(spectrum_dict.iteritems(), key=operator.itemgetter(1))[1]
#            topE=max(spectrum_dict.iteritems(), key=operator.itemgetter(1))[1]
#            bottomE=max(spectrum_dict.iteritems(), key=operator.itemgetter(1))[1]
#            E=((topE-bottomE)((Lambda-bottomLambda) /(topLambda-bottomLambda)))+bottomE            
#            Intensity=(E-EMin)/(EMax-EMin)*255
       
    return z