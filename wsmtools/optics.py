import numpy as np

import constants as c


def ccd_loop(SEDMap, beams, optics, camera, stheta, orders): 
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
       
    CCDMap = np.array([[]])
    CCDMap3D = np.array([[[]]])
    steps3D = [0,30,1,2,2,1,30]
    
    pixelSize = camera.pSize
    fLength = camera.fLength
    blazeAngle = stheta #Approximately np.arctan(2) #todo check how to make this general
#     blazeAngle = np.arctan(2)  #Approximat                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          ly np.arctan(2) #todo check how to make this general

    #Retrieves max and min lambdas for intensity calculation
    minLambda=min(SEDMap[:,c.sedMap.wavelength])
    maxLambda=max(SEDMap[:,c.sedMap.wavelength])

    
    for i in range(len(optics[:,c.optics.coords1x])):
        if (optics[i,c.optics.type]==c.opticTypes.eGrating or optics[i,c.optics.type]==c.opticTypes.VPHGrating):
            GPeriod=optics[i,c.optics.gPeriod]

    beamID = 0   
    #One per beam
    for beam in beams:
        
        #Navigates orders within the range given   
        for nOrder in orders:
#             print 'Order:',nOrder
            #loop lambda for current order
            for i in np.arange(len(SEDMap[:,c.sedMap.wavelength])): 
                
                Lambda = SEDMap[i,c.sedMap.wavelength]
                inI = SEDMap[i,c.sedMap.intensity]
#                 print 'wavelength:', Lambda            
                
                #the wavelength range is from the blaze wavelength of the next order and the blaze wavelength of the previous order
                if (abs(Lambda*(nOrder+0.5)) >= abs(2*GPeriod*np.sin(blazeAngle)) and abs(Lambda*(nOrder-0.5)) <= abs(2*GPeriod*np.sin(blazeAngle))):

                    #Computes the unit vector that results from the optical system for a given wavelength and order
                    #This is the actual tracing of the ray for each wavelength             
                    v, v3D, isValid = ray_trace_flex(beam, Lambda, nOrder, optics, blazeAngle)
                    
                    if isValid: #no errors in calculation, within 1 order of blaze wavelength and beam makes it back through the prism
#                         print Lambda, nOrder,'^^^^^^---------------------------------------------------------------------Good Vector'
                        x=v[0]*fLength*1000/pixelSize # x-coord in focal plane in pixels
                        z=v[2]*fLength*1000/pixelSize # z-coord in focal plane in pixels
                        
                        x3D = np.zeros(len(v3D[:,0]))
                        y3D = np.zeros(len(v3D[:,0]))
                        z3D = np.zeros(len(v3D[:,0]))

                        for i in range(1,len(v3D[:,0])):
                            x3D[i] = x3D[i-1]+v3D[i,0]*steps3D[i]
                            y3D[i] = y3D[i-1]+v3D[i,1]*steps3D[i]
                            z3D[i] = z3D[i-1]+v3D[i,2]*steps3D[i]
                            
                            
#                         x3D = np.array(v3D[:,0]*steps3D)
#                         y3D = np.array(v3D[:,1]*steps3D)
#                         z3D = np.array(v3D[:,2]*steps3D)
                                  
                        outI = Intensity(Lambda, minLambda, maxLambda)    
                    
                        #Add results to output arrays
                        if CCDMap.shape[1]==0:
                            CCDMap = np.append(CCDMap,[[x,z,Lambda,inI*outI,nOrder,beamID]],axis=1)
                        else:
                            CCDMap = np.append(CCDMap,[[x,z,Lambda,inI*outI,nOrder,beamID]],axis=0)

                        if CCDMap3D.shape[2]==0:
                            CCDMap3D = np.array([[x3D,y3D,z3D,np.ones(len(x3D))*Lambda,np.ones(len(x3D))*inI*outI,np.ones(len(x3D))*nOrder,np.ones(len(x3D))*beamID]])
                        else:
                            CCDMap3D = np.append(CCDMap3D,np.array([[x3D,y3D,z3D,np.ones(len(x3D))*Lambda,np.ones(len(x3D))*inI*outI,np.ones(len(x3D))*nOrder,np.ones(len(x3D))*beamID]]),axis=0)
                        
                        beamID += 1


#     print 'created ccdmap:', CCDMap.shape
    return CCDMap, CCDMap3D

def ray_trace_flex(beam, Lambda, nOrder, optics, blazeAngle):
    ''' Traces a beam through the spectrograph. 
    Spectrograph frame of reference, from the opposite end of the camera looking at the camera
    x=to the right, y=to camera, z=up
    u*=beam, n*=surface normals
    s=grating, perp to the grooves.
    l=grating, parallel to the grooves.
    d=blaze period'''
    
    v = np.array(beam[0])
    v3D = np.array([[0,0,0],beam[0]])
#     print 'input beam',v
    isValid = True

    #loops through optics array
    for i in range(len(optics[:,c.optics.coords1x])):
        #First time, incident medium is air.
        if i==0: n_i=n(Lambda,c.media.air)
        
        if isValid:
            #Optical element and refr. index for this loop.
            optType=optics[i,c.optics.type]
            n_r=n(Lambda,optics[i,c.optics.n])
            
            if optType==c.opticTypes.boundary :
                isValid=False
                surfNormal=optics[i,c.optics.coords1x:c.optics.coords1z+1] 
                v_out, isValid = Snell3D(n_i, n_r, v, surfNormal)
            elif optType==c.opticTypes.eGrating:
                isValid=False
                GPeriod=optics[i,c.optics.gPeriod]
                s=optics[i,c.optics.coords1x:c.optics.coords1z+1]
                l=optics[i,c.optics.coords2x:c.optics.coords2z+1]
                v_out, isValid = RGrating(v, s, l, nOrder, Lambda, GPeriod)
            elif optType==c.opticTypes.VPHGrating:
                isValid=False
                GPeriod=optics[i,c.optics.gPeriod]
                s=optics[i,c.optics.coords1x:c.optics.coords1z+1]
                l=optics[i,c.optics.coords2x:c.optics.coords2z+1]
                v_out, isValid = VPHGrating(v, s, l, nOrder, Lambda, GPeriod)
    
            n_i = n_r
            v = v_out
            v3D = np.append(v3D,[v],axis=0)
#             if isValid: print i,v

    return v, v3D, isValid

def RGrating(u, s, l, nOrder, Lambda, d):
    """Computes the new direction of a vector when hitting a grating."""
    
    isValid=False

    n = np.cross(s, l) #normal vector from s and l
    
    u_l = np.dot(u, l)
    u_s = np.dot(u, s) + nOrder*Lambda/d  
    
    if (1-u_l**2 -u_s**2)>=0: 
        u_n = np.sqrt(1- u_l**2 - u_s**2)
        u = u_l*l + u_s*s + u_n*n
        isValid=True    
#         print 'l,s,n,Order:',u_l,u_s,u_n, nOrder
    else:
#         print 'RGrating: Not unit vector'
        pass
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
    
    theta_i = np.arccos(np.dot(u,n))
    if np.degrees(theta_i)>90:
        n = -n
        theta_i = np.arccos(np.dot(u,n))
        
    #u_p unit vector across boundary surface. 
    #1 of the vectors that forms the basis for the result (the other is n, the normal)
    #runs in the direction of the projection of the output beam on the surface
    u_p = u - np.dot(u,n)*n
    u_p /= np.linalg.norm(u_p)
       
    isValid = False
    sin_theta_f = n_i*np.sin(theta_i)/n_r
    if sin_theta_f<=1:
        theta_f = np.arcsin(sin_theta_f)
        u = u_p*np.sin(theta_f) + n*np.cos(theta_f) #builds output vector from output angle theta_f
        isValid = True
    else:
#         print 'Snell: Total internal reflection. Incident angle(deg):', np.degrees(theta_i)
        pass
    
    return u, isValid       

def n(Lambda, medium=c.media.air, t=18, p=101325):
    
    if medium==c.media.air:
        '''Function that calculates the refractive index of air for a given wavelength'''
        n = (0.0472326 * (173.3 - (1/Lambda)**2)**(-1))+1
    elif medium==c.media.nkzfs8:  #nkzfs8 only for the moment. Should change into a sellmeier eq with material as input parameter
        '''Function that calculates the refractive index of the prism for a given wavelength'''
        x = float(Lambda)
        n = np.sqrt(1 + (1.62693651*x**2)/(x**2-0.010880863) + 0.24369876*x**2/(x**2-0.0494207753) + 1.62007141*x**2/(x**2-131.009163) )
    elif medium==c.media.NSF11:  #N-SF11 only for the moment. Should change into a sellmeier eq with material as input parameter
        '''Function that calculates the refractive index of the prism for a given wavelength'''
        x = float(Lambda)
        n = np.sqrt(1 + (1.73759695*x**2)/(x**2-0.013188707) + 0.313747346*x**2/(x**2-0.0623068142) + 1.89878101*x**2/(x**2-155.23629) )

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
    if len(inX)>0:
        
        Delta=1
        r = np.sqrt((inX-np.ones(len(inX))*Xc)**2+(inY-np.ones(len(inY))*Yc)**2)
        r /= max(r)
           
        for i in np.arange(len(K)):
            Delta += K[i] * r**((i+1)*2)
        
        outX = (inX-np.ones(len(inX))*Xc)*Delta
        outY = (inY-np.ones(len(inY))*Yc)*Delta
    
    else:
        outX, outY = inX, inY
        
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
            R = 0.0
            G = 0.0
            B = 0.0
    
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