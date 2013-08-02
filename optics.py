import numpy as np

def Grating(u, s, l, nOrder, Lambda, d):
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
