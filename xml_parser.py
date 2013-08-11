from xml.dom import minidom
import numpy as np
from constants import *
import os, time, sys
import wsmtools as wt

def read_p(specXMLFileName):
    
    p=np.zeros(14)
    
    
    xmldoc = minidom.parse(SPEC_PATH + specXMLFileName)
    
    
#    optElements = xmldoc.getElementsByTagName('optical_element') 
    #print itemlist[0].attributes['name'].value
    
    spectrograph=xmldoc.childNodes[0]
    for specElement in spectrograph.childNodes:
        if specElement.nodeType==1:
            if specElement.hasAttribute('param'): 
                p[int(specElement.attributes.getNamedItem('param').value)] = float(specElement.firstChild.data)
            for Element in specElement.childNodes:   
                        if Element.nodeType==1:    
                            if specElement.hasAttribute('param'): 
                                p[int(specElement.attributes.getNamedItem('param').value)] = float(specElement.firstChild.data)
                            for child in Element.childNodes:
                                if ((child.nodeType==1) and child.hasAttribute('param')):
                                    p[int(child.attributes.getNamedItem('param').value)] = float(child.firstChild.data)

    return p

def write_p(p, specXMLFileName):
    
    bkp_time = time.time()
    
    newSpecXMLFileName = wt.find_specXMLFileName(specXMLFileName)
    os_command = 'cp ' + SPEC_PATH + specXMLFileName + ' ' +  SPEC_PATH + newSpecXMLFileName
    os.system(os_command)
    
    specXMLFileName = newSpecXMLFileName
    
    xmldoc = minidom.parse(SPEC_PATH + specXMLFileName)
    
    spectrograph=xmldoc.childNodes[0]
    for specElement in spectrograph.childNodes:
        if specElement.nodeType==1:
            if specElement.hasAttribute('param'): 
                specElement.firstChild.data = p[int(specElement.attributes.getNamedItem('param').value)]
            for Element in specElement.childNodes:   
                        if Element.nodeType==1:    
                            if specElement.hasAttribute('param'): 
                                specElement.firstChild.data = p[int(specElement.attributes.getNamedItem('param').value)]
                            for child in Element.childNodes:
                                if ((child.nodeType==1) and child.hasAttribute('param')):
                                    child.firstChild.data = p[int(child.attributes.getNamedItem('param').value)]                     
    
    f = open(SPEC_PATH + specXMLFileName, 'w')
    xmldoc.writexml(f)
    
def read_all(specXMLFileName, p_in = []):

    p=np.zeros(14)
    Optics=np.array([])
    Beams=np.array([])
    Cameras=np.array([])
    
    xmldoc = minidom.parse('spectrographs/' + specXMLFileName)
    
    
#    optElements = xmldoc.getElementsByTagName('optical_element') 
    #print itemlist[0].attributes['name'].value
    
    spectrograph=xmldoc.childNodes[0]
    for specElement in spectrograph.childNodes:
        if specElement.nodeType==1:
            
            #pulls p values from the spectrograph level
            if specElement.hasAttribute('param'): 
                if p_in==[]:
                    p[int(specElement.attributes.getNamedItem('param').value)] = float(specElement.firstChild.data)
                else:
                    p[int(specElement.attributes.getNamedItem('param').value)] = p_in[int(specElement.attributes.getNamedItem('param').value)]
            
            
            
            #Explores the first level child nodes (beams, optical elements, cameras)          
            if specElement.nodeName=='cameras':
                for camera in specElement.childNodes:   
                    if camera.nodeType==1:
                        if camera.hasAttribute('param'):
                            if p_in==[]:
                                p[int(camera.attributes.getNamedItem('param').value)] = float(camera.firstChild.data)
                            else:
                                p[int(camera.attributes.getNamedItem('param').value)] = p_in[int(camera.attributes.getNamedItem('param').value)]                    
                        
                        for child in camera.childNodes:
                            if child.nodeType==1: 
                                if child.nodeName=='name':
                                    name=str(child.firstChild.data)     
                                elif child.nodeName=='focalLength':
                                    fLength=float(child.firstChild.data)                                     
                                elif child.nodeName=='width':
                                    width=int(child.firstChild.data)
                                elif child.nodeName=='height':
                                   height=int(child.firstChild.data)
                                elif child.nodeName=='pSize':
                                   pSize=float(child.firstChild.data)
                                elif child.nodeName=='minLambda':
                                   minLambda=float(child.firstChild.data)
                                elif child.nodeName=='maxLambda':
                                   maxLambda=float(child.firstChild.data)
                                elif child.nodeName=='distortion':
                                   distortion=float(child.firstChild.data)
                                elif child.nodeName=='distortionCenterX':
                                   distortionCenterX=float(child.firstChild.data)
                                elif child.nodeName=='distortionCenterY':
                                   distortionCenterY=float(child.firstChild.data)
                                                                                                                            
                                if child.hasAttribute('param'):
                                    if p_in==[]:
                                        p[int(child.attributes.getNamedItem('param').value)] = float(child.firstChild.data)
                                    else:
                                        p[int(child.attributes.getNamedItem('param').value)] = p_in[int(child.attributes.getNamedItem('param').value)]
                                        if child.nodeName=='name':
                                            name = p_in[int(child.attributes.getNamedItem('param').value)]                                          
                                        elif child.nodeName=='width':
                                            width = p_in[int(child.attributes.getNamedItem('param').value)]   
                                        elif child.nodeName=='height':
                                            height = p_in[int(child.attributes.getNamedItem('param').value)]   
                                        elif child.nodeName=='pSize':
                                            pSize = p_in[int(child.attributes.getNamedItem('param').value)]   
                                        elif child.nodeName=='minLambda':
                                            minLambda = p_in[int(child.attributes.getNamedItem('param').value)]   
                                        elif child.nodeName=='maxLambda':
                                            maxLambda = p_in[int(child.attributes.getNamedItem('param').value)]   
                                        elif child.nodeName=='distortion':
                                            distortion = p_in[int(child.attributes.getNamedItem('param').value)]   
                                        elif child.nodeName=='distortionCenterX':
                                            distortionCenterX = p_in[int(child.attributes.getNamedItem('param').value)]   
                                        elif child.nodeName=='distortionCenterY':
                                            distortionCenterY = p_in[int(child.attributes.getNamedItem('param').value)]   
                                
                        newCamera=[[name, fLength, width, height, pSize, minLambda, maxLambda, distortion, distortionCenterX, distortionCenterY]]
                        
                        if len(Cameras)==0:
                            Cameras = np.array(newCamera)
                        else:
                            Cameras = np.append(Cameras, newCamera, 0)
                            
            
            elif specElement.nodeName=='beams':
                for beam in specElement.childNodes:   
                    if beam.nodeType==1:
                        if beam.hasAttribute('param'):
                            if p_in==[]:
                                p[int(beam.attributes.getNamedItem('param').value)] = float(beam.firstChild.data)
                            else:
                                p[int(beam.attributes.getNamedItem('param').value)] = p_in[int(beam.attributes.getNamedItem('param').value)]                    
                        for child in beam.childNodes:
                            if child.nodeType==1: 
                                if child.nodeName=='phi':
                                    phi=child.firstChild.data                                          
                                elif child.nodeName=='theta':
                                    theta=child.firstChild.data
                                                                                                                            
                                if child.hasAttribute('param'):
                                    if p_in==[]:
                                        p[int(child.attributes.getNamedItem('param').value)] = float(child.firstChild.data)
                                    else:
                                        p[int(child.attributes.getNamedItem('param').value)] = p_in[int(child.attributes.getNamedItem('param').value)]
                                        if child.nodeName=='phi':
                                            phi = p_in[int(child.attributes.getNamedItem('param').value)]                                          
                                        elif child.nodeName=='theta':
                                            theta = p_in[int(child.attributes.getNamedItem('param').value)]   
                                
                        newBeam=[[np.cos(np.radians(float(phi)))*np.sin(np.radians(float(theta))),np.sin(np.radians(float(phi)))*np.sin(np.radians(float(theta))),np.cos(np.radians(float(theta)))]]
                        
                        if len(Beams)==0:
                            Beams = np.array(newBeam)
                        else:
                            Beams = np.append(Beams, newBeam, 0)
                            
            elif specElement.nodeName=='optical_elements':  
                for optElement in specElement.childNodes:   
                    if optElement.nodeType==1:    
                        if optElement.hasAttribute('type'): 
                            optType=optElement.attributes.getNamedItem('type').value
                            if optType=='boundary':
                                for child in optElement.childNodes:
                                    if child.nodeType==1: 
                                        if child.nodeName=='medium':
                                            medium=child.firstChild.data
                                        elif child.nodeName=='phi':
                                            phi=child.firstChild.data                                          
                                        elif child.nodeName=='theta':
                                            theta=child.firstChild.data
                                                                                                                                    
                                        if child.hasAttribute('param'):
                                            if p_in==[]:
                                                p[int(child.attributes.getNamedItem('param').value)] = float(child.firstChild.data)
                                            else:
                                                p[int(child.attributes.getNamedItem('param').value)] = p_in[int(child.attributes.getNamedItem('param').value)]
                                                if child.nodeName=='medium':
                                                    medium=p_in[int(child.attributes.getNamedItem('param').value)]
                                                elif child.nodeName=='phi':
                                                    phi=p_in[int(child.attributes.getNamedItem('param').value)]                                            
                                                elif child.nodeName=='theta':
                                                    theta=p_in[int(child.attributes.getNamedItem('param').value)]

                                
                                #Build normal vector and update Optics array
                                n=np.array([np.cos(np.radians(float(phi)))*np.sin(np.radians(float(theta))),np.sin(np.radians(float(phi)))*np.sin(np.radians(float(theta))),np.cos(np.radians(float(theta)))])

                                newOptics = [[n,[0,0,0],OpticsBoundary,medium,0]]

                            elif optType=='RGrating':
                                for child in optElement.childNodes:
                                    if child.nodeType==1: 
                                        if child.nodeName=='medium':
                                            medium=child.firstChild.data
                                        elif child.nodeName=='phi':
                                            phi=child.firstChild.data                                         
                                        elif child.nodeName=='theta':
                                            stheta=child.firstChild.data
                                            theta=child.firstChild.data
                                        elif child.nodeName=='alpha':
                                            alpha=child.firstChild.data
                                        elif child.nodeName=='bl_period':
                                            bl_period=float(child.firstChild.data)
                                                                                                                                    
                                        if child.hasAttribute('param'):
                                            if p_in==[]:
                                                p[int(child.attributes.getNamedItem('param').value)] = float(child.firstChild.data)
                                            else:
                                                p[int(child.attributes.getNamedItem('param').value)] = p_in[int(child.attributes.getNamedItem('param').value)]
                                                if child.nodeName=='medium':
                                                    medium=p_in[int(child.attributes.getNamedItem('param').value)]
                                                elif child.nodeName=='phi':
                                                    phi=p_in[int(child.attributes.getNamedItem('param').value)]                                            
                                                elif child.nodeName=='theta':
                                                    stheta=theta=p_in[int(child.attributes.getNamedItem('param').value)]
                                                elif child.nodeName=='alpha':
                                                    alpha=p_in[int(child.attributes.getNamedItem('param').value)]
                                                elif child.nodeName=='bl_period':
                                                    bl_period=p_in[int(child.attributes.getNamedItem('param').value)]
                                               
                                #Build normal and grating vector and update Optics array
                                s = np.array([np.cos(np.radians(float(phi)))*np.sin(np.radians(float(theta))),np.sin(np.radians(float(phi)))*np.sin(np.radians(float(theta))),np.cos(np.radians(float(theta)))]) #component perp to grooves   
                                #Now find two vectors (a and b) perpendicular to s:
                                a = np.array([s[1]/np.sqrt(s[0]**2 + s[1]**2), -s[0]/np.sqrt(s[0]**2 + s[1]**2), 0])
                                b = np.cross(a,s)
                                #Create l from given alpha using a and b as basis
                                l = np.cos(np.radians(float(alpha)))*a + np.sin(np.radians(float(alpha)))*b #component along grooves
                                
                                newOptics = [[s,l,OpticsRGrating,'air',bl_period]]
                                    
                            if len(Optics)==0:
                                Optics = np.array(newOptics)
                            else:
                                Optics = np.append(Optics, newOptics, 0)



    return Beams, Optics, Cameras, p, np.radians(float(stheta))