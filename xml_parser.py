from xml.dom import minidom
import numpy as np
from constants import *

def read_p(specFileName='default_spec.xml'):
    
    p=np.zeros(11)
    
    
    xmldoc = minidom.parse(SPEC_DIR + specFileName)
    
    
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

def write_p(p, specFileName='defaul_spec.xml'):
    
    p=np.zeros(11)
    
    #todo backup prev_xml_file
    xmldoc = minidom.parse('spectrographs/'+specFileName)
    
    
#    optElements = xmldoc.getElementsByTagName('optical_element') 
    #print itemlist[0].attributes['name'].value
    
#    spectrograph=xmldoc.childNodes[0]
#    for specElement in spectrograph.childNodes:
#        if specElement.nodeType==1:
#            if specElement.hasAttribute('param'): 
#                p[int(specElement.attributes.getNamedItem('param').value)] = float(specElement.firstChild.data)
#            for Element in specElement.childNodes:   
#                        if Element.nodeType==1:    
#                            if specElement.hasAttribute('param'): 
#                                p[int(specElement.attributes.getNamedItem('param').value)] = float(specElement.firstChild.data)
#                            for child in Element.childNodes:
#                                if ((child.nodeType==1) and child.hasAttribute('param')):
#                                    p[int(child.attributes.getNamedItem('param').value)] = float(child.firstChild.data)
#
#    return p

def read_all(specFileName='default_spec.xml'):
    
    p=np.zeros(11)
    Optics=np.array([])
    Beams=np.array([])
    
    xmldoc = minidom.parse('spectrographs/'+specFileName)
    
    
#    optElements = xmldoc.getElementsByTagName('optical_element') 
    #print itemlist[0].attributes['name'].value
    
    spectrograph=xmldoc.childNodes[0]
    for specElement in spectrograph.childNodes:
        if specElement.nodeType==1:
            
            #pulls p values from the spectrograph level
            if specElement.hasAttribute('param'): 
                p[int(specElement.attributes.getNamedItem('param').value)] = float(specElement.firstChild.data)
            
            #Explores the first level child nodes
            if specElement.nodeName=='focal_length':
                fLength=float(specElement.firstChild.data)
            
            elif specElement.nodeName=='beams':
                for beam in specElement.childNodes:   
                    if beam.nodeType==1:
                        if beam.hasAttribute('param'):
                            p[int(beam.attributes.getNamedItem('param').value)] = float(beam.firstChild.data)
                        for child in beam.childNodes:
                            if child.nodeType==1: 
                                if child.nodeName=='phi':
                                    phi=np.radians(float(child.firstChild.data))                                            
                                elif child.nodeName=='theta':
                                    theta=np.radians(float(child.firstChild.data))
                                                                                                                            
                                if child.hasAttribute('param'):
                                    p[int(child.attributes.getNamedItem('param').value)] = float(child.firstChild.data)
                        
                        newBeam=[np.cos(phi)*np.sin(theta),np.sin(phi)*np.sin(theta),np.cos(theta)]
                        
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
                                            phi=np.radians(float(child.firstChild.data))                                            
                                        elif child.nodeName=='theta':
                                            theta=np.radians(float(child.firstChild.data))
                                                                                                                                    
                                        if child.hasAttribute('param'):
                                            p[int(child.attributes.getNamedItem('param').value)] = float(child.firstChild.data)
                                
                                #Build normal vector and update Optics array
                                n=np.array([np.cos(phi)*np.sin(theta),np.sin(phi)*np.sin(theta),np.cos(theta)])

                                newOptics = [[n,[0,0,0],OpticsBoundary,medium,0]]

                            elif optType=='RGrating':
                                for child in optElement.childNodes:
                                    if child.nodeType==1: 
                                        if child.nodeName=='medium':
                                            medium=child.firstChild.data
                                        elif child.nodeName=='phi':
                                            phi=np.radians(float(child.firstChild.data))                                            
                                        elif child.nodeName=='theta':
                                            stheta=theta=np.radians(float(child.firstChild.data))  #todo
                                        elif child.nodeName=='alpha':
                                            alpha=np.radians(float(child.firstChild.data))
                                        elif child.nodeName=='bl_period':
                                            bl_period=float(child.firstChild.data)
                                                                                                                                    
                                        if child.hasAttribute('param'):
                                            p[int(child.attributes.getNamedItem('param').value)] = float(child.firstChild.data)
                                
                                #Build normal and grating vector and update Optics array
                                s = np.array([np.cos(phi)*np.sin(theta),np.sin(phi)*np.sin(theta),np.cos(theta)]) #component perp to grooves   
                                #Now find two vectors (a and b) perpendicular to s:
                                a = np.array([s[1]/np.sqrt(s[0]**2 + s[1]**2), -s[0]/np.sqrt(s[0]**2 + s[1]**2), 0])
                                b = np.cross(a,s)
                                #Create l from given alpha using a and b as basis
                                l = np.cos(alpha)*a + np.sin(alpha)*b #component along grooves
                                
                                newOptics = [[s,l,OpticsRGrating,'air',bl_period]]
                                    
                            if len(Optics)==0:
                                Optics = np.array(newOptics)
                            else:
                                Optics = np.append(Optics, newOptics, 0)



    return Optics, Beams, fLength, p, stheta