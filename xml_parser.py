from xml.dom import minidom
import numpy as np
import constants

def read_p(specFileName='defaul_spec.xml'):
    
    p=np.zeros(11)
    
    
    xmldoc = minidom.parse('spectrographs/'+specFileName)
    
    
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


def read_Optics(specFileName='defaul_spec.xml'):
    p=np.zeros(11)
    
    Optics=np.zeros(11)
    Beams=np.zeros(100)
    
    xmldoc = minidom.parse('spectrographs/'+specFileName)
    
    
#    optElements = xmldoc.getElementsByTagName('optical_element') 
    #print itemlist[0].attributes['name'].value
    
    spectrograph=xmldoc.childNodes[0]
    for specElement in spectrograph.childNodes:
        if specElement.nodeType==1:
            
            if specElement.hasAttribute('param'): 
                p[int(specElement.attributes.getNamedItem('param').value)] = float(specElement.firstChild.data)
            
            if specElement.nodeName=='focal_length':
                fLength=float(specElement.firstChild.data)
            elif specElement.nodeName=='beams':
                currentElement=Beams
            elif specElement.nodeName=='optical_elements':  
                for Element in specElement.childNodes:   
                    if Element.nodeType==1:    
                        if specElement.hasAttribute('type'): 
                            optType=specElement.attributes.getNamedItem('type').value
                            if optType=='Boundary':
                                [n4,[0,0,0],OpticsPrism,1,0]
                        for child in Element.childNodes:
                            if ((child.nodeType==1) and child.hasAttribute('param')):
                                p[int(child.attributes.getNamedItem('param').value)] = float(child.firstChild.data)

    return Optics