from xml.dom import minidom
import numpy as np
import os, time

import constants as c


def read_p(modelXMLFile):
    
    p=np.zeros(11)
    
    
    xmldoc = minidom.parse(modelXMLFile)
    
    
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

def write_p(p, modelXMLFile):
    
    bkp_time = time.time()
    
    newmodelXMLFile = find_modelXMLFile(modelXMLFile)
#     file_path = modelXMLFile[:-len(modelXMLFile.rsplit('/')[-1])] #take only path
    
    os_command = 'cp ' + modelXMLFile + ' ' + newmodelXMLFile
    os.system(os_command)
    
    modelXMLFile = newmodelXMLFile
    
    xmldoc = minidom.parse(modelXMLFile)
    
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
    
    f = open(modelXMLFile, 'w')
    xmldoc.writexml(f)
    print "Wrote new model in " + modelXMLFile
    return modelXMLFile

def read_all(modelXMLFile, p_try = []):

    p=np.zeros(16)
    Optics=np.array([])
    Beams=np.array([])
    Cameras=np.array([])
    
    xmldoc = minidom.parse(modelXMLFile)
    
    
#    optElements = xmldoc.getElementsByTagName('optical_element') 
    #print itemlist[0].attributes['name'].value
    
    spectrograph=xmldoc.childNodes[0]
    for specElement in spectrograph.childNodes:
        if specElement.nodeType==1:
            
            #pulls p values from the spectrograph level
            if specElement.hasAttribute('param'): 
                if p_try==[]:
                    p[int(specElement.attributes.getNamedItem('param').value)] = float(specElement.firstChild.data)
                else:
                    p[int(specElement.attributes.getNamedItem('param').value)] = p_try[int(specElement.attributes.getNamedItem('param').value)]
            
            
            
            #Explores the first level child nodes (beams, optical elements, cameras)          
            if specElement.nodeName=='cameras':
                for camera in specElement.childNodes:   
                    if camera.nodeType==1:
                        if camera.hasAttribute('param'):
                            if p_try==[]:
                                p[int(camera.attributes.getNamedItem('param').value)] = float(camera.firstChild.data)
                            else:
                                p[int(camera.attributes.getNamedItem('param').value)] = p_try[int(camera.attributes.getNamedItem('param').value)]                    
                        
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
                                elif child.nodeName=='distortion1':
                                   distortion1=float(child.firstChild.data)
                                elif child.nodeName=='distortion2':
                                   distortion2=float(child.firstChild.data)
                                elif child.nodeName=='distortion3':
                                   distortion3=float(child.firstChild.data)
                                elif child.nodeName=='distortionCenterX':
                                   distortionCenterX=float(child.firstChild.data)
                                elif child.nodeName=='distortionCenterY':
                                   distortionCenterY=float(child.firstChild.data)
                                                                                                                            
                                if child.hasAttribute('param'):
                                    if p_try==[]:
                                        p[int(child.attributes.getNamedItem('param').value)] = float(child.firstChild.data)
                                    else:
                                        p[int(child.attributes.getNamedItem('param').value)] = p_try[int(child.attributes.getNamedItem('param').value)]
                                        if child.nodeName=='name':
                                            name = p_try[int(child.attributes.getNamedItem('param').value)]                                          
                                        elif child.nodeName=='width':
                                            width = p_try[int(child.attributes.getNamedItem('param').value)]   
                                        elif child.nodeName=='height':
                                            height = p_try[int(child.attributes.getNamedItem('param').value)]   
                                        elif child.nodeName=='pSize':
                                            pSize = p_try[int(child.attributes.getNamedItem('param').value)]   
                                        elif child.nodeName=='minLambda':
                                            minLambda = p_try[int(child.attributes.getNamedItem('param').value)]   
                                        elif child.nodeName=='maxLambda':
                                            maxLambda = p_try[int(child.attributes.getNamedItem('param').value)]   
                                        elif child.nodeName=='distortion1':
                                            distortion1 = p_try[int(child.attributes.getNamedItem('param').value)]   
                                        elif child.nodeName=='distortion2':
                                            distortion2 = p_try[int(child.attributes.getNamedItem('param').value)]   
                                        elif child.nodeName=='distortion3':
                                            distortion3 = p_try[int(child.attributes.getNamedItem('param').value)]   
                                        elif child.nodeName=='distortionCenterX':
                                            distortionCenterX = p_try[int(child.attributes.getNamedItem('param').value)]   
                                        elif child.nodeName=='distortionCenterY':
                                            distortionCenterY = p_try[int(child.attributes.getNamedItem('param').value)]   
                                
                        newCamera=[[name, fLength, width, height, pSize, minLambda, maxLambda]]#, distortion1, distortion2, distortion3, distortionCenterX, distortionCenterY]]
                        
                        if len(Cameras)==0:
                            Cameras = np.array(newCamera)
                        else:
                            Cameras = np.append(Cameras, newCamera, 0)
                            
            
            elif specElement.nodeName=='beams':
                for beam in specElement.childNodes:   
                    if beam.nodeType==1:
                        
                        if beam.hasAttribute('param'):
                            if p_try==[]:
                                p[int(beam.attributes.getNamedItem('param').value)] = float(beam.firstChild.data)
                            else:
                                p[int(beam.attributes.getNamedItem('param').value)] = p_try[int(beam.attributes.getNamedItem('param').value)]                    
                        
                        if beam.hasAttribute('beamID'):
                            beamID = int(beam.attributes.getNamedItem('beamID').value)
                   
                        
                        for child in beam.childNodes:
                            if child.nodeType==1: 
                                if child.nodeName=='phi':
                                    phi=child.firstChild.data                                          
                                elif child.nodeName=='theta':
                                    theta=child.firstChild.data
                                                                                                                            
                                if child.hasAttribute('param'):
                                    if p_try==[]:
                                        p[int(child.attributes.getNamedItem('param').value)] = float(child.firstChild.data)
                                    else:
                                        p[int(child.attributes.getNamedItem('param').value)] = p_try[int(child.attributes.getNamedItem('param').value)]
                                        if child.nodeName=='phi':
                                            phi = p_try[int(child.attributes.getNamedItem('param').value)]                                          
                                        elif child.nodeName=='theta':
                                            theta = p_try[int(child.attributes.getNamedItem('param').value)]   
                                
                        newBeam=[np.cos(np.radians(float(phi)))*np.sin(np.radians(float(theta))),np.sin(np.radians(float(phi)))*np.sin(np.radians(float(theta))),np.cos(np.radians(float(theta)))]
                        
                        if len(Beams)==0:
                            Beams = np.array([[newBeam, beamID]])
                        else:
                            Beams = np.append(Beams, [[newBeam, beamID]], 0)
                            
            elif specElement.nodeName=='optical_elements':  
                for optElement in specElement.childNodes:   
                    if optElement.nodeType==1:    
                        if optElement.hasAttribute('type'): 
                            optType=int(optElement.attributes.getNamedItem('type').value)
                            if optType==c.opticTypes.boundary:
                                for child in optElement.childNodes:
                                    if child.nodeType==1: 
                                        if child.nodeName=='medium':
                                            medium=child.firstChild.data
                                        elif child.nodeName=='phi':
                                            phi=child.firstChild.data                                          
                                        elif child.nodeName=='theta':
                                            theta=child.firstChild.data
                                                                                                                                    
                                        if child.hasAttribute('param'):
                                            if p_try==[]:
                                                p[int(child.attributes.getNamedItem('param').value)] = float(child.firstChild.data)
                                            else:
                                                p[int(child.attributes.getNamedItem('param').value)] = p_try[int(child.attributes.getNamedItem('param').value)]
                                                if child.nodeName=='medium':
                                                    medium=p_try[int(child.attributes.getNamedItem('param').value)]
                                                elif child.nodeName=='phi':
                                                    phi=p_try[int(child.attributes.getNamedItem('param').value)]                                            
                                                elif child.nodeName=='theta':
                                                    theta=p_try[int(child.attributes.getNamedItem('param').value)]

                                
                                #Build normal vector and update Optics array
                                
#                                 n=np.array([np.cos(np.radians(float(phi)))*np.sin(np.radians(float(theta))),np.sin(np.radians(float(phi)))*np.sin(np.radians(float(theta))),np.cos(np.radians(float(theta)))])

                                newOptics = [[np.cos(np.radians(float(phi)))*np.sin(np.radians(float(theta))),
                                              np.sin(np.radians(float(phi)))*np.sin(np.radians(float(theta))),
                                              np.cos(np.radians(float(theta))),
                                              0,0,0,
                                              c.opticTypes.boundary,
                                              int(medium),
                                              0,
                                              0]]

                            elif optType==c.opticTypes.eGrating:
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
                                            if p_try==[]:
                                                p[int(child.attributes.getNamedItem('param').value)] = float(child.firstChild.data)
                                            else:
                                                p[int(child.attributes.getNamedItem('param').value)] = p_try[int(child.attributes.getNamedItem('param').value)]
                                                if child.nodeName=='medium':
                                                    medium=p_try[int(child.attributes.getNamedItem('param').value)]
                                                elif child.nodeName=='phi':
                                                    phi=p_try[int(child.attributes.getNamedItem('param').value)]                                            
                                                elif child.nodeName=='theta':
                                                    stheta=theta=p_try[int(child.attributes.getNamedItem('param').value)]
                                                elif child.nodeName=='alpha':
                                                    alpha=p_try[int(child.attributes.getNamedItem('param').value)]
                                                elif child.nodeName=='bl_period':
                                                    bl_period=p_try[int(child.attributes.getNamedItem('param').value)]
                                               
                                #Build normal and grating vector and update Optics array
                                s = np.array([np.cos(np.radians(float(phi)))*np.sin(np.radians(float(theta))),np.sin(np.radians(float(phi)))*np.sin(np.radians(float(theta))),np.cos(np.radians(float(theta)))]) #component perp to grooves   
                                #Now find two vectors (a and b) perpendicular to s:
                                a = np.array([s[1]/np.sqrt(s[0]**2 + s[1]**2), -s[0]/np.sqrt(s[0]**2 + s[1]**2), 0])
                                b = np.cross(a,s)
                                #Create l from given alpha using a and b as basis
                                l = np.cos(np.radians(float(alpha)))*a + np.sin(np.radians(float(alpha)))*b #component along grooves
                                
                                newOptics = [[s[0], s[1], s[2],
                                            l[0], l[1], l[2],
                                            c.opticTypes.eGrating,
                                            c.media.air,
                                            bl_period,
                                            0]]
                                    
                            elif optType==c.opticTypes.VPHGrating:
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
                                            if p_try==[]:
                                                p[int(child.attributes.getNamedItem('param').value)] = float(child.firstChild.data)
                                            else:
                                                p[int(child.attributes.getNamedItem('param').value)] = p_try[int(child.attributes.getNamedItem('param').value)]
                                                if child.nodeName=='medium':
                                                    medium=p_try[int(child.attributes.getNamedItem('param').value)]
                                                elif child.nodeName=='phi':
                                                    phi=p_try[int(child.attributes.getNamedItem('param').value)]                                            
                                                elif child.nodeName=='theta':
                                                    stheta=theta=p_try[int(child.attributes.getNamedItem('param').value)]
                                                elif child.nodeName=='alpha':
                                                    alpha=p_try[int(child.attributes.getNamedItem('param').value)]
                                                elif child.nodeName=='bl_period':
                                                    bl_period=p_try[int(child.attributes.getNamedItem('param').value)]
                                               
                                #Build normal and grating vector and update Optics array
                                s = np.array([np.cos(np.radians(float(phi)))*np.sin(np.radians(float(theta))),np.sin(np.radians(float(phi)))*np.sin(np.radians(float(theta))),np.cos(np.radians(float(theta)))]) #component perp to grooves   
                                #Now find two vectors (a and b) perpendicular to s:
                                a = np.array([s[1]/np.sqrt(s[0]**2 + s[1]**2), -s[0]/np.sqrt(s[0]**2 + s[1]**2), 0])
                                b = np.cross(a,s)
                                #Create l from given alpha using a and b as basis
                                l = np.cos(np.radians(float(alpha)))*a + np.sin(np.radians(float(alpha)))*b #component along grooves

#                                 newOptics = [[s,l,OpticsVPHGrating,'air',bl_period]]
                                newOptics = [[s[0], s[1], s[2],
                                            l[0], l[1], l[2],
                                            OpticsVPHGrating,
                                            c.media.air,
                                            bl_period,
                                            0]]

                            if Optics.shape[0]==0:
                                Optics = np.array(newOptics)
                            else:
                                Optics = np.append(Optics,newOptics,axis=0)
                            

    return Beams, Optics, Cameras, p, np.radians(float(stheta))

def find_modelXMLFile( modelXMLFile):
    
    specNameNoPath = modelXMLFile
    
    #find if it has a version number
    specName = specNameNoPath.rsplit('_v')[0] 
    if len(specName)==len(specNameNoPath):
        specName = specNameNoPath[:-4]
        
    a = os.listdir('.') # retrieves folder contents
    b = np.array(a) # converts to numpy
    b = b[(np.char.startswith(b, specName) & np.char.endswith(b,'xml'))] # removes non xml files and keeps only spectrograph spcecific
    
    b = np.char.rstrip(b, 'xml') 
    b = np.char.rstrip(b, '.') # removes .xml bit
    
    b = np.char.lstrip(b, specName+ '_v')
    
    prevVNumbers = np.zeros(len(b))
    
    for index in range(len(b)): # I don't know how to turn all members of the array to int in one line
        try:
            prevVNumbers[index] = int(b[index])
        except:
            prevVNumbers[index] = 0
            
    vNumber = max(prevVNumbers) + 1
    modelXMLFile = specName + '_v' + str(int(vNumber)) + '.xml'    #compose final filename 
    
    return modelXMLFile

