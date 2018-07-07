import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy.spatial import ConvexHull
from skimage import data
#from skimage.filters import threshold_otsu
#from skimage.segmentation import clear_border
from skimage.measure import label, regionprops
from skimage.morphology import closing, square
from skimage.color import label2rgb
from cross2Ellipse import cross2Ellipse

from PIL import Image

"""
---------------------------------Use cross2Ellipse-----------------------------
"""
def analyzeCrossImage(fname):
    """
    Take a tif file with many crosses representing structures and returns number of viruses, conc and lengths 
    """
    iname = fname[:-7] + ".jpg"
    pname = fname[:-7] + "-analyzed.png"
    i = Image.open(iname)
    
    #Run cross2Ellipse
    infoReg = cross2Ellipse(fname)
    
    #Define pixel constants needed for conversion 
    pxSize = 19.3 #19.3nm/px
    pxSqrd = 19.3**2 #nm^2 per pixel
    
    #Define what parameters a virus is categorized by
    #Length in nm/(nm/px) = px
    minLength = 900/pxSize
    maxLength = 1800/pxSize
    #Ellipticity = length/width
    minEllipticity = 1.4 

    #Create lists for all information
    vLengths = []
    allLengths = []
    vWidths = []
    allWidths = []
    allPositions = []
    allMidPoints = []
    viruses = []
    
    #Extract information from each region
    for structure in infoReg:
        allPositions.append(structure[3])
        allMidPoints.append(structure[2])
        allLengths.append(structure[0])
        allWidths.append(structure[1])

        if structure[0] > minLength:
            if structure[0] < maxLength:
                if structure[0]/structure[1] > minEllipticity:
                    viruses.append(structure)
                    vLengths.append(structure[0])
                    vWidths.append(structure[1])

    #Extract information about each virus
    vMidPoints = []
    vPositions = []
    for virus in viruses:
        vPositions.append(virus[3])
        vMidPoints.append(virus[2])
    
    vPoints = np.array(vPositions)
    allPoints = np.array(allPositions)
    
    plt.figure(figsize=(6,6), dpi=100)
    plt.imshow(i)
    
    #Create a ConvexHull which can be used to calculate area
    aAmoeba = ConvexHull(allPoints)
    nVirus = len(viruses)
    plt.plot(allPoints[:,0], allPoints[:,1], '.')
    for simplex in aAmoeba.simplices:
        plt.plot(allPoints[simplex, 0], allPoints[simplex, 1], 'k-')
    if nVirus > 2:  
        hull = ConvexHull(vPoints)
        plt.plot(vPoints[:,0], vPoints[:,1], 'r.')
        for simplex in hull.simplices:
            plt.plot(vPoints[simplex, 0], vPoints[simplex, 1], 'r-')

    plt.savefig(pname)
    print('Saved png to: %s' % pname)
    
    #Define area in px and sq μm    
    aPx = aAmoeba.volume
    aMikroM = aPx *pxSqrd/1000000 #Scales 19.3 each axis. Turns nm^2 to μm^2
    
    #Return values
    nVirus = len(viruses)
    conc = nVirus/aMikroM
    return nVirus, conc, allLengths 