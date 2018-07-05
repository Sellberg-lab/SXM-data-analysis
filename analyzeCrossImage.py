"""
==============================Analyze cross image=============================
Opens image with crosses indicating viral particles in x-ray images.
Extracts data about size, position etc. Uses cross2Ellipse.
Plots retults together with original image and crosses.

ans = 168

Alexander Woxström, 2018-06-26
#==============================================================================
"""
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
#fname = input("Enter filename and position :")
fname = "/Users/alexanderwoxstrom/Forskningsprojektet Rays/sofie+alex/cells/20180110_amoeba22-01.tif"
iname = fname[:-7] + ".jpg"
i = Image.open(iname)
infoReg = cross2Ellipse(fname)

# Problematiska bilder:
# /Users/alexanderwoxstrom/Forskningsprojektet Rays/sofie+alex/cells/20171025_amoeba5-01.tif

"""
------------------------Sort marked objects according to size------------------
"""
pxSize = 19.3 #19.3nm/px
pxSqrd = 19.3**2 #nm^2 per pixel

"""
-------------------------------SET PARAMETERS HERE-----------------------------
"""
#Length in nm/(nm/px) = px
minLength = 750/pxSize
maxLength = 1500/pxSize
#Ellipticity = length/width
minEllipticity = 1.2 
"""
-----------------Display result together with original image-------------------
"""
lengths = []
ordLengths = sorted(lengths)
widths = []
ordWidths = sorted(widths)

allPositions = []
allMidPoints = []
viruses = []
for structure in infoReg:
    allPositions.append(structure[3])
    allMidPoints.append(structure[2])

    if structure[0] > minLength:
        if structure[0] < maxLength:
            if structure[0]/structure[1] > minEllipticity:
                viruses.append(structure)
                lengths.append(structure[0])
                widths.append(structure[1])
          
nVirus = len(viruses)

vMidPoints = []
vPositions = []
for virus in viruses:
    vPositions.append(virus[3])
    vMidPoints.append(virus[2])

vPoints = np.array(vPositions)
allPoints = np.array(allPositions)

plt.imshow(i)
aAmoeba = ConvexHull(allPoints)
#plt.plot(allPoints[:,0], allPoints[:,1], 'o')
#plt.axis([2049,0,0,2049])
for simplex in aAmoeba.simplices:
    plt.plot(allPoints[simplex, 0], allPoints[simplex, 1], 'k-')

hull = ConvexHull(vPoints)
plt.plot(vPoints[:,0], vPoints[:,1], 'ro')
#plt.axis([2049,0,0,2049])
for simplex in hull.simplices:
    plt.plot(vPoints[simplex, 0], vPoints[simplex, 1], 'r-')

plt.show()

aPx = aAmoeba.volume
aMikroM = aPx *pxSqrd/1000000 #Scales 19.3 each axis. Turns nm^2 to μm^2

conc = nVirus/aMikroM
print(aPx)
print(aMikroM)
print(conc)


