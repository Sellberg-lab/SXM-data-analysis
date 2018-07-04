import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
#from crossInfo import crossInfo
from skimage import data
from skimage.io import imread
from skimage.transform import radon, rescale

#from skimage.filters import threshold_otsu
#from skimage.segmentation import clear_border
from skimage.measure import label, regionprops
from skimage.morphology import closing, square
from skimage.color import label2rgb

from PIL import Image

def crossInfo(region):
    """
    crossInfo takes the logical image of a cross and analyzes its length, width,
    angle and position in the xy-plane. 
    """
    image = region.image
    scale= 4
    image = rescale(image, scale, multichannel=True)
    #img = image[:,:,0]
    lim = 0.6*image.max()
    img = image[:,:]
    image = (img > lim).astype(np.int)

    
    #fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8, 4.5))
    #ax1.set_title("Original")
    #ax1.imshow(image, cmap=plt.cm.Greys_r)
    
    theta = np.linspace(0., 180., max(image.shape), endpoint=False)
    radIm = radon(image, theta=theta, circle=False)
    indMax = radIm.argmax()
    
    lenRadIm = len(radIm) 
    radMaxPos = np.unravel_index(radIm.argmax(), radIm.shape)
    
    dTheta = theta[1]-theta[0]
    ang1 = (dTheta)*radMaxPos[1]
    if ang1 <= 90:
        ang2 = ang1 + 90
    else:
        ang2 = ang1 - 90
    
    length = radIm[radMaxPos]/scale
    width = radIm[:, np.round(ang2/dTheta).astype(np.int)].max()/scale
    #b = radon(image, theta=[ang1], circle=False)
    #length2 = b.max()*pxSize/scale
    b = radon(image, theta=[ang2], circle=False)
    width2 = b.max()/scale
    
    xPos = image.shape[1]/2
    yPos = image.shape[0]/2
    #ax2.set_title("Radon transform\n(RadonIm)")
    #ax2.set_xlabel("Projection angle (deg)")
    #ax2.set_ylabel("Projection position (pixels)")
    #ax2.imshow(radIm, cmap=plt.cm.Greys_r,
    #           extent=(0, 180, 0, radIm.shape[0]), aspect='auto')
    
    #fig.tight_layout()
    #plt.show()
    print(width, width2)
    print ([length, width, xPos,yPos, ang1])
    return [length, width, xPos,yPos, ang1]

