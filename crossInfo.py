import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
#from crossInfo import crossInfo
from skimage import data
from skimage.io import imread
from skimage.transform import radon, rescale, resize

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
    
    #Original image of region 
    orgImage = region.image
    
    #Set scale, OBS: Exponentially slows program
    scale = 2.5
    image = rescale(orgImage, scale, multichannel=False)
   
    #Binarize image
    lim = 0.6*image.max()
    img = image[:,:]
    image = (img > lim).astype(np.int)
 
    #Begin Radon Transform and define parameters for it
    theta = np.linspace(0., 180., max(image.shape), endpoint=False)
    radIm = radon(image, theta=theta, circle=False)
    radMaxPos = np.unravel_index(radIm.argmax(), radIm.shape)
    
    #Set angles
    dTheta = theta[1]-theta[0]
    ang1 = (dTheta)*radMaxPos[1]
    if ang1 <= 90:
        ang2 = ang1 + 90
    else:
        ang2 = ang1 - 90
    
    #Define lengths and widths, OBS: 2 different ways. Current setup works.
    length = radIm[radMaxPos]/scale
    #width = radIm[:, np.round(ang2/dTheta).astype(np.int)].max()/scale
    #b = radon(image, theta=[ang1], circle=False)
    #length2 = b.max()*pxSize/scale
    b = radon(image, theta=[ang2], circle=False)
    width = b.max()/scale
    
    #Define the position of the region
    minr, minc, maxr, maxc = region.bbox
    xMid = orgImage.shape[1]/2
    yMid = orgImage.shape[0]/2
    xPos = xMid + minc
    yPos = yMid + minr
    middle = [xMid, yMid]
    pos = [xPos, yPos]
    
    #Return values for region
    return [length, width, middle, pos, ang1]

 