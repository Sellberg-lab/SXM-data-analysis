import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from crossInfo import crossInfo
from skimage import data
#from skimage.filters import threshold_otsu
#from skimage.segmentation import clear_border
from skimage.measure import label, regionprops
from skimage.morphology import closing, square
from skimage.color import label2rgb

from PIL import Image

def cross2Ellipse(crossImage):
    """
    cross2Ellipse takes an image with black crosses and interprets them as the 
    dimention of elliptical particles. A list containing the lengths, width and positions of the 
    ellipses is returned.
    """
    #Open image
    im = Image.open(crossImage)    
    image = np.array(im.getdata()).reshape(im.size[1], im.size[0], 4)[:,:,3]
    
    #Set threshold
    threshold = 0.1
    cleared = image > threshold*image.max()
    
    # label image regions
    label_image = label(cleared)
    
    #Create list for information about regions.
    infoReg = []
    numTimes = 0
    
    #Use regionprops and extract data about each region, append info to infoReg
    for region in regionprops(label_image):
        if region.area >= 10:
            numTimes +=1
            if numTimes == 1:
                print(numTimes, "/", label_image.max())
            if numTimes%200 == 0:
                print(numTimes, "/", label_image.max())
            infoReg.append(crossInfo(region))
            
    #Return infoReg = [lengths, widths, xPositions, yPositions, angles]        
    return infoReg
 