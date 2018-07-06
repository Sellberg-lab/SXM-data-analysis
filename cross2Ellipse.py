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
    
    im = Image.open(crossImage)    
    image = np.array(im.getdata()).reshape(im.size[1], im.size[0], 4)[:,:,3]
    
    threshold = 0.1
    cleared = image > threshold*image.max()
    
    # label image regions
    label_image = label(cleared)
    #image_label_overlay = label2rgb(label_image, image=image)
    #fig, ax = plt.subplots(figsize=(10, 6))
    #ax.imshow(image_label_overlay, interpolation='nearest')
    #plt.colorbar()
    
    infoReg = []
    numTimes = 0
    for region in regionprops(label_image):
        #plt.imshow(region.image)
        # take regions with large enough areas
        
        if region.area >= 10:
            # draw rectangle around segmented coins
            #minr, minc, maxr, maxc = region.bbox
#            rect = mpatches.Rectangle((minc, minr), maxc - minc, maxr - minr,
#                                      fill=False, edgecolor='red', linewidth=2)
#            ax.add_patch(rect)
            #if numTimes == 1:
               #plt.imshow(region.image)
               #break
            numTimes +=1
            if numTimes == 1:
                print(numTimes, "/", label_image.max())
            if numTimes%200 == 0:
                print(numTimes, "/", label_image.max())
            infoReg.append(crossInfo(region))
            
            
    N = len(infoReg)
    return infoReg
   # return [lengths, widths, xPositions, yPositions, angles]

#im = imread("/Users/alexanderwoxstrom/Forskningsprojektet Rays/sofie+alex/cells/20180110_amoeba22-01.tif", as_gray=False)
#infoReg = cross2Ellipse("/Users/alexanderwoxstrom/Forskningsprojektet Rays/sofie+alex/cells/20180110_amoeba22-01.tif")

