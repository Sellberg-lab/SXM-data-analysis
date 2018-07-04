"""
==============================Analyze cross image=============================
Opens image with crosses indicating viral particles in x-ray images.
Extracts data about size, position etc. Uses cross2Ellipse.
Plots retults together with original image and crosses.

Alexander Woxström, 2018-06-26
#==============================================================================
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

from skimage import data
#from skimage.filters import threshold_otsu
#from skimage.segmentation import clear_border
from skimage.measure import label, regionprops
from skimage.morphology import closing, square
from skimage.color import label2rgb
from cross2Ellipse import cross2Ellipse

from PIL import Image

"""
----------------------------------Import images--------------------------------
"""
crossImg = ""
"""
---------------------------------Use cross2Ellipse-----------------------------
"""
a = cross2Ellipse(crossImg)



"""
------------------------Sort marked objects according to size------------------
"""
N = len(lengths)
pixelSize = 19.3

"""
-------------------------------SET PARAMETERS HERE-----------------------------
"""
#Length in nm/(nm/px) = px
minLength = 750/pixelSize
maxLength = 1500/pizelSize
#Ellipticity = length/width
minEllipticity = 1.2 


"""
-----------------Display result together with original image-------------------
"""





