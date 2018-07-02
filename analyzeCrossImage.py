"""
==============================Analyze cross image=============================
Opens image with crosses indicating viral particles in x-ray images.
Extracts data about size, position etc. Uses cross2Ellipse.
Plots retults together with original image and crosses.

Alexander Woxström, 2018-06-26
#==============================================================================
"""
import numpy as np
from matplotlib import pyplot as plt

"""
----------------------------------Import images--------------------------------
"""
#File1 och 2 variablerna har hela pathwayn till bilderna
#Använd pil för  att lagra tiff filen som t2


#Skapa variabel crossImg som är bilden som en två dim numpyarray med numeriska
#värden från 0-255. Uppdelad så att 
"""
---------------------------------Use cross2Ellipse-----------------------------
"""
from cross2Ellipse import cross2Ellipse
a = cross2Ellipse(crossImg)
N = 


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





