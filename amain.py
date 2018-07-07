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
from analyzeCrossImage import analyzeCrossImage
from PIL import Image

infTimes = [0, 6, 12, 18, 24] # in hours, B4 (healthy) is set to -12 h
pxSize = 19.3 #19.3nm/px

#Run analyzeCrossImage on all images at T0 and collect the data needed for analysis
print("Starting to analyze T0")     
lengthsT0 = []
concT0 = []
nVirusT0 = []
fname6 = "/Users/alexanderwoxstrom/Forskningsprojektet Rays/sofie+alex/cells/" + "20171207_amoeba2-01.tif"
fname7 = "/Users/alexanderwoxstrom/Forskningsprojektet Rays/sofie+alex/cells/" + "20171207_amoeba8-01.tif"
fname8 = "/Users/alexanderwoxstrom/Forskningsprojektet Rays/sofie+alex/cells/" + "20171207_amoeba12b-01.tif"
fname9 = "/Users/alexanderwoxstrom/Forskningsprojektet Rays/sofie+alex/cells/" + "20171207_amoeba12c-01.tif"
fname10 = "/Users/alexanderwoxstrom/Forskningsprojektet Rays/sofie+alex/cells/" + "20171207_amoeba12b-01.tif"
fnamesT0 = [fname6,fname7,fname8,fname9,fname10]
for cell in fnamesT0:
     [nVirus, conc, allLengths] = analyzeCrossImage(cell)
     lengthsT0.append(allLengths)
     concT0.append(conc)
     nVirusT0.append(nVirus)
lengthsT0 = np.array(lengthsT0)
concT0 = np.array(concT0)
nVirusT0 = np.array(nVirusT0)

#Run analyzeCrossImage on all images at T6 and collect the data needed for analysis
print("Starting to analyze T6")     
lengthsT6 = []
concT6 = []
nVirusT6 = []
fname11 = "/Users/alexanderwoxstrom/Forskningsprojektet Rays/sofie+alex/cells/" + "20171211_amoeba6-01.tif"
fname12 = "/Users/alexanderwoxstrom/Forskningsprojektet Rays/sofie+alex/cells/" + "20180117_amoeba2-01.tif"
fname13 = "/Users/alexanderwoxstrom/Forskningsprojektet Rays/sofie+alex/cells/" + "20180117_amoeba3-01.tif"
fname14 = "/Users/alexanderwoxstrom/Forskningsprojektet Rays/sofie+alex/cells/" + "20180117_amoeba5-01.tif"
fname15 = "/Users/alexanderwoxstrom/Forskningsprojektet Rays/sofie+alex/cells/" + "20180117_amoeba6-01.tif"
fnamesT6 = [fname11,fname12,fname13,fname14,fname15]
for cell in fnamesT6:
     [nVirus, conc, allLengths] = analyzeCrossImage(cell)
     lengthsT6.append(allLengths)
     concT6.append(conc)
     nVirusT6.append(nVirus)
lengthsT6 = np.array(lengthsT6)
concT6 = np.array(concT6)
nVirusT6 = np.array(nVirusT6)

#Run analyzeCrossImage on all images at T12 and collect the data needed for analysis
print("Starting to analyze T12")     
lengthsT12 = []
concT12 = []
nVirusT12 = []
fname16 = "/Users/alexanderwoxstrom/Forskningsprojektet Rays/sofie+alex/cells/" + "20180110_amoeba22-01.tif"
fname17 = "/Users/alexanderwoxstrom/Forskningsprojektet Rays/sofie+alex/cells/" + "20180110_amoeba24-01.tif"
fname18 = "/Users/alexanderwoxstrom/Forskningsprojektet Rays/sofie+alex/cells/" + "20180110_amoeba25-01.tif"
fname19 = "/Users/alexanderwoxstrom/Forskningsprojektet Rays/sofie+alex/cells/" + "20180125_amoeba9-01.tif"
fnamesT12 = [fname16,fname17,fname18,fname19]
for cell in fnamesT12:
     [nVirus, conc, allLengths] = analyzeCrossImage(cell)
     lengthsT12.append(allLengths)
     concT12.append(conc)
     nVirusT12.append(nVirus)
lengthsT12 = np.array(lengthsT12)
concT12 = np.array(concT12)
nVirusT12 = np.array(nVirusT12)

#Run analyzeCrossImage on all images at T18 and collect the data needed for analysis
print("Starting to analyze T18")     
lengthsT18 = []
concT18 = []
nVirusT18 = []
fname21 = "/Users/alexanderwoxstrom/Forskningsprojektet Rays/sofie+alex/cells/" + "20171212_amoeba8-01.tif"
fname22 = "/Users/alexanderwoxstrom/Forskningsprojektet Rays/sofie+alex/cells/" + "20171212_amoeba9-01.tif"
fname23 = "/Users/alexanderwoxstrom/Forskningsprojektet Rays/sofie+alex/cells/" + "20171221_amoeba4-01.tif"
fname24 = "/Users/alexanderwoxstrom/Forskningsprojektet Rays/sofie+alex/cells/" + "20171221_amoeba7a-01.tif"
fname25 = "/Users/alexanderwoxstrom/Forskningsprojektet Rays/sofie+alex/cells/" + "20171221_amoeba7b-01.tif"
fnamesT18 = [fname21,fname22,fname23,fname24,fname25]
for cell in fnamesT18:
     [nVirus, conc, allLengths] = analyzeCrossImage(cell)
     lengthsT18.append(allLengths)
     concT18.append(conc)
     nVirusT18.append(nVirus)
lengthsT18 = np.array(lengthsT18)
concT18 = np.array(concT18)
nVirusT18 = np.array(nVirusT18)

#Run analyzeCrossImage on all images at T24 and collect the data needed for analysis
print("Starting to analyze T24")     
lengthsT24 = []
concT24 = []
nVirusT24 = []
fname26 = "/Users/alexanderwoxstrom/Forskningsprojektet Rays/sofie+alex/cells/" + "20171114_amoeba3-01.tif"
fname27 = "/Users/alexanderwoxstrom/Forskningsprojektet Rays/sofie+alex/cells/" + "20171114_amoeba12-01.tif"
fname28 = "/Users/alexanderwoxstrom/Forskningsprojektet Rays/sofie+alex/cells/" + "20171208_amoeba8-01.tif"
fnamesT24 = [fname26,fname27,fname28]
for cell in fnamesT24:
     [nVirus, conc, allLengths] = analyzeCrossImage(cell)
     lengthsT24.append(allLengths)
     concT24.append(conc)
     nVirusT24.append(nVirus)
lengthsT24 = np.array(lengthsT24)
concT24 = np.array(concT24)
nVirusT24 = np.array(nVirusT24)

print("Analysis Done!")


#Create the data that needs to be plotted 
#nVirus plots:
nVirusVsTimeMeanPlusStd =  [nVirusT0.mean()+nVirusT0.std(),nVirusT6.mean()+nVirusT6.std(), nVirusT12.mean()+nVirusT12.std(),nVirusT18.mean()+nVirusT18.std(), nVirusT24.mean()+nVirusT24.std()]
nVirusVsTimeMeanMinusStd = [nVirusT0.mean()-nVirusT0.std(),nVirusT6.mean()-nVirusT6.std(), nVirusT12.mean()-nVirusT12.std(),nVirusT18.mean()-nVirusT18.std(), nVirusT24.mean()-nVirusT24.std()]
nVirusVsTimeMean = [nVirusT0.mean(),nVirusT6.mean(), nVirusT12.mean(),nVirusT18.mean(), nVirusT24.mean()]
nVirusVsTimeStd = [nVirusT0.std(), nVirusT6.std(), nVirusT12.std(), nVirusT18.std(), nVirusT24.std()]

#conc plots:
concVsTimeMeanPlusStd =  [concT0.mean()+concT0.std(),concT6.mean()+concT6.std(), concT12.mean()+concT12.std(),concT18.mean()+concT18.std(), concT24.mean()+concT24.std()]
concVsTimeMeanMinusStd =  [concT0.mean()-concT0.std(),concT6.mean()-concT6.std(), concT12.mean()-concT12.std(),concT18.mean()-concT18.std(), concT24.mean()-concT24.std()]
concVsTimeMean = [concT0.mean(), concT6.mean(), concT12.mean(), concT18.mean(), concT24.mean()]
concVsTimeStd = [concT0.std(), concT6.std(), concT12.std(), concT18.std(), concT24.std()]

#Calculate Hist - T0
lMin = 50 # nm
lMax = 1650 # nm
dL = 100 # nm
plt.show()
hist_bins = np.linspace(lMin, lMax, num=np.int((lMax - lMin)/dL + 1))
hist, hist_bins = np.histogram(np.concatenate(lengthsT0)*pxSize, bins=hist_bins)
hist_bins_center = np.array([(hist_bins[i] + hist_bins[i+1])/2 for i in range(len(hist_bins) - 1)])
plt.bar(hist_bins_center, hist, 0.8*dL)
plt.xlabel('Different sizes (nm)')
plt.ylabel('Number of particles')
plt.savefig('/Users/alexanderwoxstrom/Forskningsprojektet Rays/alex programmering/Plots/histT0')
plt.show()

#Calculate Hist - T12
lMin = 50 # nm
lMax = 1650 # nm
dL = 100 # nm
plt.show()
hist_bins = np.linspace(lMin, lMax, num=np.int((lMax - lMin)/dL + 1))
hist, hist_bins = np.histogram(np.concatenate(lengthsT12)*pxSize, bins=hist_bins)
hist_bins_center = np.array([(hist_bins[i] + hist_bins[i+1])/2 for i in range(len(hist_bins) - 1)])
plt.bar(hist_bins_center, hist, 0.8*dL)
plt.xlabel('Different sizes (nm)')
plt.ylabel('Number of particles')
plt.savefig('/Users/alexanderwoxstrom/Forskningsprojektet Rays/alex programmering/Plots/histT12')
plt.show()

#Calculate Hist - T24
lMin = 50 # nm
lMax = 1650 # nm
dL = 100 # nm
plt.show()
hist_bins = np.linspace(lMin, lMax, num=np.int((lMax - lMin)/dL + 1))
hist, hist_bins = np.histogram(np.concatenate(lengthsT24)*pxSize, bins=hist_bins)
hist_bins_center = np.array([(hist_bins[i] + hist_bins[i+1])/2 for i in range(len(hist_bins) - 1)])
plt.bar(hist_bins_center, hist, 0.8*dL)
plt.xlabel('Different sizes (nm)')
plt.ylabel('Number of particles')
plt.savefig('/Users/alexanderwoxstrom/Forskningsprojektet Rays/alex programmering/Plots/histT24')
plt.show()

#Plot ConcVsTimeMean with standarddeviation added, subtracted and without
plt.errorbar(infTimes, concVsTimeMean, concVsTimeStd, 0,"-ob", capsize=5)
plt.xlabel('Time (h)')
plt.ylabel('Concentration of viruses (virus$/μm^2$)')
plt.savefig('/Users/alexanderwoxstrom/Forskningsprojektet Rays/alex programmering/Plots/ConcVsTime')
plt.show()

#Plot nVirus with standarddeviation added, subtracted and without
plt.errorbar(infTimes, nVirusVsTimeMean, nVirusVsTimeStd, 0,"-ok", capsize=5)
plt.xlabel('Time (h)')
plt.ylabel('Number of viruses per cell')
plt.savefig('/Users/alexanderwoxstrom/Forskningsprojektet Rays/alex programmering/Plots/nVirusVsTime')
plt.show()
 