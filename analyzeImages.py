#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
#import matplotlib.patches as mpatches
import scipy.constants as constants
from scipy.spatial import ConvexHull
#from skimage import data
#from skimage.measure import label, regionprops
#from skimage.morphology import closing, square
#from skimage.color import label2rgb
from PIL import Image
import os,sys,time
import subprocess
import argparse
import csv
import h5py

import cross2Ellipse

source_dir = os.path.dirname(os.path.realpath(cross2Ellipse.__file__))
analysis_dir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(analysis_dir)

# Commandline arguments
def parse_cmdline_args():
    parser = argparse.ArgumentParser(description='SXM data analysis script for virus-infected cells')

    required_grp = parser.add_argument_group('required arguments')
    
    required_grp.add_argument('-o', '--output-directory', metavar='/output/directory/', required=True,
                        help="output directory", type=str)
    required_grp.add_argument('-i', '--image-number', metavar='image_number', required=True,
                        help="image number, can also be a series of images, for example: 1,3,5-10,20,22", type=str)
    parser.add_argument('-f', '--csv-file', default="sxm_data_analysis.csv", metavar='csv_file',
                        help="filename of parameter sheet (default: sxm_data_analysis.csv)", type=str)
    parser.add_argument('-y', '--image-type', metavar='image_type', required=False,
                        help="image type of the images that shall be processed (default: Cell,Virus)", type=str, default="Cell,Virus")
    parser.add_argument('-s', '--sorters', metavar='sorters', required=False,
                        help="sorters (in lowercase) who shall be processed (default: gunnar,komang)", type=str, default="gunnar,komang")
    parser.add_argument('-t', '--infection-time', metavar='infection_time', required=False,
                        help="only process single infection time (default: off)", type=int, default=None)
    parser.add_argument('-l', '--output-level', metavar='output_level',
                        help="output level defines how much data per event will be stored (default: 2; 0: no output (\"dry-image\"), 1: only scalar output (Number of particles, number of viruses, etc.), 2: scalar output and images, 3: scalar output, images and HDF5 file of processed data)", type=int, default=2)
    parser.add_argument('-r', '--re-process', default=0, metavar='re_process',
                        help="re-process already processed images according to .csv sheet (default: 0)", type=int)
    parser.add_argument('-v', '--verbose', default=0, metavar='verbose',
                        help="add additional output to screen (default: 0)", type=int)
    parser.add_argument('--sha-version-check', default=0, metavar='sha_version_check',
                        help="Perform version check for SXM data analysis git repo whenever script is executed (default: 0)", type=int)
    
    if(len(sys.argv) == 1):
        parser.print_help()
        
    return parser.parse_args()

def analyzeCrossImage(fname, output_dir, save_png=True):
    """
    Take a tif file with many crosses representing structures and returns number of particles, viruses, conc and lengths
    """
    if "komang" in fname:
        iname = fname[:-7] + ".jpg"
        pname = fname.split('/')[-1][:-7] + "-analyzed.png"
        pdir = output_dir + '/komang/'
    elif "gunnar" in fname:
        iname = fname[:-11] + ".jpg"
        pname = fname.split('/')[-1][:-11] + "-analyzed.png"
        pdir = output_dir + '/gunnar/'
    else:
        print("Unknown sorter (Gunnar or Komang), aborting..")
        sys.exit(-1)
    i = Image.open(iname)
    
    #Run cross2Ellipse
    infoReg = cross2Ellipse.cross2Ellipse(fname)
    
    #Define pixel constants needed for conversion 
    pxSize = 19.3 #19.3nm/px
    pxSqrd = 19.3**2 #nm^2 per pixel
    
    #Define what parameters a virus is categorized by
    #Length in nm/(nm/px) = px
    minLength = 700/pxSize
    maxLength = 1300/pxSize
    #Ellipticity = length/width
    minEllipticity = 1.6 

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
    
    if save_png:
        if not os.path.exists(pdir):
            os.mkdir(pdir)
            print("Created directory: %s" % pdir)
        plt.savefig(pdir + pname)
        print('Saved png to: %s' % pname)
    
    #Define area in px and sq um    
    aPx = aAmoeba.volume
    aMicroM = aPx*pxSqrd/1000000 #Scales 19.3 each axis. Turns nm^2 to um^2
    
    #Return values
    nVirus = len(viruses)
    conc = nVirus/aMicroM
    return aMicroM, allLengths, allWidths, allPositions, vLengths, vWidths, vPositions

if __name__ == "__main__":
    args = parse_cmdline_args()
    
    if args.output_directory is None:
        print "ERROR: Output directory must to be provided. Abort."
        sys.exit(1)
    if args.image_number is None:
        print "ERROR: Image number must to be provided. Abort."
        sys.exit(1)
    
    if args.sha_version_check:
        git_sha_analysis = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=source_dir).split("\n")[-2]
        if len(subprocess.check_output(['git', 'diff'], cwd=source_dir)) > 0:
            print "There are uncommitted changes in the working tree of the loaded Git repository \'SXM-data-analysis\' (%s)." % source_dir
            print "For reproducibility reasons please commit first before submitting jobs for imagening the pre-processing routine."
            sys.exit(1)
    
    # Read image numbers from CSV table
    filename_csv = '%s/%s' % (source_dir, args.csv_file)
    print "reading filenames from: %s" % filename_csv
    with open(filename_csv, "r") as f:
        _reader = csv.DictReader(f, delimiter='\t')
        fieldNames = _reader.fieldnames
        reader = [r for r in _reader]
    all_images = [int(row["ImageNr"]) for row in reader]
    all_types = [row["ImageType"].strip() for row in reader]
    all_files = [row["FileName"].strip() for row in reader]
    all_infection_times = [int(row["InfectionTime"]) for row in reader]
    type_images = [int(row["ImageNr"]) for row in reader if row["ImageType"] in args.image_type.split(",")]
    analyzed_images = [int(row["ImageNr"]) for row in reader if int(row["VirusesTot"]) >= 0]
    if args.infection_time is not None:
        infection_time_images = [int(row["ImageNr"]) for row in reader if int(row["InfectionTime"]) == args.infection_time]
    else:
        infection_time_images = all_images
    
    # Construct list of images that shall be processed
    if args.image_number == "all":
        images = type_images
        imageTypes = all_types
        imageFiles = all_files
        infectionTimes = all_infection_times
    else:
        images_tmp = []
        tmp = args.image_number
        for s in tmp.split(','):
            if "-" in s:
                imin, imax = s.split('-')
                images_tmp += range(int(imin), int(imax)+1)
            else:
                images_tmp += [int(s)]
        # Arrays with CSV data of images to process
        images = [] # col1
        imageTypes = [] # col2
        imageFiles = [] # col3
        infectionTimes = [] # col4
        # Check whether the selected image numbers are all listed in the CSV file, if not skip
        for image in images_tmp:
            if image not in all_images:
                print "WARNING: Requested image number %i is not listed in %s and is therefore skipped for processing." % (image, filename_csv)
            elif image not in type_images:
                print "NOTE: Requested image number %i has wrong type and is therefore skipped for processing." % (image)
            elif image in analyzed_images and not args.re_process:
                print "NOTE: Requested image number %i has already been analyzed (toggle with: --re-process 1) and is therefore skipped for processing." % (image)
            elif image not in infection_time_images:
                print "NOTE: Requested image number %i does not correspond to the chosen infection time (%d h) and is therefore skipped for processing." % (image, args.infection_time)
            else:
                images.append(image)
                imageTypes.append(np.array(all_types)[np.where(np.array(all_images) == image)][0])
                imageFiles.append(np.array(all_files)[np.where(np.array(all_images) == image)][0])
                infectionTimes.append(np.array(all_infection_times)[np.where(np.array(all_images) == image)][0])
    
    # The way how the output level argument is interpreted
    save_text = args.output_level > 0 # TODO: save CSV
    save_png = args.output_level >= 2 # TODO: save pngs
    save_all = args.output_level >= 3 # TODO: save HDF5 file with processed data (histograms, etc.)
    sorters = args.sorters.split(",")
    
    # Setup common arrays for all images
    aMicroM = [] # out[0]
    allLengths = [] # out[1]
    allWidths = [] # out[2]
    allPositions = [] # out[3]
    vLengths = [] # out[4]
    vWidths = [] # out[5]
    vPositions = [] # out[6]
    vLengthsTot = [] # computed below
    vWidthsTot = [] # computed below
    vPositionsTot = [] # computed below
    counter = 0

    print("Processing %d images.." % len(imageFiles))
    if args.verbose:
        print fieldNames
    # Process single image
    for f in imageFiles:
        data = []
        for s in sorters:
            if s == "komang":
                image_path = source_dir + '/%s/cells/' % s + f + '-01.tif'
            elif s == "gunnar":
                image_path = source_dir + '/%s/cells/' % s + f + '_marked.tif'
            else:
                image_path = source_dir + '/%s/cells/' % s + f + '.tif'
            print("analyzing image: %s" % image_path)
            data.append(analyzeCrossImage(image_path, args.output_directory, save_png=save_png))
        if args.verbose > 1:
            #print data
            print len(data) # N_sorters
            print len(data[0]), len(data[1]) # N_data_types
            print len(data[0][6]), len(data[1][6]) # N_data_points
            print len(data[0][6][0]), len(data[1][6][0]) # N_data_dimensions (1D or 2D)
        vPosTemp = []
        vLenTemp = []
        vWidTemp = []
        assert (len(data) == 2)
        for i in range(len(data[0][6])):
            for j in range(len(data[1][6])):
                # Compare overlap of two virus particles from different sorters, ignore orientation of virusparticle (assume sphere with diameter vLength/2
                if (data[0][6][i][0] + data[0][4][i]/2 > data[1][6][j][0] and # 1_x + 1_dx > 2_x
                    data[0][6][i][1] + data[0][4][i]/2 > data[1][6][j][1] and # 1_y + 1_dy > 2_y
                    data[0][6][i][0] - data[0][4][i]/2 < data[1][6][j][0] and # 1_x - 1_dx < 2_x
                    data[0][6][i][1] - data[0][4][i]/2 < data[1][6][j][1] and # 1_y - 1_dy < 2_y
                    data[1][6][j][0] + data[1][4][j]/2 > data[0][6][i][0] and # 2_x + 2_dx > 1_x
                    data[1][6][j][1] + data[1][4][j]/2 > data[0][6][i][1] and # 2_y + 2_dy > 1_y
                    data[1][6][j][0] - data[1][4][j]/2 < data[0][6][i][0] and # 2_x - 2_dx < 1_x
                    data[1][6][j][1] - data[1][4][j]/2 < data[0][6][i][1]):   # 2_y - 2_dy < 1_y
                    vPosTemp.append((np.array(data[0][6][i]) + np.array(data[1][6][j]))/2.)
                    vLenTemp.append((data[0][5][i] + data[1][5][j])/2.)
                    vWidTemp.append((data[0][4][i] + data[1][4][j])/2.)
                    break # one virus particle is only allowed to overlap with one other virus particle
        if args.verbose:
            outTemp = [images[counter], imageTypes[counter], f, infectionTimes[counter]]
            for s in range(len(sorters)):
                outTemp.append(len(data[s][1]))
                outTemp.append(len(data[s][4]))
                outTemp.append(data[s][0])
            #print vPositionsTot
            #print len(vPosTemp)
            outTemp.append(-1) # TODO: ParticlesTot
            outTemp.append(len(vPosTemp))
            print outTemp
        # Add image data to common arrays
        aMicroM.append([np.array(data[s][0]) for s in range(len(sorters))])
        allLengths.append([np.array(data[s][1]) for s in range(len(sorters))])
        allWidths.append([np.array(data[s][2]) for s in range(len(sorters))])
        allPositions.append([np.array(data[s][3]) for s in range(len(sorters))])
        vLengths.append([np.array(data[s][4]) for s in range(len(sorters))])
        vWidths.append([np.array(data[s][5]) for s in range(len(sorters))])
        vPositions.append([np.array(data[s][6]) for s in range(len(sorters))])
        vLengthsTot.append(np.array(vLenTemp))
        vWidthsTot.append(np.array(vWidTemp))
        vPositionsTot.append(np.array(vPosTemp))
        counter += 1
    
    # TODO: Plot histograms
    
    # Save data
    if save_text and len(images) > 0:
        # Overwrite existing csv with updated table (reader dict)
        print "writing parameters to: %s" % filename_csv
        with open(filename_csv, "w") as f:
            _writer = csv.DictWriter(f, fieldNames, delimiter='\t')
            assert (len(sorters) == 2)
            assert (sorters[0] == 'gunnar')
            assert (sorters[1] == 'komang')
            # Make data into 2D list
            #data = [images, imageTypes, imageFiles, infectionTimes, np.array(allLengths)[:,0].shape[0], np.array(vLengths)[:,0].shape[0], np.array(aMicroM)[:,0], np.array(allLengths)[:,1].shape[0], np.array(vLengths)[:,1].shape[0], np.array(aMicroM)[:,1], [-1 for i in images], np.array(vLengthsTot).shape # BUG: shape[0] doesn't work as intended
            #data = [images, imageTypes, imageFiles, infectionTimes, [l[0].shape[0] for l in allLengths], [v[0].shape[0] for v in vLengths], [a[0] for a in aMicroM], [l[1].shape[0] for l in allLengths], [v[1].shape[0] for v in vLengths], [a[1] for a in aMicroM], [-1 for i in images], [v.shape[0] for v in vLengthsTot]] # JAS: this should work
            data = [images, imageTypes, imageFiles, infectionTimes]
            for s in range(len(sorters)): # JAS: this is prettier
                data += [[l[s].shape[0] for l in allLengths]]
                data += [[v[s].shape[0] for v in vLengths]]
                data += [[a[s].item() for a in aMicroM]] # JAS: item() gets rid of np.array(scalar), why is it needed?
            data += [[-1 for i in images]]
            data += [[v.shape[0] for v in vLengthsTot]]
            #print data
            #print aMicroM, aMicroM[0], aMicroM[0][0]
            # Update reader dict
            for i in range(len(images)):
                #print i
                row = {fieldNames[j]: data[j][i] for j in range(len(fieldNames))}
                reader[row['ImageNr']-1] = row
            _writer.writeheader()
            _writer.writerows(reader)
        print "%d rows in %s were successfully updated!" % (len(images), filename_csv)
    
