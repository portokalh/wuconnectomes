
import os
import sys
import numpy as np
import re

if np.size(sys.argv)>1:
    attribute_file = sys.argv[1]
else:
    attribute_file = "./mainlaunch.rtf"

if not os.path.exists(attribute_file):
    print("The file could not be found at "+attribute_file+", end of scrypt run")

verbose=True
with open(attribute_file, 'rb') as source:
    if verbose: print('INFO    : Extracting pipeline parameters')
    #header_size = source.read(4)
    #header_size = struct.unpack('I', header_size)
    #if verbose: print('INFO    : Header size = ', int(header_size[0]))
    i = 0
    stopsign = 200
    pattern1 = 'Diffusion folder'
    rx1 = re.compile(pattern1, re.IGNORECASE | re.MULTILINE | re.DOTALL)
    pattern2 = 'Tract folder'
    rx2 = re.compile(pattern2, re.IGNORECASE | re.MULTILINE | re.DOTALL)
    pattern3 = 'Figures folder'
    rx3 = re.compile(pattern3, re.IGNORECASE | re.MULTILINE | re.DOTALL)
    pattern4 = 'Atlas file'
    rx4 = re.compile(pattern4, re.IGNORECASE | re.MULTILINE | re.DOTALL)
    pattern5 = 'Subject list'
    rx5 = re.compile(pattern5, re.IGNORECASE | re.MULTILINE | re.DOTALL)
    pattern6 = 'Step size'
    rx6 = re.compile(pattern5, re.IGNORECASE | re.MULTILINE | re.DOTALL)
    pattern7 = 'ROI targets'
    rx7 = re.compile(pattern5, re.IGNORECASE | re.MULTILINE | re.DOTALL)
    pattern8 = 'Denoise'
    rx8 = re.compile(pattern5, re.IGNORECASE | re.MULTILINE | re.DOTALL)
    pattern9 = 'Ratio of streams'
    rx9 = re.compile(pattern5, re.IGNORECASE | re.MULTILINE | re.DOTALL)
    pattern10 = 'Overwrite'
    rx10 = re.compile(pattern5, re.IGNORECASE | re.MULTILINE | re.DOTALL)
    patternmax = 'Max Processors'
    rxmax = re.compile(pattern5, re.IGNORECASE | re.MULTILINE | re.DOTALL)
    patternverb = 'Verbose'
    rxverb = re.compile(pattern5, re.IGNORECASE | re.MULTILINE | re.DOTALL)

    for line in source:
        i += 1
        if i == stopsign:
            print("hi")
        for a in rx1.findall(str(line)):
            diff_fol = str(line).split(':')[1]
            diff_fol = diff_fol.split('\\')[0]
        for a in rx2.findall(str(line)):
            trk_fol = str(line).split(':')[1]
            trk_fol = trk_fol.split('\\')[0]
        for a in rx3.findall(str(line)):
            fig_fol = str(line).split(':')[1]
            fig_fol = fig_fol.split('\\')[0]
        for a in rx4.findall(str(line)):
            atlas_file = str(line).split(':')[1]
            atlas_file = atlas_file.split('\\')[0]
        for a in rx5.findall(str(line)):
            subj_list = str(line).split(':')[1]
            subj_list = subj_list.split('\\')[0]
        for a in rx6.findall(str(line)):
            stepsize = str(line).split(':')[1]
            stepsize = float(stepsize.split('\\')[0])
        for a in rx7.findall(str(line)):
            targetrois = str(line).split(':')[1]
            targetrois = targetrois.split('\\')[0]
        for a in rx8.findall(str(line)):
            denoise = str(line).split(':')[1]
            denoise = denoise.split('\\')[0]
        for a in rx9.findall(str(line)):
            ratio = str(line).split(':')[1]
            ratio = int(ratio.split('\\')[0])
        for a in rx10.findall(str(line)):
            overwrite = str(line).split(':')[1]
            overwrite = overwrite.split('\\')[0]
        for a in rxverb.findall(str(line)):
            verbose = str(line).split(':')[1]
            verbose = verbose.split('\\')[0]
        for a in rxmax.findall(str(line)):
            max_processors = str(line).split(':')[1]
            max_processors = int(max_processors.split('\\')[0])

if max_processors not in locals():
    max_processors = 1
if verbose not in locals():
    verbose = True
if overwrite not in locals():
    overwrite = False
if ratio not in locals():
    ratio = 1
if denoise not in locals():
    ratio = "None"
if targetrois not in locals():
    targetrois = "Whole Brain"
if stepsize not in locals():
    stepsize = 2
if subj_list or atlas_file or fig_fol or trk_fol or diff_fol not in locals():
    txt = "Could not find variable of "
    if diff_fol not in locals():
        txt += pattern1 + ", "
    if trk_fol not in locals():
        txt += pattern2 + ", "
    if fig_fol not in locals():
        txt += pattern3 + ", "
    if atlas_file not in locals():
        txt += pattern4 + ", "
    if subj_list not in locals():
        txt += pattern5 + ", "
    txt += " in the txt folder located at " + (os.path.abspath(attribute_file))

subject_processes = np.size(subj)
if max_processors < subject_processes:
    subject_processes = max_processors
function_processes = np.int(max_processors/subject_processes)

dwi_results = []
tract_results = []



trkroi = targetrois
if len(trkroi)==1:
    roistring = "_" + trkroi[0] + "_"
elif len(trkroi)>1:
    roistring="_"
    for roi in trkroi:
        roistring = roistring + roi[0:4]
    roistring = roistring + "_"

if ratio == 1:
    str_identifier = "ratio_" + str(ratio) + '_stepsize_' + str(stepsize)
else:
    str_identifier = "ratio_" + str(ratio) + '_stepsize_' + str(stepsize)

if len(targetrois)==1:
    targetroistring = "_" + targetrois[0] + "_"
elif len(targetrois)>1:
    roistring="_"
    for roi in rois:
        targetoistring = roistring + roi[0:4]
    targetroistring = roistring + "_"

labelslist=[]
for roi in targetrois:
    rslt_df = df.loc[df['Structure'] == roi.lower()]
    if roi.lower() == "wholebrain" or roi.lower() == "brain":
        labelslist=None
    else:
        labelslist=np.concatenate((labelslist,np.array(rslt_df.index2)))
print(labelslist)
if isempty(labelslist) and roi.lower() != "wholebrain" and roi.lower() != "brain":
    txt = "Warning: Unrecognized roi, will take whole brain as ROI. The roi specified was: " + roi
    print(txt)

bvec_orient = [-2, 1, 3]

savefa = True

if subject_processes>1:
    if function_processes>1:
        pool = MyPool(subject_processes)
    else:
        pool = mp.Pool(subject_processes)

    tract_results = pool.starmap_async(create_tracts, [(diff_fol, trk_fol, subject, stepsize, function_processes, str_identifier,
                                                            ratio, savefa, labelslist, bvec_orient, verbose) for subject in
                                                           l]).get()
#    tract_results = pool.starmap_async(evaluate_tracts, [(dwipath, outtrkpath, subject, stepsize, saved_streamlines,
#                                                          labelslist, outpathpickle, figspath, function_processes,
#                                                          doprune, display, verbose) for subject in l]).get()
    pool.close()
else:
    for subject in subj_list:
        dwi_results.append(dwi_preprocessing(dwipath, outpath, subject, bvec_orient, denoise, savefa, function_processes,
                          labelslist, str_identifier, verbose=False)
        tract_results.append(create_tracts(dwipath, outtrkpath, subject, stepsize, function_processes, strproperty,
                                         saved_streamlines, savefa, labelslist, bvec_orient, verbose))
