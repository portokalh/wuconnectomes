
"""
Created by Jacques
Get full set of parameters given a set of trks
Part of project for TN-PCA using Whe
"""

from tract_manager import get_tract_params
import xlsxwriter
import numpy as np
import multiprocessing as mp

excel_path = "/Users/alex/jacques/connectomes_testing/tracts_summary.xlsx"
workbook = xlsxwriter.Workbook(excel_path)
worksheet = workbook.add_worksheet()

worksheet.write(0, 1, "Number of Tracts")
worksheet.write(0, 2, "Minimum Tract length")
worksheet.write(0, 3, "Maximum Tract length")
worksheet.write(0, 4, "Mean Tract length")
worksheet.write(0, 5, "Standard Deviancy")

trkpath = "/Volumes/Data/Badea/Lab/mouse/C57_JS/TRK_RAS_40subj/N57442_small_fimbria_stepsize_2_pruned.trk"
trkpath = "/Volumes/Data/Badea/Lab/mouse/C57_JS/TRK_RAS_40subj/"

l = ['N57437', 'N57442', 'N57446', 'N57447','N57449','N57451','N57496','N57498','N57500','N57502','N57504','N57513','N57515','N57518', 'N57520','N57522','N57546','N57447','N57449','N57451','N57496','N57498','N57500','N57502','N57504','N57513','N57515','N57518','N57520','N57522','N57546','N57548', 'N57550', 'N57552', 'N57554', 'N57559', 'N57580', 'N57582', 'N57584', 'N57587', 'N57590', 'N57692', 'N57694', 'N57700', 'N57702', 'N57709']
l = ['N57700']

max_processors = 1
subject_processes = np.size(l)
if max_processors < subject_processes:
    subject_processes = max_processors
function_processes = np.int(max_processors/subject_processes)

str_identifier = "_wholebrain_all_stepsize_2_pruned"
verbose = True


tract_params = []

if subject_processes>1:
    if function_processes>1:
        pool = MyPool(subject_processes)
    else:
        pool = mp.Pool(subject_processes)

    tract_params = pool.starmap_async(get_tract_params, [(trkpath, subject, str_identifier, verbose) for subject in
                                                           l]).get()
#    tract_results = pool.starmap_async(evaluate_tracts, [(dwipath, outtrkpath, subject, stepsize, saved_streamlines,
#                                                          labelslist, outpathpickle, figspath, function_processes,
#                                                          doprune, display, verbose) for subject in l]).get()
    pool.close()
else:
    for subject in l:
        tract_params.append(get_tract_params(trkpath, subject, str_identifier, verbose))

for i in np.arange(1,np.shape(tract_params)[0]+1):
    for j in np.arange(0,6):
        worksheet.write(i, j, tract_params[i-1][j])

workbook.close()