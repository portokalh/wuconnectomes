from os import path
from dipy.io.streamline import load_trk
import csv
import numpy as np
from dipy.io.image import load_nifti, save_nifti
from tract_manager import convert_labelmask
from tract_save import save_trk_heavy_duty
from dipy.io.utils import (create_tractogram_header)
from dif_to_trk import check_for_fa
from tract_handler import get_connectome_attributes
import pandas as pd

def intarray(x):
  return np.int(x)

from tract_eval import launch_quickbundles
import sys

#2yr
subj = "H29060"
#subj = "H29410"

#init
#subj = "H26966"
#subj = "H26637"

orig_ratio = 10
if orig_ratio == 1:
    oldratio = '_all_'
else:
    oldratio = "_ratio_"+str(orig_ratio)+"_"
reduce_ratio = 10
if reduce_ratio == 1:
    newratio = '_all_'
else:
    newratio = "_ratio_"+str(reduce_ratio)+"_"

#tract_path = "C:\\Users\\Jacques Stout\\Documents\\Work\\VBM_whiston_data\\H21593_stepsize_2_all_wholebrain_pruned.trk"

tract_dir = "/Volumes/Data/Badea/Lab/mouse/C57_JS/VBM_whiston_QA_new/"
#tract_dir = "/Volumes/dusom_dibs_ad_decode/all_staff/VBM_whiston_QA/"
diff_dir = "/Volumes/Data/Badea/Lab/mouse/VBM_19BrainChAMD01_IITmean_RPI_with_2yr-results/connectomics"
labels_output_dir = "/Volumes/dusom_dibs_ad_decode/all_staff/VBM_whiston_labels/"
connectomes_dir = "/Volumes/Data/Badea/Lab/mouse/C57_JS/VBM_whiston_Figs_inclusive_new/"
ROI_excel = "/Users/alex/jacques/whiston_test/IITmean_RPI_index.xlsx"
textfilepath = "/Users/alex/jacques/whiston_test/myresults_initpaired.csv"


textfilepath = "/Users/alex/jacques/whiston_test/myresults_diff.csv"
textfilepath = None

BIGGUS_DISKUS = "/Volumes/Data/Badea/Lab/mouse"
fapath = BIGGUS_DISKUS + "/C57_JS/diff_whiston_preprocessed/"
~, fa = check_for_fa(fapath, subj, getdata=True)

verbose = True

stepsize = 2
masktype = "FA"
masktype = "binary"
trkroi = ["wholebrain"]
if len(trkroi)==1:
    roistring = "_" + trkroi[0] #+ "_"
elif len(trkroi)>1:
    roistring="_"
    for roi in trkroi:
        roistring = roistring + roi[0:4]
    roistring = roistring #+ "_"
if reduce_ratio == 1:
    saved_ratiostreamlines = "_all"
else:
    saved_streamlines = "_ratio_" + str(reduce_ratio)


stepsize = 2
trkroi = ["wholebrain"]
if len(trkroi)==1:
    roistring = "_" + trkroi[0] #+ "_"
elif len(trkroi)>1:
    roistring="_"
    for roi in trkroi:
        roistring = roistring + roi[0:4]
    roistring = roistring #+ "_"
if reduce_ratio == 1:
    saved_ratiostreamlines = "_all"
else:
    saved_streamlines = "_ratio_" + str(reduce_ratio)

if masktype == "FA":
    maskuse = "_fa"
else:
    maskuse = "_binary"

str_identifier = '_stepsize_' + str(stepsize) + maskuse + roistring + saved_streamlines

streamlines_grouping_path = connectomes_dir + subj + str_identifier + "_grouping.xlsx"
anat_path = diff_dir + subj + "/" + subj + "_nii4D_masked_isotropic.nii.gz"
anat_path = path.join(diff_dir, subj, subj + "_nii4D_masked_isotropic.nii.gz")

outpath = "/Volumes/Data/Badea/Lab/mouse/C57_JS/Whiston_figures_files/" + subj + "*/"
import glob
dir = glob.glob(outpath)
outpath = dir[0]
save_trk = True

tract_path = tract_dir + subj + str_identifier + "_pruned.trk"
labelspath = path.join(diff_dir, subj, subj + "_IITmean_RPI_labels.nii.gz")
labels_convert_path = labels_output_dir + subj + "_IITmean_RPI_labels_convert.nii.gz"


### Change this to a list with append that is converted to a matrix instead of this
if textfilepath is not None:
    roi_selection_matrix = np.zeros((100,2))
    i = 0
    with open(textfilepath, newline='', encoding='utf-8-sig') as csvfile:
        matrixread = csv.reader(csvfile, delimiter=' ', quotechar='|')
        for row in matrixread:
            roi_selection_matrix[i, :] = row[0].split(',')
            i += 1
    print(roi_selection_matrix)

if not path.exists(labels_convert_path):
    labelmask, affine_labels = load_nifti(labelspath)
    if np.size(np.shape(labelmask)) == 1:
        labelmask = labelmask[0]
    if np.size(np.shape(labelmask)) == 4:
        labelmask = labelmask[:, :, :, 0]
    print("Mask shape is " + str(np.shape(labelmask)))


    labelmask_new = convert_labelmask(ROI_excel, labelmask)

    save_nifti(labels_convert_path, labelmask_new, affine_labels)
    labelmask_new = labelmask

else:
    labelmask, affine_labels = load_nifti(labels_convert_path)




grouping = pd.read_excel(streamlines_grouping_path)
group_matrix = grouping.values
select_streamgroups = 1
intarray2 = np.vectorize(intarray)
ROI_names = []
ROI_streamlines_all = []

if textfilepath is None:

    trkdata = load_trk(tract_path, "same")
    trkdata.to_vox()
    trkstreamlines = trkdata.streamlines
    affine = trkdata._affine
    if hasattr(trkdata, 'space_attribute'):
        header = trkdata.space_attribute
    elif hasattr(trkdata, 'space_attributes'):
        header = trkdata.space_attributes


    for i in np.arange(np.size(group_matrix),1)
        for j in np.arange(np.size(group_matrix),1)
            ROI_tuple = intarray2(i,j)
            select_streamlines = group_matrix[ROI_tuple[0] - 1, ROI_tuple[1]]
            select_new = []
            select_streamlines = select_streamlines.replace('[', '')
            select_streamlines = select_streamlines.replace(']', '')
            select_streamlines = (select_streamlines.split(','))
            try:
                select_streamlines.remove("")
            except ValueError:
                pass
            try:
                select_streamlines.remove(" ")
            except ValueError:
                pass

            ROI_streamlines = [int(i) for i in select_streamlines]
            ROI_streamlines_all.append(ROI_streamlines)
            ROI_name = ([group_matrix[ROI_tuple[0] - 1, 0], group_matrix[ROI_tuple[1] - 1, 0]])
            ROI_names.append([group_matrix[ROI_tuple[0] - 1, 0], group_matrix[ROI_tuple[1] - 1, 0]])

            trk_ROI_streamlines = []
            for ROI_stream in ROI_streamlines[i]:
                trk_ROI_streamlines.append(trkstreamlines[ROI_stream])
            # trk_ROI_streamlines = trkstreamlines[ROI_streamlines]
            pathfile_name = outpath + subj + str_identifier + ROI_name[0][0:11] + "_" + ROI_name[1][0:11] + '.trk'
            ROI_sl = lambda: (s for s in trk_ROI_streamlines)
            myheader = create_tractogram_header(pathfile_name, *header)
            save_trk_heavy_duty(pathfile_name, streamlines=ROI_sl, affine=affine, header=myheader)

            affine_streams = np.eye(4)
            favals, mdvals, numtracts, minlength, maxlength, meanlength, stdlength = get_connectome_attributes(ROI_sl, fa, md=None, verbose)
