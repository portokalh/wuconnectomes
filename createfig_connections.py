from dipy.io.utils import create_tractogram_header
from tract_save import save_trk_heavy_duty
from os import path
from dipy.io.streamline import load_trk
import csv
import numpy as np
from dipy.io.streamline import load_tractogram
from dipy.tracking._utils import (_mapping_to_voxel, _to_voxel_coordinates)

from tract_manager import reducetractnumber
from tract_eval import connectivity_selection, connectivity_selection_getsl
from dipy.io.image import load_nifti, save_nifti
from tract_manager import convert_labelmask
from figures_handler import connective_streamlines_figuremaker
from tract_save import save_trk_heavy_duty

import pandas as pd

def intarray(x):
  return np.int(x)

from tract_eval import launch_quickbundles
import sys

subj = "H29060"
orig_ratio = 100
if orig_ratio == 1:
    ratio = '_all_'
else:
    ratio = "_ratio_"+str(orig_ratio)+"_"
reduce_ratio = 1

#tract_path = "C:\\Users\\Jacques Stout\\Documents\\Work\\VBM_whiston_data\\H21593_stepsize_2_all_wholebrain_pruned.trk"
tract_path = "C:\\Users\\Jacques Stout\\Documents\\Work\\VBM_whiston_data\\" + subj +"_stepsize_2_ratio_100_wholebrain_pruned.trk"
tract_newpath = "C:\\Users\\Jacques Stout\\Documents\\Work\\VBM_whiston_data\\" + subj + "_stepsize_2_ratio_1000_wholebrain_pruned.trk"
labelspath = "C:\\Users\\Jacques Stout\\Documents\\Work\\VBM_whiston_data\\" + subj + "_IITmean_RPI_labels.nii.gz"
labels_convert_path = "C:\\Users\\Jacques Stout\\Documents\\Work\\VBM_whiston_data\\" + subj + "_IITmean_RPI_labels_convert.nii.gz"

tract_path = "/Volumes/dusom_dibs_ad_decode/all_staff/VBM_whiston_QA/" + subj + "_stepsize_2_ratio_100_wholebrain_pruned.trk"
tract_newpath = "//Volumes/dusom_dibs_ad_decode/all_staff/VBM_whiston_QA/" + subj + "_stepsize_2_ratio_100_wholebrain_pruned.trk"
labelspath = "/Volumes/Data/Badea/Lab/mouse/VBM_19BrainChAMD01_IITmean_RPI_with_2yr-results/connectomics/" + subj + "/" + subj + "_IITmean_RPI_labels.nii.gz"
labels_convert_path = "/Volumes/dusom_dibs_ad_decode/all_staff/VBM_whiston_labels/" + subj + "_IITmean_RPI_labels_convert.nii.gz"


tract_path = "/Volumes/dusom_dibs_ad_decode/all_staff/VBM_whiston_QA/" + subj + "_stepsize_2"+ratio+"wholebrain_pruned.trk"
tract_newpath = "//Volumes/dusom_dibs_ad_decode/all_staff/VBM_whiston_QA/" + subj + "_stepsize_2"+ratio+"wholebrain_pruned.trk"
labelspath = "/Volumes/Data/Badea/Lab/mouse/VBM_19BrainChAMD01_IITmean_RPI_with_2yr-results/connectomics/" + subj + "/" + subj + "_IITmean_RPI_labels.nii.gz"
labels_convert_path = "/Volumes/dusom_dibs_ad_decode/all_staff/VBM_whiston_labels/" + subj + "_IITmean_RPI_labels_convert.nii.gz"

#labels_convert_vtk = "/Users/alex/jacques/whiston_test/" + subj + "_IITmean_RPI_labels_convert.nii.gz"
streamlines_list_path = "/Volumes/Data/Badea/Lab/mouse/C57_JS/VBM_whiston_Figs_inclusive/" + subj + "_stepsize_2"+ratio+"wholebrain_grouping.xlsx"
anat_path = "/Volumes/Data/Badea/Lab/mouse/VBM_19BrainChAMD01_IITmean_RPI_with_2yr-results/connectomics/" + subj + "/" + subj + "_nii4D_masked_isotropic.nii.gz"

ROI_excel = 'C:\\Users\\Jacques Stout\\Documents\\Work\\VBM_whiston_data\\IITmean_RPI_index.xlsx'
ROI_excel = "/Users/alex/jacques/whiston_test/IITmean_RPI_index.xlsx"

textfilepath = "C:\\Users\\Jacques Stout\\Documents\\Work\\results\\myresults_initpaired.csv"
textfilepath = "/Users/alex/jacques/whiston_test/myresults_initpaired.csv"

outpath = '/Users/alex/jacques/whiston_test/'

save_trk = True

#bigtracts = load_tractogram(tract_path, 'same', bbox_valid_check=False)
if not path.exists(tract_newpath):
    trkstreamlines, affine = reducetractnumber(tract_path, tract_newpath, getdata=True, ratio=reduce_ratio, return_affine=True,
                                               verbose=False)
else:
    trkdata = load_trk(tract_newpath, "same")
    trkstreamlines = trkdata.streamlines
    affine = trkdata._affine
    if hasattr(trkdata, 'space_attribute'):
        header = trkdata.space_attribute
    elif hasattr(trkdata, 'space_attributes'):
        header = trkdata.space_attributes


matrix = np.zeros((100,2))
i = 0
with open(textfilepath, newline='', encoding='utf-8-sig') as csvfile:
    matrixread = csv.reader(csvfile, delimiter=' ', quotechar='|')
    for row in matrixread:
        matrix[i, :] = row[0].split(',')
        i += 1
print(matrix)

#tinystreamlines = True
#if tinystreamlines:
#    trkstreamlines = trkstreamlines[0:4]

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
"""
    imageFileName = labels_convert_path

    Dimension = 2
    PixelType = itk.UC
    ImageType = itk.Image[PixelType, Dimension]

    reader = itk.ImageFileReader[ImageType].New()
    reader.SetFileName(imageFileName)

    itkToVtkFilter = itk.ImageToVTKImageFilter[ImageType].New()
    itkToVtkFilter.SetInput(reader.GetOutput())

    itkToVtkFilter.Update()
    myvtkImageData = itkToVtkFilter.GetOutput()
    print(myvtkImageData)
"""
#launch_quickbundles(trkstreamlines, outpath, ROIname="all", labelmask=labelmask, affine=affine_labels, interactive=True)
##matrix[0,0] = 0
#matrix[0,1] = 30

grouping = pd.read_excel(streamlines_list_path)
group_matrix = grouping.values
select_streamgroups = 2
intarray2 = np.vectorize(intarray)
ROI_names = []
ROI_streamlines = []
for i in np.arange(select_streamgroups):
    ROI_tuple = intarray2(matrix[i,:])
    #ROI_tuple = np.array([4,10])
    select_streamlines = group_matrix[ROI_tuple[0]-1, ROI_tuple[1]]
    select_new = []
    select_streamlines = select_streamlines.replace('[', '')
    select_streamlines = select_streamlines.replace(']','')
    select_streamlines = (select_streamlines.split(','))
    select_streamlines = [int(i) for i in select_streamlines]
    ROI_streamlines.append(select_streamlines)
    ROI_names.append([group_matrix[ROI_tuple[0]-1,0], group_matrix[ROI_tuple[1]-1,0]])

if save_trk:
    for i in np.arange(select_streamgroups):
        pathfile_name = outpath + ROI_names[i][0][0:10] + ROI_names[i][1][0:10] + '.trk'
        save_trk_heavy_duty(pathfile_name, streamlines = Streamlines(ROI_streamlines), affine = affine, header = myheader)

connective_streamlines_figuremaker(trkstreamlines, ROI_streamlines, ROI_names, anat_path)
#M, grouping = connectivity_selection(trkstreamlines, affine, labelmask, matrix[0,:], symmetric = True, return_mapping=True,
#                                     mapping_as_streamlines=True)
#M, grouping = connectivity_selection_getsl(trkstreamlines, affine, labelmask, symmetric = True, return_mapping=True,
#                                     mapping_as_streamlines=True)

