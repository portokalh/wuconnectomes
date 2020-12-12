from dipy.io.utils import create_tractogram_header
from tract_save import save_trk_heavy_duty
from os import path
from dipy.io.streamline import load_trk
import csv
import numpy as np
from dipy.io.streamline import load_tractogram
from dipy.tracking._utils import (_mapping_to_voxel, _to_voxel_coordinates)

from tract_manager import reducetractnumber
from tract_eval import connectivity_selection
from dipy.io.image import load_nifti, save_nifti
from tract_manager import convert_labelmask

from tract_eval import launch_quickbundles
import itk
import sys

#tract_path = "C:\\Users\\Jacques Stout\\Documents\\Work\\VBM_whiston_data\\H21593_stepsize_2_all_wholebrain_pruned.trk"
tract_path = "C:\\Users\\Jacques Stout\\Documents\\Work\\VBM_whiston_data\\H29056_stepsize_2_ratio_100_wholebrain_pruned.trk"
tract_newpath = "C:\\Users\\Jacques Stout\\Documents\\Work\\VBM_whiston_data\\H29056_stepsize_2_ratio_1000_wholebrain_pruned.trk"
labelspath = "C:\\Users\\Jacques Stout\\Documents\\Work\\VBM_whiston_data\\H29056_IITmean_RPI_labels.nii.gz"

tract_path = "/Users/alex/jacques/whiston_test/H29056_stepsize_2_ratio_100_wholebrain_pruned.trk"
tract_newpath = "/Users/alex/jacques/whiston_test/H29056_stepsize_2_ratio_1000_wholebrain_pruned.trk"
labelspath = "/Users/alex/jacques/whiston_test/H29056_IITmean_RPI_labels.nii.gz"
labels_convert_path = "/Users/alex/jacques/whiston_test/H29056_IITmean_RPI_labels_convert.nii.gz"
labels_convert_vtk = "/Users/alex/jacques/whiston_test/H29056_IITmean_RPI_labels_convert.nii.gz"

#bigtracts = load_tractogram(tract_path, 'same', bbox_valid_check=False)
if not path.exists(tract_newpath):
    trkstreamlines, affine = reducetractnumber(tract_path, tract_newpath, getdata=True, ratio=10, return_affine = True, verbose=False)
else:
    trkdata = load_trk(tract_newpath, "same")
    trkstreamlines = trkdata.streamlines
    affine = trkdata._affine

outpath = '/Users/alex/jacques/whiston_test/'

textfilepath = "C:\\Users\\Jacques Stout\\Documents\\Work\\results\\myresults_initpaired.csv"
textfilepath = "/Users/alex/jacques/whiston_test/myresults_initpaired.csv"
matrix = np.zeros((2,2))
i = 0
with open(textfilepath, newline='', encoding='utf-8-sig') as csvfile:
    matrixread = csv.reader(csvfile, delimiter=' ', quotechar='|')
    for row in matrixread:
        matrix[i, :] = row[0].split(',')
        i += 1
print(matrix)

tinystreamlines = True
if tinystreamlines:
    trkstreamlines = trkstreamlines[0:4]

if not path.exists(labels_convert_path):
    labelmask, affine_labels = load_nifti(labelspath)
    if np.size(np.shape(labelmask)) == 1:
        labelmask = labelmask[0]
    if np.size(np.shape(labelmask)) == 4:
        labelmask = labelmask[:, :, :, 0]
    print("Mask shape is " + str(np.shape(labelmask)))

    ROI_excel = 'C:\\Users\\Jacques Stout\\Documents\\Work\\VBM_whiston_data\\IITmean_RPI_index.xlsx'
    ROI_excel = "/Users/alex/jacques/whiston_test/IITmean_RPI_index.xlsx"
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
matrix[0,0] = 0
matrix[0,1] = 30
M, grouping = connectivity_selection(trkstreamlines, affine, labelmask, matrix[0,:], symmetric = True, return_mapping=True,
                                     mapping_as_streamlines=True)
