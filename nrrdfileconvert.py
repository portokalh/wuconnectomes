import os
from glob import glob
import nrrd #pip install pynrrd, if pynrrd is not already installed
import nibabel as nib #pip install nibabel, if nibabel is not already installed
import numpy as np

"""
baseDir = os.path.normpath('path/to/file/')
files = glob(baseDir+'/*.nrrd')

for file in files:
#load nrrd
  _nrrd = nrrd.read(file)
  data = _nrrd[0]
  header = _nrrd[1]
"""

file = "/Users/alex/Downloads/average_template_10.nrrd"
mynrrd = nrrd.read(file)
data = mynrrd[0]
header = mynrrd[1]
img = nib.Nifti1Image(data, np.eye(4))
nib.save(img,"/Users/alex/jacques/atlastest/average_template_10.nii.gz")