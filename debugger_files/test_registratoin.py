from os.path import expanduser, join
import numpy as np
from dipy.viz import regtools
#from dipy.data import fetch_stanford_hardi
from dipy.data.fetcher import fetch_syn_data
from dipy.io.image import load_nifti
from dipy.align.imaffine import (transform_centers_of_mass,
                                 AffineMap,
                                 MutualInformationMetric,
                                 AffineRegistration)
from dipy.align.transforms import (TranslationTransform3D,
                                   RigidTransform3D,
                                   AffineTransform3D)

home = 'D:\\Personal\\School\\Duke\\2020-2021\\Extracurriculurs\\Bass Connections'
home = "/Users/alex/jacques/Nikhiltest/"
dname = join(home, 'Data', 'test_images')
dname = "/Users/alex/jacques/Nikhiltest/"

#Specifying data file locations
fdwi = join(dname, 'N57437_nii4D.nii')

fbval = join(dname, 'N57437_bvals.txt')

fbvec = join(dname, 'N57437_bvec.txt')

#Specifying second image (moving)
bdwi = join(dname, 'N57442_nii4D.nii')

bbval = join(dname, 'N57442_bvals.txt')

bbvec = join(dname, 'N57442_bvec.txt')
fdwi = "/Volumes/Data/Badea/Lab/mouse/C57_JS/Nikhil_temp/N57437_nii4D_RAS.nii"
fdwi = "/Volumes/Data/Badea/Lab/mouse/C57_JS/Nikhil_temp/niftifiles/N57437_nii4D_RAS.nii"
fdwi = "/Volumes/Data/Badea/Lab/mouse/C57_JS/Nikhil_temp2/niftifiles/N57437_nii4D_RAS.nii"
#Loading static image
static_data, static_affine, static_img = load_nifti(fdwi, return_img=True)
#static_data, static_affine, static_img = load_nifti(fdwi, return_img=True)

static = np.squeeze(static_data)[..., 0]
static_grid2world = static_affine


#Loading moving image
moving_data, moving_affine, moving_img = load_nifti(bdwi, return_img=True)

moving = moving_data
moving_grid2world = moving_affine


#Visualize images
identity = np.eye(4)
affine_map = AffineMap(identity,
                       static.shape, static_grid2world,
                       moving.shape, moving_grid2world)
resampled = affine_map.transform(moving)

regtools.overlay_slices(static, resampled, None, 0,
                        "Static", "Moving", "resampled_0.png")
regtools.overlay_slices(static, resampled, None, 1,
                        "Static", "Moving", "resampled_1.png")
regtools.overlay_slices(static, resampled, None, 2,
                        "Static", "Moving", "resampled_2.png")
