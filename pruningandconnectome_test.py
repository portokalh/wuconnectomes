
from tract_handler import prune_streamlines
from dipy.io.utils import (create_tractogram_header)
import numpy as np
from dipy.io.streamline import load_trk
from dipy.tracking.streamline import Streamlines
from dipy.io.image import load_nifti
import tract_save
from dipy.tracking import utils
from tract_manager import connectomes_to_excel

labels_nii = load_nifti("/Volumes/Data/Badea/Lab/mouse/VBM_19BrainChAMD01_IITmean_RPI_with_2yr-results/connectomics/H21593/H21593_IITmean_RPI_labels.nii.gz")
anat_nii = load_nifti("/Volumes/Data/Badea/Lab/mouse/VBM_19BrainChAMD01_IITmean_RPI_with_2yr-results/connectomics/H21593/H21593_nii4D_masked_isotropic.nii.gz")
trkfile = "/Volumes/Data/Badea/Lab/mouse/C57_JS/VBM_whiston_QA/H21593_stepsize_2_all_wholebrain.trk"
#trkfile = "/Volumes/Data/Badea/Lab/mouse/C57_JS/Testzone/H21593_wholebrain_ratio_100_pruned.trk"
trkprunepath = "/Volumes/Data/Badea/Lab/mouse/C57_JS/Testzone/H21593_stepsize_2_all_wholebrain_pruned.trk"
ROI_excel = '/Volumes/Data/Badea/Lab/atlases/IITmean_RPI/IITmean_RPI_lookup.xlsx'

verbose = True

affine_streams = np.eye(4)

trkdata = load_trk(trkfile, "same")
affine = trkdata._affine
trkdata.to_vox()
trkstreamlines = trkdata.streamlines
#trkprunepath = os.path.dirname(trkfilepath) + '/' + subject + str_identifier + '_pruned.trk'
cutoff = 4
pruned_streamlines = prune_streamlines(list(trkstreamlines), anat_nii[0][:,:,:,0], cutoff=cutoff, verbose=verbose)
pruned_streamlines_SL = Streamlines(pruned_streamlines)
if hasattr(trkdata, 'space_attribute'):
    header = trkdata.space_attribute
elif hasattr(trkdata, 'space_attributes'):
    header = trkdata.space_attributes
myheader = create_tractogram_header(trkprunepath, *header)
prune_sl = lambda: (s for s in pruned_streamlines)
tract_save.save_trk_heavy_duty(trkprunepath, streamlines=prune_sl, affine=affine, header=myheader)
del (prune_sl, pruned_streamlines, trkdata)

affine_streams = np.eye(4)

M, grouping = utils.connectivity_matrix(pruned_streamlines_SL, affine_streams, labels_nii[0],
                                        return_mapping=True,
                                        mapping_as_streamlines=True)

excel_path = "/Volumes/Data/Badea/Lab/mouse/C57_JS/Testzone/H21593_connectomes_100test.xlsx"
connectomes_to_excel(M, ROI_excel, excel_path)
