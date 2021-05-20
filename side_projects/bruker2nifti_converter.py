
import os
from bruker2nifti.converter import Bruker2Nifti

data_folder = '/Users/alex/jacques/Testing_Scans_zone/FA_catch/N57433_m00'
pfo_study_in = os.path.join(data_folder, 'raw', 'McGill_Orientation', 'a20130329_APM_DEV_Orient.j71')
pfo_study_in = data_folder
destination_folder = '/Users/alex/jacques/Testing_Scans_zone/FA_catch/'
pfo_study_out = destination_folder

bru = Bruker2Nifti(pfo_study_in, pfo_study_out, study_name='my_study')

bru.verbose = 2
bru.correct_slope = False
bru.get_acqp = False
bru.get_method = False
bru.get_reco = False
bru.nifti_version = 1
bru.qform_code = 1
bru.sform_code = 2
bru.save_human_readable = True
bru.save_b0_if_dwi = True

print(bru.scans_list)
print(bru.list_new_name_each_scan)

bru.convert()
