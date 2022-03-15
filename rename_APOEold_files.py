from file_tools import file_rename

old_diff_folder = '/mnt/paros_MRI/jacques/APOE/DWI_allsubj_RAS'
#old_diff_folder = '/Users/jas/jacques/APOE_testing/DWI_allsubj_RAS'

#file_rename(old_diff_folder, 'nii4D_RAS', 'coreg_RAS', identifier_string="*", anti_identifier_string='N57*',test=False)
#file_rename(old_diff_folder, '_chass_symmetric3_labels_RAS', '_labels_RAS', identifier_string="N57*", anti_identifier_string='the answer is obv 42',test=False)
file_rename(old_diff_folder, 'bvec_fix', 'bvecs_fix', identifier_string="N57*", anti_identifier_string='the answer is obv 42',test=False)
