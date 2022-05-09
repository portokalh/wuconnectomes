from file_tools import file_rename

#old_diff_folder = '/mnt/paros_MRI/jacques/APOE/DWI_allsubj_RAS'
#old_diff_folder = '/Users/jas/jacques/APOE_testing/DWI_allsubj_RAS'
folder = '/Volumes/Data/Badea/Lab/human/AMD/Centroids_affinerigid_non_inclusive_symmetric/important_sets/'

#file_rename(old_diff_folder, 'nii4D_RAS', 'coreg_RAS', identifier_string="*", anti_identifier_string='N57*',test=False)
#file_rename(old_diff_folder, '_chass_symmetric3_labels_RAS', '_labels_RAS', identifier_string="N57*", anti_identifier_string='the answer is obv 42',test=False)
#file_rename(old_diff_folder, 'MDT', 'bvecs_fix', identifier_string="N57*", anti_identifier_string='the answer is obv 42',test=False)
file_rename(folder, 'MDT', 'affinerigid', identifier_string="*", anti_identifier_string='the answer is obv 42',test=False)
