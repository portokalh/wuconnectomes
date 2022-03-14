import os, glob
from nifti_handler import fix_bvals_bvecs

path = '/mnt/paros_MRI/jacques/APOE/DWI_allsubj_RAS'
#path = '/Users/jas/jacques/APOE_testing/bvalues_test'
subjects_files = glob.glob(os.path.join(path,'*bvals.txt'))
subjects = []
for subject_file in subjects_files:
    subjects.append(os.path.basename(subject_file)[:6])

for subject in subjects:
    fbvals = os.path.join(path, subject+'_bvals.txt')
    fbvecs = os.path.join(path, subject+'_bvecs.txt')
    fbvals_new = os.path.join(path, subject+'_bvals_fix.txt')
    fbvecs_new = os.path.join(path, subject+'_bvecs_fix.txt')
    if not os.path.exists(fbvals_new) or not os.path.exists(fbvecs_new):
        fix_bvals_bvecs(fbvals,fbvecs)