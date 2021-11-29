from transform_handler import img_transform_exec
import os, socket

subjects_list = ['N58355']

computer_name = socket.gethostname()

if 'santorini' in computer_name:
    DWI_folder = '/Users/alex/jacques/APOE_test/DWI_allsubj'
    labels_folder = '/Users/alex/jacques/APOE_test/DWI_allsubj'
    output_folder = '/Users/alex/jacques/APOE_test/DWI_allsubj_RAS/'

if 'samos' in computer_name:
    DWI_folder = '/mnt/paros_MRI/jacques/APOE/DWI_allsubj/'
    labels_folder = '/mnt/paros_MRI/jacques/APOE/DWI_allsubj/'
    output_folder = '/mnt/paros_MRI/jacques/APOE/DWI_allsubj_RAS/'

for subject in subjects_list:
    coreg_file = os.path.join(DWI_folder,f'{subject}_subjspace_coreg.nii.gz')
    labels_file = os.path.join(labels_folder,f'{subject}_labels.nii.gz')
    labels_lr_file = os.path.join(labels_folder,f'{subject}_labels_lr_ordered.nii.gz')
    coreg_RAS_file = os.path.join(output_folder,f'{subject}_coreg_RAS.nii.gz')
    labels_RAS_file = os.path.join(output_folder,f'{subject}_labels_RAS.nii.gz')
    labels_RAS_lr_file = os.path.join(output_folder,f'{subject}_labels_lr_ordered_RAS.nii.gz')
    img_transform_exec(coreg_file,'ARI','RAS',coreg_RAS_file)
    img_transform_exec(labels_file,'ARI','RAS',labels_RAS_file)
    img_transform_exec(labels_lr_file,'ARI','RAS',labels_RAS_lr_file)