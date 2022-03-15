from transform_handler import img_transform_exec
from file_tools import mkcdir
import os, socket, glob

subjects_list = ["N58214", "N58215",
     "N58216", "N58217", "N58218", "N58219", "N58221", "N58222", "N58223", "N58224",
                "N58225", "N58226", "N58228",
                "N58229", "N58230", "N58231", "N58232", "N58633", "N58634", "N58635", "N58636", "N58649", "N58650",
                "N58651", "N58653", "N58654",
                'N58408', 'N58714', 'N58740', 'N58477', 'N58734', 'N58309', 'N58792', 'N58302',
                'N58784', 'N58706', 'N58361', 'N58355', 'N58712', 'N58790', 'N58606', 'N58350', 'N58608',
                'N58779', 'N58500', 'N58604', 'N58749', 'N58510', 'N58394', 'N58346', 'N58344', 'N58788', 'N58305',
                'N58514', 'N58794', 'N58733', 'N58655', 'N58735', 'N58310', 'N58400', 'N58708', 'N58780', 'N58512',
                'N58747', 'N58303', 'N58404', 'N58751', 'N58611', 'N58745', 'N58406', 'N58359', 'N58742', 'N58396',
                'N58613', 'N58732', 'N58516', 'N58402']
subjects_list = ["N58831","N59022","N59026","N59033","N59035","N59039","N59041","N59065","N59066","N59072","N59076","N59078","N59080","N59097","N59099","N59109","N59116","N59118","N59120"]
#removed 'N58398' could not find
computer_name = socket.gethostname()

if 'santorini' in computer_name:
    DWI_folder = '/Users/alex/jacques/APOE_test/DWI_allsubj'
    labels_folder = '/Users/alex/jacques/APOE_test/DWI_allsubj'
    output_folder = '/Users/alex/jacques/APOE_test/DWI_allsubj_RAS/'

if 'samos' in computer_name:
    DWI_folder = '/mnt/paros_MRI/jacques/APOE/DWI_allsubj/'
    labels_folder = '/mnt/paros_MRI/jacques/APOE/DWI_allsubj/'
    output_folder = '/mnt/paros_MRI/jacques/APOE/DWI_allsubj_RAS/'

subjects_all = glob.glob(os.path.join(DWI_folder,'*_subjspace_coreg.nii.gz'))
subjects_list = []
for subject in subjects_all:
    subject_name = os.path.basename(subject)
    subjects_list.append(subject_name[:6])
subjects_list.sort()
subjects_list = subjects_list[:]
print(subjects_list)
mkcdir(output_folder)


for subject in subjects_list:
    coreg_file = os.path.join(DWI_folder,f'{subject}_subjspace_coreg.nii.gz')
    labels_file = os.path.join(labels_folder,f'{subject}_labels.nii.gz')
    labels_lr_file = os.path.join(labels_folder,f'{subject}_labels_lr_ordered.nii.gz')
    mask_file = os.path.join(DWI_folder, f'{subject}_subjspace_mask.nii.gz')
    coreg_RAS_file = os.path.join(output_folder,f'{subject}_coreg_RAS.nii.gz')
    labels_RAS_file = os.path.join(output_folder,f'{subject}_labels_RAS.nii.gz')
    labels_RAS_lr_file = os.path.join(output_folder,f'{subject}_labels_lr_ordered_RAS.nii.gz')
    mask_RAS_file = os.path.join(output_folder, f'{subject}_RAS_mask.nii.gz')
    transferred=0
    if not os.path.exists(coreg_RAS_file):
        img_transform_exec(coreg_file,'ARI','RAS',coreg_RAS_file)
        transferred=1
    if not os.path.exists(labels_RAS_file):
        img_transform_exec(labels_file,'ARI','RAS',labels_RAS_file)
        transferred=1
    if not os.path.exists(labels_RAS_lr_file): 
        if os.path.exists(labels_lr_file):
            img_transform_exec(labels_lr_file,'ARI','RAS',labels_RAS_lr_file)
            transferred=1
    if not os.path.exists(mask_RAS_file):
        img_transform_exec(mask_file,'ARI','RAS',mask_RAS_file)
        transferred=1
    if transferred:
        print(f'transferred subject {subject}')
    else:
        print(f'already transferred subject {subject}')
