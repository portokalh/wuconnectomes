from pathlib import Path
import os
import shutil
from file_tools import buildlink
import glob
import numpy as np
import warnings

SAMBA_prep_folder = "/Volumes/Data/Badea/Lab/mouse/ADDeccode_symlink_pool/"
SAMBA_label_folder = "/Volumes/Data/Badea/Lab/mouse/VBM_21ADDecode03_IITmean_RPI_fullrun-results/connectomics/"
DTC_DWI_folder = "/Volumes/Data/Badea/ADdecode.01/Analysis/DWI/"
DTC_labels_folder = "/Volumes/Data/Badea/ADdecode.01/Analysis/DWI/"

bash_label_maker = "/Volumes/Data/Badea/Lab/mouse/create_backported_labels_for_fiber_looper.bash"

#copytype = [shortcut, truecopy]
copytype="truecopy"
overwrite=False
file_ids=["coreg", "fa", "b0", "bval", "bvec", "mask"]
subjects = ["S02654", "S02666",  "S02670",  "S02686", "S02690", "S02695", "S02720", "S02737", "S02753", "S02765", "S02781", "S02802",
            "S02804", "S02813", "S02817", "S02840", "S02877", "S02898", "S02926", "S02938", "S02939", "S02954", "S02967", "S02987",
            "S02987", "S03010", "S03017", "S03033", "S03034", "S03045", "S03048"]

for filename in os.listdir(SAMBA_prep_folder):
    if any(x in filename for x in file_ids) and any(x in filename for x in subjects):
        filepath=os.path.join(SAMBA_prep_folder,filename)
        if Path(filepath).is_symlink():
            filename=Path('A').resolve()
        filenewpath = os.path.join(DTC_DWI_folder, filename)
        if not os.path.isfile(filenewpath) or overwrite:
            print('shouldnt happen')
            if copytype=="shortcut":
                buildlink(filepath, filenewpath)
            elif copytype=="truecopy":
                shutil.copy(filepath, filenewpath)


for subject in subjects:
    subjectpath = glob.glob(os.path.join(SAMBA_label_folder, f'{subject}/'))
    if np.size(subjectpath)==1:
        subjectpath = subjectpath[0]
    elif np.size(subjectpath)>1:
        raise Exception('Too many subject folders')
    else:
        subjectpath = SAMBA_label_folder

    labelspath = glob.glob(os.path.join(subjectpath, f'{subject}*_labels.nii*'))
    if np.size(labelspath) == 1:
        labelspath = labelspath[0]
    else:
        warnings.warn("'This has not been debugged yet and will probably not work, placeholder for a better script'")
        cmd = f'bash {bash_label_maker} {subject}'
        os.system(cmd)
        labelspath = glob.glob(os.path.join(SAMBA_label_folder, f'{subject}*_labels.nii*'))
        if np.size(labelspath) == 1:
            labelspath = labelspath[0]
        else:
            raise Exception(f'Bash script {bash_label_maker} did not work, to debug')

    newlabelspath = os.path.join(DTC_labels_folder,f'{subject}_labels.nii.gz')
    if not os.path.exists(newlabelspath) or overwrite:
        shutil.copy(labelspath, newlabelspath)


