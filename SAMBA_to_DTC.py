from pathlib import Path
import os
import shutil
from file_tools import buildlink

SAMBA_prep_folder = "/Volumes/Data/Badea/Lab/mouse/ADDeccode_symlink_pool/"
DTC_DWI_folder = "/Volumes/Data/Badea/ADdecode.01/Analysis/DWI/"

#copytype = [shortcut, truecopy]
copytype="truecopy"
overwrite=False
file_ids=["coreg", "fa", "b0", "bval", "bvec", "mask"]
subjects = ["02654", "02666",  "02670",  "02686", "02690", "02695", "02720", "02737", "02753", "02765", "02781", "02802",
            "02804", "02813", "02817", "02840", "02877", "02898", "02926", "02938", "02939", "02954", "02967", "02987",
            "02987", "03010", "03017", "03033", "03034", "03045", "03048"]

#if you want everybody, just add S or 0 or something that is common to all subject id
for filename in os.listdir(SAMBA_prep_folder):
    if any(x in filename for x in file_ids) and any(x in filename for x in subjects):
        filepath=os.path.join(SAMBA_prep_folder,filename)
        if Path(filepath).is_symlink():
            filename=Path('A').resolve()
        filenewpath = os.path.join(DTC_DWI_folder, filename)
        if not os.path.isfile(filenewpath) or overwrite:
            if copytype=="shortcut":
                buildlink(filepath, filenewpath)
            elif copytype=="truecopy":
                shutil.copy(filepath, filenewpath)
