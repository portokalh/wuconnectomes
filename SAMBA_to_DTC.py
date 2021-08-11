from pathlib import Path
import os
import shutil

SAMBA_prep_folder = ""
DTC_DWI_folder = ""

#copytype = [shortcut, truecopy]
copytype="shortcut"

for filename in os.listdir(SAMBA_prep_folder):
    if "_coreg" in filename:
        filebasename=os.name(filename)
        if Path(filename).is_symlink():
            filename=Path('A').resolve()
        filenewpath = os.path.join(DTC_DWI_folder, filebasename)
        if copytype=="shortcut":

        elif copytype=="truecopy":
            shutil.copy(filename, )
