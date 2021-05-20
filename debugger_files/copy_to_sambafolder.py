import glob
from bvec_handler import checkbxh
import os
import shutil

subjects = ["01912", "02110", "02224", "02227", "02231", "02266", "02289", "02320", "02361", "02363", "02373", "02386", "02390", "02402", "02410", "02421", "02424", "02446", "02451", "02469", "02473", "02485", "02491", "02506"]
subjects = ["02110", "02224", "02227", "02231", "02266", "02289", "02320", "02361", "02363", "02373", "02386", "02390", "02402", "02410", "02421", "02424", "02446", "02451", "02469", "02473", "02485", "02491", "02506"]

inpath = "/Volumes/Data/Badea/ADdecode.01/Data/Anat/"
outpath = "/Volumes/Data/Badea/ADdecode.01/Data/SAMBA_input"

for subject in subjects:
    subjfolder = glob.glob(os.path.join(inpath, "*"+subject+"*"))[0]
    subjbxh = glob.glob(os.path.join(subjfolder, "*.bxh"))
    for bxhfile in subjbxh:
        dwidone = False
        bxhtype = checkbxh(bxhfile, False)
        if bxhtype == "dwi" and not dwidone:
            bxhname = os.path.basename(bxhfile)
            dwi_fname = bxhname.replace(".bxh", ".nii.gz")
            dwi_outputname = subject + "_dwi.nii.gz"
            bxh_outputname = subject + "_dwi.bxh"
            shutil.copyfile(os.path.join(subjfolder, bxhname), os.path.join(outpath, bxh_outputname))
            shutil.copyfile(os.path.join(subjfolder, dwi_fname), os.path.join(outpath, dwi_outputname))
            dwidone = True
        elif bxhtype =="T1":
            bxhname = os.path.basename(bxhfile)
            T1_fname = bxhname.replace(".bxh", ".nii.gz")
            T1_outputname = subject + "_T1.nii.gz"
            bxh_outputname = subject + "_T1.bxh"
            shutil.copyfile(os.path.join(subjfolder, bxhname), os.path.join(outpath, bxh_outputname))
            shutil.copyfile(os.path.join(subjfolder, T1_fname), os.path.join(outpath, T1_outputname))

