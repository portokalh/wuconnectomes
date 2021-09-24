from pathlib import Path
import os
import shutil
from file_tools import buildlink
import glob
import numpy as np
import warnings
from create_backported_labels import create_backport_labels
import os
import getpass

#project = ["AD_Decode", "APOE"]
project = "APOE"
verbose = True
mainpath = "/Volumes/Data/Badea/Lab/"
#mainpath = "/mnt/munin6/Badea/Lab/"
SAMBA_headfile_dir = os.path.join(mainpath, "samba_startup_cache")

if project == "AD_Decode":

    SAMBA_mainpath = os.path.join(mainpath, "mouse")
    SAMBA_projectname = "VBM_21ADDecode03_IITmean_RPI_fullrun"
    SAMBA_headfile = os.path.join(SAMBA_headfile_dir, "jas297_SAMBA_ADDecode.headfile")
    gunniespath = "~/gunnies/"
    recenter = 0
    #SAMBA_prep_folder = os.path.join(SAMBA_mainpath, SAMBA_projectname+"-inputs")
    SAMBA_prep_folder = os.path.join(mainpath, "mouse","ADDeccode_symlink_pool")
    atlas_labels = os.path.join(mainpath, "atlas","IITmean_RPI","IITmean_RPI_labels.nii.gz")

    DTC_DWI_folder = os.path.join(mainpath, "..","ADdecode.01","Analysis","DWI")
    DTC_labels_folder = os.path.join(mainpath, "..","ADdecode.01","Analysis","DWI")
    DTC_transforms = os.path.join(mainpath, "..","ADdecode.01","Analysis","Transforms")

    SAMBA_label_folder = os.path.join(SAMBA_mainpath, SAMBA_projectname+"-results", "connectomics")
    orient_string = os.path.join(SAMBA_prep_folder,"relative_orientation.txt")
    superpose=True
    copytype="truecopy"
    overwrite=False
    file_ids=["coreg", "fa", "b0", "bval", "bvec", "mask"]
    subjects = ["S02654", "S02666",  "S02670",  "S02686", "S02690", "S02695", "S02720", "S02737", "S02753", "S02765", "S02781", "S02802",
                "S02804", "S02813", "S02817", "S02840", "S02877", "S02898", "S02926", "S02938", "S02939", "S02954", "S02967",
                "S02987", "S03010", "S03017", "S03033", "S03034", "S03045", "S03048", "S02771", "S02842", "S02812", "S02871", "S02715","S03069"]

elif project == "APOE":

    SAMBA_mainpath = os.path.join(mainpath, "mouse")
    SAMBA_projectname = "VBM_20APOE01_chass_symmetric3_allAPOE"
    SAMBA_headfile = os.path.join(SAMBA_headfile_dir, "jas297_SAMBA_APOE.headfile")
    gunniespath = "~/gunnies/"
    recenter = 0
    # SAMBA_prep_folder = os.path.join(SAMBA_mainpath, SAMBA_projectname+"-inputs")
    SAMBA_prep_folder = os.path.join(mainpath, "19abb14/")
    atlas_labels = os.path.join(mainpath,"atlases","chass_symmetric3","chass_symmetric3_labels.nii.gz")

    DTC_DWI_folder = os.path.join(mainpath,"mouse","APOE_series","DWI")
    DTC_labels_folder = os.path.join(mainpath,"mouse","APOE_series","DWI")
    DTC_DWI_folder = "samos.dhe.duke.edu:/mnt/paros_MRI/jacques/APOE/DWI_allsubj/"
    DTC_labels_folder = "samos.dhe.duke.edu:/mnt/paros_MRI/jacques/APOE/DWI_allsubj/"

    SAMBA_label_folder = os.path.join(SAMBA_mainpath, SAMBA_projectname + "-results", "connectomics")
    SAMBA_work_folder = os.path.join(SAMBA_mainpath, SAMBA_projectname + "-work")
    orient_string = os.path.join(SAMBA_prep_folder, "relative_orientation.txt")
    superpose = False
    copytype = "truecopy"
    overwrite = False
    file_ids = ["coreg", "fa", "b0", "bval", "bvec", "mask"]
    subjects = ['N57437', 'N57442', 'N57446', 'N57447','N57449','N57451','N57496','N57498','N57500','N57502','N57504', 'N57513',
                'N57515','N57518','N57520','N57522','N57546','N57548','N57550','N57552','N57554','N57559','N57580','N57582','N57584',
                'N57587','N57590','N57692','N57694','N57700','N57500','N57702','N57709',
                "N58214", "N58215", "N58216", "N58217", "N58218", "N58219", "N58221", "N58222", "N58223", "N58224",
                "N58225", "N58226", "N58228",
                "N58229", "N58230", "N58231", "N58232", "N58633", "N58634", "N58635", "N58636", "N58649", "N58650",
                "N58651", "N58653", "N58654",
                'N58408', 'N58398', 'N58714', 'N58740', 'N58477', 'N58734', 'N58309', 'N58792', 'N58302',
                'N58784', 'N58706', 'N58361', 'N58355', 'N58712', 'N58790', 'N58606', 'N58350', 'N58608',
                'N58779', 'N58500', 'N58604', 'N58749', 'N58510', 'N58394', 'N58346', 'N58344', 'N58788', 'N58305',
                'N58514', 'N58794', 'N58733', 'N58655', 'N58735', 'N58310', 'N58400', 'N58708', 'N58780', 'N58512',
                'N58747', 'N58303', 'N58404', 'N58751', 'N58611', 'N58745', 'N58406', 'N58359', 'N58742', 'N58396',
                'N58613', 'N58732', 'N58516', 'N58402']


    subjects = ["N58214", "N58215", "N58216", "N58217", "N58218", "N58219", "N58221", "N58222", "N58223", "N58224",
                "N58225", "N58226", "N58228",
                "N58229", "N58230", "N58231", "N58232", "N58633", "N58634", "N58635", "N58636", "N58649", "N58650",
                "N58651", "N58653", "N58654",
                'N58408', 'N58398', 'N58714', 'N58740', 'N58477', 'N58734', 'N58309', 'N58792', 'N58302',
                'N58784', 'N58706', 'N58361', 'N58355', 'N58712', 'N58790', 'N58606', 'N58350', 'N58608',
                'N58779', 'N58500', 'N58604', 'N58749', 'N58510', 'N58394', 'N58346', 'N58344', 'N58788', 'N58305',
                'N58514', 'N58794', 'N58733', 'N58655', 'N58735', 'N58310', 'N58400', 'N58708', 'N58780', 'N58512',
                'N58747', 'N58303', 'N58404', 'N58751', 'N58611', 'N58745', 'N58406', 'N58359', 'N58742', 'N58396',
                'N58613', 'N58732', 'N58516', 'N58402']
    # subject 'N58610' 'N58612' 'N58813' retired, back on SAMBA_prep, to investigate
else:
    raise Exception("Unknown project name")

for subject in subjects:
    create_backport_labels(subject, SAMBA_mainpath, SAMBA_projectname, atlas_labels, orient_string, headfile = SAMBA_headfile, overwrite=overwrite)

overwrite=False

remote=False

if "." and ":" in DTC_DWI_folder:
    import paramiko
    if "@" in DTC_DWI_folder:
        DTC_DWI_folder_split = DTC_DWI_folder.split("@")
        username = DTC_DWI_folder_split[0]
        server = DTC_DWI_folder_split[1].split(".")[0]
        password = getpass.getpass()
    else:
        server = DTC_DWI_folder.split(".")[0]
        username = getpass.getuser()
        password = getpass.getpass()
        DTC_DWI_folder_split = username + "@" + DTC_DWI_folder
    DTC_DWI_folder = DTC_DWI_folder.split(":")[1]
    ssh = paramiko.SSHClient()
    ssh.load_host_keys(os.path.expanduser(os.path.join("~", ".ssh", "known_hosts")))
    ssh.connect(server, username=username, password=password)
    DTC_labels_folder = DTC_DWI_folder
    remote=True

if remote:
    sftp = ssh.open_sftp()

overwrite=True
for filename in os.listdir(SAMBA_prep_folder):
    if any(x in filename for x in file_ids) and any(x in filename for x in subjects):
        filepath=os.path.join(SAMBA_prep_folder,filename)
        if Path(filepath).is_symlink():
            filename=Path('A').resolve()
        filenewpath = os.path.join(DTC_DWI_folder, filename)
        if not os.path.isfile(filenewpath) or overwrite:
            if copytype=="shortcut":
                if remote:
                    raise Exception("Can't build shortcut to remote path")
                else:
                    buildlink(filepath, filenewpath)
            elif copytype=="truecopy":
                if remote:
                    if not overwrite:
                        try:
                            sftp.stat(filenewpath)
                            if verbose:
                                print(f'file at {filenewpath} exists')
                        except IOError:
                            if verbose:
                                print(f'copying file {filepath} to {filenewpath}')
                            sftp.put(filepath, filenewpath)
                    else:
                        if verbose:
                            print(f'copying file {filepath} to {filenewpath}')
                        sftp.put(filepath, filenewpath)
                else:
                    shutil.copy(filepath, filenewpath)


template_type_prefix = os.path.basename(os.path.dirname(glob.glob(os.path.join(SAMBA_work_folder,"dwi","SyN*/"))[0]))
template_runs = glob.glob((os.path.join(SAMBA_work_folder,"dwi",template_type_prefix,"*/")))
mymax=-1
for template_run in template_runs:
    if "NoNameYet" in template_run and template_run[-4:-2]=="_i":
        if int(template_run[-2])>mymax:
            mymax=int(template_run[-2])
            final_template_run=template_run
if mymax==-1:
    raise Exception(f"Could not find template runs in {os.path.join(mainpath, f'{SAMBA_projectname}-work','dwi',template_type_prefix)}")

for subject in subjects:
    subjectpath = glob.glob(os.path.join(SAMBA_label_folder, f'{subject}/'))
    if np.size(subjectpath) == 1:
        subjectpath = subjectpath[0]
    elif np.size(subjectpath) > 1:
        raise Exception('Too many subject folders')
    else:
        subjectpath = SAMBA_label_folder

    labelspath = glob.glob(os.path.join(subjectpath, f'{subject}*_preprocess_labels.nii*'))
    if np.size(labelspath) == 1:
        labelspath = labelspath[0]
    else:
        raise Exception(f"Could not find file at {os.path.join(subjectpath, f'{subject}*_labels.nii*')}")
    newlabelspath = os.path.join(DTC_labels_folder,f'{subject}_labels.nii.gz')
    if not os.path.exists(newlabelspath) or overwrite:
        if remote:
            if not overwrite:
                try:
                    sftp.stat(newlabelspath)
                    if verbose:
                        print(f'file at {newlabelspath} exists')
                except IOError:
                    if verbose:
                        print(f'copying file {labelspath} to {newlabelspath}')
                    sftp.put(labelspath, newlabelspath)
            else:
                sftp.put(filepath, filenewpath)
        else:
            shutil.copy(labelspath, newlabelspath)
    else:
        print(f"File already exists at {newlabelspath}")
if remote:
    sftp.close()
    ssh.close()


for subject in subjects:
    trans = os.path.join(SAMBA_work_folder, "preprocess", "base_images", "translation_xforms",
                         f"{subject}_0DerivedInitialMovingTranslation.mat")
    rigid = os.path.join(SAMBA_work_folder, "dwi", f"{subject}_rigid.mat")
    affine = os.path.join(SAMBA_work_folder, "dwi", f"{subject}_affine.mat")
    runno_to_MDT = os.path.join(final_template_run, "reg_diffeo", f"{subject}_to_MDT_warp.nii.gz")

    transform_files = [trans, rigid, affine, runno_to_MDT]

    for file in transform_files:
        if any(x in filename for x in file_ids) and any(x in filename for x in subjects):
            filepath=os.path.join(SAMBA_prep_folder,filename)
            if Path(filepath).is_symlink():
                filename=Path('A').resolve()
            filenewpath = os.path.join(DTC_transforms, filename)
            if not os.path.isfile(filenewpath) or overwrite:
                if copytype=="shortcut":
                    if remote:
                        raise Exception("Can't build shortcut to remote path")
                    else:
                        buildlink(filepath, filenewpath)
                elif copytype=="truecopy":
                    if remote:
                        if not overwrite:
                            try:
                                sftp.stat(filenewpath)
                                if verbose:
                                    print(f'file at {filenewpath} exists')
                            except IOError:
                                if verbose:
                                    print(f'copying file {filepath} to {filenewpath}')
                                sftp.put(filepath, filenewpath)
                        else:
                            if verbose:
                                print(f'copying file {filepath} to {filenewpath}')
                            sftp.put(filepath, filenewpath)
                    else:
                        shutil.copy(filepath, filenewpath)


# bash_label_maker = "/Volumes/Data/Badea/Lab/mouse/create_backported_labels_for_fiber_looper.bash"
# mainpath = "/Volumes/Data/Badea/Lab/mouse/"
# gunniespath = "/Users/alex/bass/gitfolder/wuconnectomes/gunnies/"
# project_name = "VBM_21ADDecode03_IITmean_RPI_fullrun"
# atlas_labels = "/Volumes/Data/Badea/Lab/atlas/IITmean_RPI/IITmean_RPI_labels.nii.gz"

#copytype = [shortcut, truecopy]