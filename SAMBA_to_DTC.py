from pathlib import Path
import os
import shutil
from file_tools import buildlink, mkcdir, getfromfile
import glob
import numpy as np
import warnings
from create_backported_labels import create_backport_labels
import os
import getpass
from computer_nav import get_mainpaths, checkfile_exists_remote
from convert_atlas_mask import convert_labelmask, atlas_converter

#project = ["AD_Decode", "APOE"]
#project = "APOE"
#project = "AMD"
project = 'APOE'
#project = 'AMD'
#project = 'Chavez'
verbose = True
mainpath = "/Volumes/Data/Badea/Lab/"
#mainpath = "/mnt/munin6/Badea/Lab/"
SAMBA_headfile_dir = os.path.join(mainpath, "samba_startup_cache")
file_ids = ["coreg", "subjspace_fa", "subjspace_b0", "bval", "bvec", "subjspace_mask", "reference", "subjspace_dwi", "relative_orientation"]
file_ids = ['coreg']
#if project == 'AMD':
#    file_ids = ["relative_orientation"]


if project == 'Chavez':
    SAMBA_mainpath = os.path.join(mainpath, "mouse")

    #SAMBA_projectname = "VBM_20APOE01_chass_symmetric3_allAPOE"
    SAMBA_projectname = "VBM_21Chavez01_chass_symmetric3_all"

    SAMBA_headfile = os.path.join(SAMBA_headfile_dir, "rja20_hm190_SAMBA_Chavez.headfile")

    gunniespath = "~/gunnies/"
    recenter = 0
    # SAMBA_prep_folder = os.path.join(SAMBA_mainpath, SAMBA_projectname+"-inputs")
    #SAMBA_prep_folder = os.path.join(mainpath, "APOE_symlink_pool")
    #SAMBA_prep_folder = os.path.join(mainpath, '19abb14')
    SAMBA_prep_folder = os.path.join(mainpath, "mouse","Chavez_symlink_pool_allfiles")

    atlas_labels = os.path.join(mainpath,"atlases","chass_symmetric3","chass_symmetric3_labels.nii.gz")
    atlas_legends = os.path.join(mainpath,'atlases/CHASSSYMM3AtlasLegends.xlsx')
    #DTC_DWI_folder = "samos.dhe.duke.edu:/mnt/paros_DB/Projects/21.chavez.01/Analysis/DWI_allsubj_2/"
    #DTC_labels_folder = "samos.dhe.duke.edu:/mnt/paros_DB/Projects/21.chavez.01/Analysis/DWI_allsubj_2/"
    DTC_DWI_folder = "DWI_allsubj/"
    DTC_labels_folder = "DWI_allsubj/"

    SAMBA_label_folder = os.path.join(SAMBA_mainpath, SAMBA_projectname + "-results", "connectomics")
    SAMBA_work_folder = os.path.join(SAMBA_mainpath, SAMBA_projectname + "-work")
    orient_string = os.path.join(SAMBA_prep_folder, "relative_orientation.txt")
    superpose = False
    copytype = "truecopy"
    overwrite = False
    preppath = None
    subjects = ['C_20220124_001', 'C_20220124_002', 'C_20220124_003', 'C_20220124_004', 'C_20220124_005', 'C_20220124_006', 'C_20220124_007']
    subjects = ['C_20220124_001', 'C_20220124_002', 'C_20220124_003', 'C_20220124_004', 'C_20220124_005', 'C_20220124_006', 'C_20220124_007']

    #subjects_folders = glob.glob(os.path.join(SAMBA_work_folder, '*affine.mat/'))

    """
    subjects_all = glob.glob(os.path.join(SAMBA_work_folder, 'preprocess','*_dwi_masked.nii.gz'))
    subjects = []
    for subject in subjects_all:
        subject_name = os.path.basename(subject)
        subjects.append(subject_name[:6])
    """

    #removed_list = ['C_20220124_004'] #004 had a terrible orientation and would have to be individually treated in order to be on par with the others
    removed_list = []
    for remove in removed_list:
        if remove in subjects:
            subjects.remove(remove)


elif project == "AD_Decode":

    SAMBA_mainpath = os.path.join(mainpath, "mouse")
    SAMBA_projectname = "VBM_21ADDecode03_IITmean_RPI_fullrun"
    SAMBA_headfile = os.path.join(SAMBA_headfile_dir, "jas297_SAMBA_ADDecode.headfile")
    gunniespath = "~/gunnies/"
    recenter = 0
    #SAMBA_prep_folder = os.path.join(SAMBA_mainpath, SAMBA_projectname+"-inputs")
    SAMBA_prep_folder = os.path.join(mainpath, "human","ADDeccode_symlink_pool_allfiles")
    atlas_labels = os.path.join(mainpath, "atlas","IITmean_RPI","IITmean_RPI_labels.nii.gz")
    atlas_legends = os.path.join(mainpath, "/atlases/IITmean_RPI/IITmean_RPI_index.xlsx")
    DTC_DWI_folder = os.path.join(mainpath, "..","ADdecode.01","Analysis","DWI")
    DTC_labels_folder = os.path.join(mainpath, "..","ADdecode.01","Analysis","DWI")
    DTC_transforms = os.path.join(mainpath, "..","ADdecode.01","Analysis","Transforms")

    #DTC_DWI_folder = "samos.dhe.duke.edu:/mnt/paros_MRI/jacques/AD_Decode/Analysis/DWI"
    #DTC_labels_folder = "samos.dhe.duke.edu:/mnt/paros_MRI/jacques/AD_Decode/Analysis/DWI"

    DTC_DWI_folder = "DWI"
    DTC_labels_folder = "DWI"

    SAMBA_work_folder = os.path.join(SAMBA_mainpath, SAMBA_projectname + "-work")
    SAMBA_label_folder = os.path.join(SAMBA_mainpath, SAMBA_projectname+"-results", "connectomics")
    orient_string = os.path.join(SAMBA_prep_folder,"relative_orientation.txt")
    superpose=True
    copytype="truecopy"
    overwrite=False
    preppath = None
    subjects = ["S02654", "S02666",  "S02670",  "S02686", "S02690", "S02695",  "S02715", "S02720", "S02737", "S02753", "S02765", "S02771", "S02781", "S02802",
                "S02804", "S02813", "S02817", "S02840", "S02877", "S02898", "S02926", "S02938", "S02939", "S02954", "S02967",
                "S02987", "S03010", "S03017", "S03033", "S03034", "S03045", "S03048"]
    subjects = ["S02802","S01912", "S02110", "S02224", "S02227", "S02230", "S02231", "S02266", "S02289", "S02320", "S02361", "S02363", "S02373", "S02386", "S02390", "S024S02", "S02410", "S02421", "S02424", "S02446", "S02451", "S02469", "S02473", "S02485", "S02490", "S02491", "S02506","S02524","S02535","S02690","S02715","S02771","S02804","S02817", "S02840","S02877","S02898","S02926","S02938","S02939","S02954", "S03017", "S03028", "S03048","S02524","S02535","S02690","S02715","S02771","S02804","S02812","S02817", "S02840","S02871","S02877","S02898","S02926","S02938","S02939","S02954", "S03017", "S03028", "S03048", "S03069","S02817"]

    subjects = ['S01912', 'S02110', 'S02224', 'S02227', 'S02230', 'S02231', 'S02266', 'S02289', 'S02320', 'S02361', 'S02363',
                'S02373', 'S02386', 'S02390', 'S02402', 'S02410', 'S02421', 'S02424', 'S02446', 'S02451', 'S02469', 'S02473',
                'S02485', 'S02491', 'S02490', 'S02506', 'S02523', 'S02524', 'S02535', 'S02654', 'S02666', 'S02670', 'S02686',
                'S02690', 'S02695', 'S02715', 'S02720', 'S02737', 'S02745', 'S02753', 'S02765', 'S02771', 'S02781', 'S02802',
                'S02804', 'S02813', 'S02812', 'S02817', 'S02840', 'S02842', 'S02871', 'S02877', 'S02898', 'S02926', 'S02938',
                'S02939', 'S02954', 'S02967', 'S02987', 'S03010', 'S03017', 'S03028', 'S03033', 'S03034', 'S03045', 'S03048',
                'S03069', 'S03225', 'S03265', 'S03293', 'S03308', 'S03321', 'S03343', 'S03350', 'S03378', 'S03391', 'S03394']

    subjects = ["S03343", "S03350", "S03378", "S03391", "S03394"]

    #subjects = ["S02695", "S02686", "S02670", "S02666", "S02654","S03010"]
    subjects = ["S01912"]
    subjects = ['S01912', 'S02110', 'S02224', 'S02227', 'S02230', 'S02231', 'S02266', 'S02289', 'S02320', 'S02361',
                'S02363',
                'S02373', 'S02386', 'S02390', 'S02402', 'S02410', 'S02421', 'S02424', 'S02446', 'S02451', 'S02469',
                'S02473',
                'S02485', 'S02491', 'S02490', 'S02506']
    removed_list = ['S02230', 'S02490', 'S02745']
    for remove in removed_list:
        if remove in subjects:
            subjects.remove(remove)

    #oldsubjects = ["01912", "02110", "02224", "02227", "02230", "02231", "02266", "02289", "02320", "02361",
    # #           "02363", "02373", "02386", "02390", "02402", "02410", "02421", "02424", "02446", "02451",
    #            "02469", "02473", "02485", "02490", "02491", "02506"]
    #subjects = []
    #for subject in oldsubjects:
    #    subjects.append('S'+subject)
    # 02842, 03028 has apparently a 92 stack ? to investigate
    #"S02812, , "S02871""  "S03069" has a problem? not prepped

elif project == "APOE":

    SAMBA_mainpath = os.path.join(mainpath, "mouse")
    SAMBA_projectname = "VBM_20APOE01_chass_symmetric3_allAPOE"
    SAMBA_headfile = os.path.join(SAMBA_headfile_dir, "jas297_SAMBA_APOE.headfile")
    gunniespath = "~/gunnies/"
    recenter = 0
    # SAMBA_prep_folder = os.path.join(SAMBA_mainpath, SAMBA_projectname+"-inputs")
    #SAMBA_prep_folder = os.path.join(mainpath, "APOE_symlink_pool")
    #SAMBA_prep_folder = os.path.join(mainpath, '19abb14')
    SAMBA_prep_folder = os.path.join(mainpath, "mouse","APOE_symlink_pool_allfiles")

    atlas_labels = os.path.join(mainpath,"atlases","chass_symmetric3","chass_symmetric3_labels.nii.gz")
    atlas_legends = os.path.join(mainpath,'atlases/CHASSSYMM3AtlasLegends.xlsx')

    DTC_DWI_folder = os.path.join(mainpath,"mouse","APOE_series","DWI")
    DTC_labels_folder = os.path.join(mainpath,"mouse","APOE_series","DWI")
    #DTC_DWI_folder = "samos.dhe.duke.edu:/mnt/paros_MRI/jacques/APOE/DWI_allsubj/"
    #DTC_labels_folder = "samos.dhe.duke.edu:/mnt/paros_MRI/jacques/APOE/DWI_allsubj/"
    DTC_DWI_folder = "DWI_allsubj/"
    DTC_labels_folder = "DWI_allsubj/"

    SAMBA_label_folder = os.path.join(SAMBA_mainpath, SAMBA_projectname + "-results", "connectomics")
    SAMBA_work_folder = os.path.join(SAMBA_mainpath, SAMBA_projectname + "-work")
    orient_string = os.path.join(SAMBA_prep_folder, "relative_orientation.txt")
    superpose = False
    copytype = "truecopy"
    overwrite = False
    preppath = None
    subjects = ['N57437', 'N57442', 'N57446', 'N57447', 'N57449', 'N57451', 'N57496', 'N57498', 'N57500', 'N57502', 'N57504', 'N57513', 'N57515', 'N57518', 'N57520', 'N57522', 'N57546', 'N57548', 'N57550', 'N57552', 'N57554', 'N57559', 'N57580', 'N57582', 'N57584', 'N57587', 'N57590', 'N57692', 'N57694', 'N57700', 'N57702', 'N57709', 'N58302', ' N58303', ' N58305', ' N58309', ' N58310', ' N58344', 'N58346', 'N58350', 'N58355', 'N58359', 'N58361', 'N58394', 'N58396', 'N58398', 'N58400', 'N58402', 'N58404', 'N58406', 'N58408', 'N58477', 'N58500', 'N58510', 'N58512', 'N58514', 'N58516', 'N58604', 'N58606', 'N58608', 'N58610', 'N58611', 'N58613', 'N58706', 'N58708', 'N58712', ' N58214', 'N58215', 'N58216', 'N58217', 'N58218', 'N58219', 'N58221', 'N58222', 'N58223', 'N58224', 'N58225', 'N58226', 'N58228', 'N58229', 'N58230', 'N58231', 'N58232', 'N58633', 'N58634', 'N58635', 'N58636', 'N58649', 'N58650', 'N58651', 'N58653', 'N58654', ' N58398', ' N58714', ' N58740', ' N58477', ' N58734', ' N58792', ' N58302', ' N58784', ' N58706', ' N58361', ' N58355', ' N58712', ' N58790', ' N58606', ' N58350', ' N58608', ' N58779', ' N58500', ' N58604', ' N58749', ' N58510', ' N58394', ' N58346', ' N58788', ' N58514', ' N58794', ' N58733', ' N58655', ' N58735', ' N58400', ' N58708', ' N58780', ' N58512', ' N58747', ' N58404', ' N58751', ' N58611', ' N58745', ' N58406', ' N58359', ' N58742', ' N58396', ' N58613', ' N58732', ' N58516', ' N58402', ' N58935', ' N59003', ' N58819', ' N58909', ' N58919', ' N58889', ' N59010', ' N58859', ' N58917', ' N58815', ' N58997', ' N58999', ' N58881', ' N58853', ' N58995', ' N58877', ' N58883', ' N58885', ' N58906', ' N58821', ' N58855', ' N58861', ' N58857', ' N58851', ' N58887', ' N58879', ' N58829', ' N58913', ' N58831', ' N58941', ' N58813', ' N58952', ' N59022', ' N59026', ' N59033', ' N59035', ' N59039', ' N59041', ' N59065', ' N59066', ' N59072', ' N59076', ' N59078', ' N59080', ' N59097', ' N59099', ' N59109', ' N59116', ' N59118', ' N59120']
    #subjects_folders = glob.glob(os.path.join(SAMBA_work_folder, '*affine.mat/'))
    subjects_all = glob.glob(os.path.join(SAMBA_work_folder, 'preprocess','*_dwi_masked.nii.gz'))
    subjects = []
    for subject in subjects_all:
        subject_name = os.path.basename(subject)
        subjects.append(subject_name[:6])

    subjects = ['N58952', 'N59022', 'N59026', 'N59033', 'N59035', 'N59039', 'N59041', 'N59065', 'N59066', 'N59072', 'N59076', 'N59080', 'N59097', 'N59099', 'N59109', 'N59116', 'N59118', 'N59120']
    removed_list = ['N58610', 'N58613', 'N58732']
    for remove in removed_list:
        if remove in subjects:
            subjects.remove(remove)


    # subject 'N58610' 'N58612' 'N58813' retired, back on SAMBA_prep, to investigate

elif project == "AMD":

    SAMBA_mainpath = os.path.join(mainpath, "mouse")
    SAMBA_projectname = "VBM_19BrainChAMD01_IITmean_RPI_with_2yr"
    SAMBA_headfile = os.path.join(SAMBA_headfile_dir, "rja20_BrainChAMD.01_with_2yr_SAMBA_startup.headfile")
    gunniespath = "~/gunnies/"
    recenter = 0
    SAMBA_prep_folder = os.path.join(SAMBA_mainpath, "whitson_symlink_pool_allfiles")
    atlas_labels = os.path.join(mainpath, "atlas","IITmean_RPI","IITmean_RPI_labels.nii.gz")

    #DTC_DWI_folder = "samos.dhe.duke.edu:/mnt/paros_MRI/jacques/AMD/DWI_v2/"
    #DTC_labels_folder = "samos.dhe.duke.edu:/mnt/paros_MRI/jacques/AMD/DWI_v2/"
    DTC_DWI_folder = "DWI_v2"
    DTC_labels_folder = "DWI_v2"

    SAMBA_label_folder = os.path.join(SAMBA_mainpath, SAMBA_projectname + "-results", "connectomics")
    SAMBA_work_folder = os.path.join(SAMBA_mainpath, SAMBA_projectname + "-work")
    orient_string = os.path.join(SAMBA_prep_folder, "relative_orientation.txt")
    superpose = False
    copytype = "truecopy"
    overwrite = False

    subjects = ["H29056", "H26578", "H29060", "H26637", "H29264", "H26765", "H29225", "H26660", "H29304", "H26890", "H29556",
                "H26862", "H29410", "H26966", "H29403", "H26841", "H21593", "H27126", "H29618", "H27111", "H29627", "H27164",
                "H29502", "H27100", "H27381", "H21836", "H27391", "H21850", "H27495", "H21729", "H27488", "H21915", "H27682",
                "H21956", "H27686", "H22331", "H28208", "H21990", "H28955", "H29878", "H27719", "H22102", "H27841", "H22101",
                "H27842", "H22228", "H28029", "H22140", "H27852", "H22276", "H27999", "H22369", "H28115", "H22644", "H28308",
                "H22574", "H28377", "H22368", "H28325", "H22320", "H28182", "H22898", "H28748", "H22683", "H28373", "H22536",
                "H28433", "H22825", "H28662", "H22864", "H28698", "H23143", "H28861", "H23157", "H28820", "H23028", "H29002",
                "H23210", "H29020", "H23309", "H29161", "H26841", "H26862", "H26949", "H26966", "H27100", "H27126", "H27163",
                "H27246", "H27488", "H27682", "H27686", "H27719", "H27841", "H27842", "H27852", "H27869", "H27999", "H28029",
                "H28068", "H28208", "H28262", "H28325", "H28820", "H28856", "H28869", "H28955", "H29002", "H29044", "H29089",
                "H29127", "H29161", "H29242", "H29254", "H26578", "H26637", "H26660", "H26745", "H26765", "H26850", "H26880",
                "H26890", "H26958", "H26974", "H27017", "H27111", "H27164", "H27381", "H27391", "H27495", "H27610", "H27640",
                "H27680", "H27778", "H27982", "H28115", "H28308", "H28338", "H28373", "H28377", "H28433", "H28437", "H28463",
                "H28532", "H28662", "H28698", "H28748", "H28809", "H28857", "H28861", "H29013", "H29020", "H29025"]
    #subjects = ["H27999"]
else:
    raise Exception("Unknown project name")

"""
if "." and ":" in DTC_DWI_folder:
    import paramiko
    if "@" in DTC_DWI_folder:
        DTC_DWI_folder_split = DTC_DWI_folder.split("@")
        username = DTC_DWI_folder_split[0]
        username = 'alex'
        server = DTC_DWI_folder_split[1].split(".")[0]
        password = getpass.getpass()
    else:
        server = DTC_DWI_folder.split(".")[0]
        username = getpass.getuser()
        username = 'alex'
        password = getpass.getpass()
        DTC_DWI_folder_split = username + "@" + DTC_DWI_folder
    DTC_DWI_folder = DTC_DWI_folder.split(":")[1]
    ssh = paramiko.SSHClient()
    ssh.load_host_keys(os.path.expanduser(os.path.join("~", ".ssh", "known_hosts")))
    #ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    ssh.connect(server, username=username, password=password)
    DTC_labels_folder = DTC_DWI_folder
    DTC_transforms = os.path.join(DTC_DWI_folder,'../Transforms')
    remote=True
if remote:
    sftp = ssh.open_sftp()
"""

remote=True
if remote:
    username, passwd = getfromfile('/Users/jas/samos_connect.rtf')
else:
    username = None
    passwd = None


_, outpath, _, sftp = get_mainpaths(remote,project = project, username=username,password=passwd)

DTC_DWI_folder = os.path.join(outpath,DTC_DWI_folder)
DTC_labels_folder = os.path.join(outpath,DTC_labels_folder)
DTC_transforms = os.path.join(DTC_DWI_folder,'../Transforms')
mkcdir([DTC_DWI_folder,DTC_labels_folder,DTC_transforms],sftp)

"""
for subject in subjects:
    labelspath_remote = os.path.join(DTC_labels_folder, f'{subject}_labels.nii.gz')
    #if not checkfile_exists_remote(labelspath_remote,sftp) or overwrite:
    create_backport_labels(subject, SAMBA_mainpath, SAMBA_projectname, SAMBA_prep_folder, atlas_labels, headfile = SAMBA_headfile, overwrite=overwrite)
"""

mkcdir([DTC_DWI_folder,DTC_labels_folder],sftp)

for filename in os.listdir(SAMBA_prep_folder):
    if any(x in filename for x in file_ids) and any(x in filename for x in subjects):
        filepath=os.path.join(SAMBA_prep_folder,filename)
        if 'N59010' in filename:
            print('hi')
        if Path(filepath).is_symlink():
            filepath=Path(filepath).resolve()
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
                        try:
                            sftp.put(filepath, filenewpath)
                        except:
                            print('test')
                            os.remove(filepath)
                else:
                    shutil.copy(filepath, filenewpath)


template_type_prefix = os.path.basename(os.path.dirname(glob.glob(os.path.join(SAMBA_work_folder,"dwi","SyN*/"))[0]))
template_runs = glob.glob((os.path.join(SAMBA_work_folder,"dwi",template_type_prefix,"*/")))
mymax = -1
for template_run in template_runs:
    if "NoNameYet" in template_run and template_run[-4:-2]=="_i":
        if int(template_run[-2])>mymax:
            mymax=int(template_run[-2])
            final_template_run=template_run
if mymax==-1:
    for template_run in template_runs:
        if "dwiMDT_Control_n72" in template_run and template_run[-4:-2]=="_i":
            if int(template_run[-2])>mymax:
                mymax=int(template_run[-2])
                final_template_run=template_run
if mymax == -1:
    raise Exception(f"Could not find template runs in {os.path.join(mainpath, f'{SAMBA_projectname}-work','dwi',template_type_prefix)}")

if project != "AMD":
    for subject in subjects:
        subjectpath = glob.glob(os.path.join(SAMBA_label_folder, f'{subject}/'))
        if np.size(subjectpath) == 1:
            subjectpath = subjectpath[0]
        elif np.size(subjectpath) > 1:
            raise Exception('Too many subject folders')
        else:
            subjectpath = SAMBA_label_folder

        labelspath = glob.glob(os.path.join(subjectpath, f'{subject}*labels.nii*'))
        if np.size(labelspath) == 1:
            labelspath = labelspath[0]
        else:
            warnings.warn(f"Could not find file at {os.path.join(subjectpath, f'{subject}*_labels.nii*')}")
            continue
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
                    sftp.put(labelspath, newlabelspath)
                    if verbose:
                        print(f'copying file {labelspath} to {newlabelspath}')

            else:
                shutil.copy(labelspath, newlabelspath)
                if verbose:
                    print(f'copying file {labelspath} to {newlabelspath}')
        else:
            if verbose:
                print(f"File already exists at {newlabelspath}")

        newlabelspath_ordered = newlabelspath.replace('_labels','_labels_lr_ordered')

        if not os.path.exists(newlabelspath_ordered) or overwrite:
            if remote:
                if not overwrite:
                    try:
                        sftp.stat(newlabelspath_ordered)
                        if verbose:
                            print(f'file at {newlabelspath_ordered} exists')
                    except IOError:
                        if verbose:
                            print(f'creating file {newlabelspath_ordered} from {labelspath}')
                        converter_lr, converter_comb, index_to_struct_lr, index_to_struct_comb = atlas_converter(
                            atlas_legends)
                        convert_labelmask(labelspath, converter_lr, atlas_outpath=newlabelspath_ordered, sftp=sftp)
                else:
                    if verbose:
                        print(f'creating file {newlabelspath_ordered} from {labelspath}')
                    converter_lr, converter_comb, index_to_struct_lr, index_to_struct_comb = atlas_converter(
                        atlas_legends)
                    convert_labelmask(labelspath, converter_lr, atlas_outpath=newlabelspath_ordered, sftp=sftp)

            else:
                converter_lr, converter_comb, index_to_struct_lr, index_to_struct_comb = atlas_converter(
                    atlas_legends)
                convert_labelmask(labelspath, converter_lr, atlas_outpath=newlabelspath_ordered)
                if verbose:
                    print(f'creating file {newlabelspath_ordered} from {labelspath}')
        else:
            if verbose:
                print(f"File already exists at {newlabelspath_ordered}")


elif project == "AMD":
    for subject in subjects:
        subjectpath = glob.glob(os.path.join(SAMBA_label_folder, f'{subject}/'))
        if np.size(subjectpath) == 1:
            subjectpath = subjectpath[0]
        elif np.size(subjectpath) > 1:
            raise Exception('Too many subject folders')
        else:
            subjectpath = SAMBA_label_folder

        labelspath = glob.glob(os.path.join(subjectpath, f'{subject}*_IITmean_RPI_labels.nii*'))
        coreg_prepro =  glob.glob(os.path.join(subjectpath, f'{subject}*_nii4D_masked_isotropic.nii*'))
        if np.size(labelspath) == 1:
            labelspath = labelspath[0]
        else:
            warnings.warn(f"Could not find file at {os.path.join(subjectpath, f'{subject}*_labels.nii*')}")
            continue

        if np.size(coreg_prepro) == 1:
            coreg_prepro = coreg_prepro[0]
        else:
            warnings.warn(f"Could not find file at {os.path.join(subjectpath, f'{subject}*_labels.nii*')}")
            continue

        newlabelspath = os.path.join(DTC_labels_folder, f'{subject}_labels.nii.gz')
        new_coregpath = os.path.join(DTC_labels_folder, f'{subject}_coreg_diff.nii.gz')

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
                    sftp.put(labelspath, newlabelspath)
                    if verbose:
                        print(f'copying file {labelspath} to {newlabelspath}')

            else:
                shutil.copy(labelspath, newlabelspath)
                if verbose:
                    print(f'copying file {labelspath} to {newlabelspath}')
        else:
            if verbose:
                print(f"File already exists at {newlabelspath}")

        if not os.path.exists(new_coregpath) or overwrite:
            if remote:
                if not overwrite:
                    try:
                        sftp.stat(new_coregpath)
                        if verbose:
                            print(f'file at {new_coregpath} exists')
                    except IOError:
                        if verbose:
                            print(f'copying file {coreg_prepro} to {new_coregpath}')
                        sftp.put(coreg_prepro, new_coregpath)
                else:
                    sftp.put(coreg_prepro, new_coregpath)
                    if verbose:
                        print(f'copying file {coreg_prepro} to {new_coregpath}')

            else:
                shutil.copy(coreg_prepro, new_coregpath)
                if verbose:
                    print(f'copying file {coreg_prepro} to {new_coregpath}')
        else:
            if verbose:
                print(f"File already exists at {new_coregpath}")

overwrite=False
for subject in subjects:
    trans = os.path.join(SAMBA_work_folder, "preprocess", "base_images", "translation_xforms",
                         f"{subject}_0DerivedInitialMovingTranslation.mat")
    rigid = os.path.join(SAMBA_work_folder, "dwi", f"{subject}_rigid.mat")
    affine = os.path.join(SAMBA_work_folder, "dwi", f"{subject}_affine.mat")
    runno_to_MDT = os.path.join(final_template_run, "reg_diffeo", f"{subject}_to_MDT_warp.nii.gz")

    burn_dir = os.path.join(SAMBA_mainpath, "burn_after_reading")
    affine_mat_path = os.path.join(burn_dir, f'{subject}_affine.txt')
    if not os.path.exists(affine_mat_path) or overwrite:
        cmd = f'ConvertTransformFile 3 {affine} {affine_mat_path} --matrix'
        os.system(cmd)

    transform_files = [trans, rigid, affine, affine_mat_path, runno_to_MDT]

    for filepath in transform_files:
        if os.path.exists(filepath):
            if Path(filepath).is_symlink():
                filepath=Path(filepath).resolve()
            filename = os.path.basename(filepath)
            filenewpath = os.path.join(DTC_transforms, filename)
            try:
                sftp.chdir(DTC_transforms)
            except IOError:
                sftp.mkdir(DTC_transforms)
            if not os.path.isfile(filenewpath) or overwrite:
                if copytype=="shortcut":
                    if remote:
                        raise Exception("Can't build shortcut to remote path")
                    else:
                        buildlink(filepath, filenewpath)
                        if verbose:
                            print(f'Built link for {filepath} at {filenewpath}')
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
                        if verbose:
                            print(f'copying file {filepath} to {filenewpath}')
        else:
            print(f'Could not find {filepath}')

MDT_refs = ['fa', 'md', 'rd', 'ad', 'b0']
for subject in subjects:
    for MDT_ref in MDT_refs:
        filepath = os.path.join(final_template_run, 'reg_images',f'{subject}_{MDT_ref}_to_MDT.nii.gz')
        if not os.path.exists(filepath):
            txt = f'Could not find {filepath} in {final_template_run} for subject {subject} in project {project}'
            warnings.warn(txt)
        else:
            if Path(filepath).is_symlink():
                filepath=Path(filepath).resolve()
            filename = os.path.basename(filepath)
            filenewpath = os.path.join(DTC_labels_folder, filename)
            if not os.path.isfile(filenewpath) or overwrite:
                if copytype=="shortcut":
                    if remote:
                        raise Exception("Can't build shortcut to remote path")
                    else:
                        buildlink(filepath, filenewpath)
                        if verbose:
                            print(f'Built link for {filepath} at {filenewpath}')
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
                        if verbose:
                            print(f'copying file {filepath} to {filenewpath}')
                        shutil.copy(filepath, filenewpath)

if remote:
    sftp.close()
# bash_label_maker = "/Volumes/Data/Badea/Lab/mouse/create_backported_labels_for_fiber_looper.bash"
# mainpath = "/Volumes/Data/Badea/Lab/mouse/"
# gunniespath = "/Users/alex/bass/gitfolder/wuconnectomes/gunnies/"
# project_name = "VBM_21ADDecode03_IITmean_RPI_fullrun"
# atlas_labels = "/Volumes/Data/Badea/Lab/atlas/IITmean_RPI/IITmean_RPI_labels.nii.gz"

#copytype = [shortcut, truecopy]

"""
    APOE
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

"""
