import glob, shutil
import numpy as np
import os, re, struct, xlsxwriter, warnings

def get_subjnum(headfile):

    with open(headfile, 'rb') as source:
        if verbose: print('INFO    : Extracting animal number')
        header_size = source.read(4)
        header_size = struct.unpack('I', header_size)
        if verbose: print('INFO    : Header size = ', int(header_size[0]))
        for line in source:
            pattern1 = 'U_specid'
            rx1 = re.compile(pattern1, re.IGNORECASE | re.MULTILINE | re.DOTALL)

            for a in rx1.findall(str(line)):
                animal_ID = str(line).split('=')[1]
                animal_ID = animal_ID.split('\\')[0]

    return(animal_ID)


rawdata_paths = ["/Volumes/dusom_civm-atlas/19.abb.14/","/Volumes/dusom_civm-atlas/20.abb.15/", "/Volumes/dusom_civm-atlas/18.abb.11/"]
rawdata_path = '/Volumes/dusom_civm-atlas/19.abb.14/'
rawdata_res_path = os.path.join(rawdata_path,'research')
outpath = '/Users/jas/jacques/APOE_runno_info/'
outpath_excel = os.path.join(outpath, 'runno_to_animalID_19abb.xlsx')
verbose= True
overwrite=True

subjects_folders = glob.glob(os.path.join(rawdata_res_path,'diffusion*/'))
subjects = []
for subject_folder in subjects_folders:
    subjects.append(subject_folder.split('diffusion')[1][:6])

subjects = sorted(subjects)
if not overwrite and os.path.exists(outpath_excel):
    raise Exception(f'Already created the outpath excel file {outpath_excel}')
if os.path.exists(outpath_excel):
    os.remove(outpath_excel)
workbook = xlsxwriter.Workbook(outpath_excel)
worksheet = workbook.add_worksheet()

i=0
for subject in subjects:
    max_size = 0
    subjectpath = glob.glob(os.path.join(os.path.join(rawdata_res_path, "diffusion*" + subject + "*")))[0]
    subject_headfile_temp = os.path.join(outpath, f'{subject}.headfile')
    if not os.path.exists(subject_headfile_temp):
        subject_headfile = glob.glob(os.path.join(subjectpath, f'archived_diffusion*{subject}*.headfile'))
        if np.size(subject_headfile)>0:
            subject_headfile = subject_headfile[0]
        else:
            subject_headfile = glob.glob(os.path.join(subjectpath, f'diffusion*{subject}*.headfile'))
            if np.size(subject_headfile)>0:
                subject_headfile = subject_headfile[0]
            else:
                txt = f'did not find subject {subject}'
        subject_headfile_temp = os.path.join(outpath, f'{subject}.headfile')
        shutil.copy(subject_headfile, subject_headfile_temp)
    animal_ID = get_subjnum(subject_headfile_temp)
    #os.remove(subject_headfile_temp)
    worksheet.write(i,0,subject)
    worksheet.write(i,1,animal_ID)
    i+=1

workbook.close()
