
import os
import glob
import re
import shutil
import gzip


topfolder="/Volumes/dusom_dibs_ad_decode/all_staff/APOE_temp/diffusion_prep_locale/"
destination="/Volumes/Data/Badea/Lab/19abb14/"
contrasts = ["dwi","fa"]

for subdir, dirs, files in os.walk(topfolder):
    if "diffusion_prep" in os.path.basename(subdir):
        #print(os.path.basename(subdir))
        skip=False
        temp = re.findall(r'\d+', subdir)
        subj = str(list(map(int, temp))[0])
        mainfile = glob.glob(os.path.join(subdir, 'LPCA_' + subj + '_nii4D.nii.gz'))
        newlocation = os.path.join(destination, 'N'+subj+'_nii4D.nii')
        if len(mainfile)>0 and not os.path.exists(newlocation):
            with gzip.open(mainfile[0], 'rb') as f_in:
                with open(newlocation, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
                    print("Moved the mainfile for subject "+subj+" from "+mainfile[0] + " to "+newlocation)
        elif len(mainfile)==0:
            print("could not find mainfile for subject "+subj+ " in subdir "+str(subdir))
            skip=True
        elif os.path.exists(newlocation):
            print("Subject "+subj+" already saved")
        if not skip:
            for contrast in contrasts:
                #adfile = glob.glob(os.path.join(subdir,'*ad.nii.gz'))
                #rdfile = glob.glob(os.path.join(subdir,'*rd.nii.gz'))
                #b0file = glob.glob(os.path.join(subdir,'*b0.nii.gz'))
                #dwifile = glob.glob(os.path.join(subdir, '*nii4D.nii.gz'))
                #mdfile = glob.glob(os.path.join(subdir, '*md.nii.gz'))
                #maskfile = glob.glob(os.path.join(subdir, '*mask.nii.gz'))
                #dwifile = glob.glob(os.path.join(subdir,'*dwi.nii.gz'))
                #fafile = glob.glob(os.path.join(subdir, '*fa.nii.gz'))
                contrastfile = glob.glob(os.path.join(subdir,str(subj)+'_'+contrast+'.nii.gz'))
                contrastfiletrue = os.readlink(contrastfile[0])
                contrastfiletrue = os.path.join(subdir, contrastfiletrue)
                newlocation = os.path.join(destination, 'N'+str(subj)+'_'+contrast+'.nii.gz')
                if not os.path.exists(newlocation):
                    shutil.copyfile(contrastfiletrue, newlocation)
                    print("Moved the file for contrast "+contrast+" of subject "+subj)
                else:
                    print("Already moved file at "+newlocation)