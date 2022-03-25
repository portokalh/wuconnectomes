
import socket, os, getpass, paramiko
from dipy.io.image import load_nifti
import fnmatch

def get_mainpaths(remote=False, project='any',username=None,password=None):
    computer_name = socket.gethostname()

    samos = False
    if 'samos' in computer_name:
        inpath = '/mnt/paros_MRI/jacques/'
        outpath = '/mnt/paros_MRI/jacques/'
        atlas_folder = '/mnt/paros_MRI/jacques/atlases/'
        #ROI_legends = "/mnt/paros_MRI/jacques/atlases/IITmean_RPI/IITmean_RPI_index.xlsx"
    elif 'santorini' in computer_name or 'hydra' in computer_name:
        # mainpath = '/Users/alex/jacques/'
        inpath = '/Volumes/Data/Badea/Lab/human/'
        outpath = '/Volumes/Data/Badea/Lab/human/'
        atlas_folder = '/Volumes/Data/Badea/ADdecode.01/Analysis/atlases/'
        #ROI_legends = "/Volumes/Data/Badea/ADdecode.01/Analysis/atlases/IITmean_RPI/IITmean_RPI_index.xlsx"
    elif 'blade' in computer_name:
        inpath = '/mnt/munin6/Badea/Lab/human/'
        outpath = '/mnt/munin6/Badea/Lab/human/'
        atlas_folder = '/mnt/munin6/Badea/Lab/atlases/'
       # ROI_legends = "/mnt/munin6/Badea/Lab/atlases/IITmean_RPI/IITmean_RPI_index.xlsx"
    else:
        raise Exception('No other computer name yet')

    sftp = None
    if remote:
        if not 'samos' in computer_name:
            inpath = 'samos.dhe.duke.edu:/mnt/paros_MRI/jacques/'
            if "@" in inpath:
                inpath = inpath.split("@")
                username = inpath[0]
                server = inpath[1].split(".")[0]
                password = getpass.getpass()
            else:
                server = inpath.split(".")[0]
                if username is None or password is None:
                    username = input("Username:")
                    password = getpass.getpass("Password for " + username + ":")
                inpath = username + "@" + inpath
            inpath = inpath.split(":")[1]
            ssh = paramiko.SSHClient()
            ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
            ssh.load_host_keys(os.path.expanduser(os.path.join("~", ".ssh", "known_hosts")))
            ssh.connect(server, username=username, password=password)
            sftp = ssh.open_sftp()
        else:
            inpath = inpath.split(":")[1]
        outpath=inpath
    if project == 'AD_Decode':
        outpath = os.path.join(outpath, project, 'Analysis')
        inpath = os.path.join(inpath, project, 'Analysis')
    else:
        outpath = os.path.join(outpath, project)
        inpath = os.path.join(inpath, project)

    return inpath, outpath, atlas_folder, sftp


def get_atlas(atlas_folder, atlas_type):
    if atlas_type == 'IIT':
        index_path = os.path.join(atlas_folder,'IITmean_RPI','IITmean_RPI_index.xlsx')
    else:
        raise Exception('unknown atlas')
    return index_path


def nii_load_remote(niipath, sftp):
    temp_path = f'{os.path.join(os.path.expanduser("~"), os.path.basename(niipath))}'
    sftp.get(niipath, temp_path)
    import nibabel as nib
    try:
        img = nib.load(temp_path)
        os.remove(temp_path)
    except Exception as e:
        os.remove(temp_path)
        raise Exception(e)
    return img


def load_nifti_remote(niipath, sftp):
    temp_path = f'{os.path.join(os.path.expanduser("~"), os.path.basename(niipath))}'
    sftp.get(niipath, temp_path)
    try:
        data, affine = load_nifti(temp_path)
        os.remove(temp_path)
    except Exception as e:
        os.remove(temp_path)
        raise Exception(e)
    return data, affine


def read_bvals_bvecs_remote(fbvals, fbvecs, sftp):
    from dipy.io.gradients import read_bvals_bvecs
    temp_path_bval = f'{os.path.join(os.path.expanduser("~"), os.path.basename(fbvals))}'
    temp_path_bvec = f'{os.path.join(os.path.expanduser("~"), os.path.basename(fbvecs))}'
    sftp.get(fbvals, temp_path_bval)
    sftp.get(fbvecs, temp_path_bvec)
    try:
        bvals, bvecs = read_bvals_bvecs(temp_path_bval, temp_path_bvec)
        os.remove(temp_path_bval)
        os.remove(temp_path_bvec)
    except Exception as e:
        os.remove(temp_path_bval)
        os.remove(temp_path_bvec)
        raise Exception(e)
    return bvals, bvecs


def load_trk_remote(trkpath,reference,sftp):
    from dipy.io.streamline import load_trk
    temp_path = f'{os.path.join(os.path.expanduser("~"), os.path.basename(trkpath))}'
    sftp.get(trkpath, temp_path)
    try:
        trkdata = load_trk(temp_path, reference)
        os.remove(temp_path)
    except Exception as e:
        os.remove(temp_path)
        raise Exception(e)
    return trkdata


def glob_remote(path, sftp):
    dirpath = os.path.dirname(path)
    allfiles = sftp.listdir(dirpath)
    match_files = []
    for filepath in allfiles:
        if fnmatch.fnmatch(os.path.basename(filepath), os.path.basename(path)):
            match_files.append(filepath)
    return(match_files)
