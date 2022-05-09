
import socket, os, getpass, paramiko, glob
from dipy.io.image import load_nifti
import fnmatch
import numpy as np
import pickle
import nibabel as nib
from nifti_handler import get_reference_info


def get_mainpaths(remote=False, project='any',username=None,password=None):
    computer_name = socket.gethostname()
    project_rename = {'Chavez':'21.chavez.01','AD_Decode':'AD_Decode','APOE':'APOE','AMD':'AMD'}
    samos = False
    if 'samos' in computer_name:
        inpath = '/mnt/paros_MRI/jacques/'
        outpath = '/mnt/paros_MRI/jacques/'
        atlas_folder = '/mnt/paros_MRI/jacques/atlases/'
        remote=False

        #ROI_legends = "/mnt/paros_MRI/jacques/atlases/IITmean_RPI/IITmean_RPI_index.xlsx"
    elif 'santorini' in computer_name or 'hydra' in computer_name:
        # mainpath = '/Users/alex/jacques/'
        inpath = '/Volumes/Data/Badea/Lab/human/'
        outpath = '/Volumes/Data/Badea/Lab/human/'
        atlas_folder = '/Volumes/Data/Badea/Lab/atlases/'
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
            if project == 'Chavez':
                inpath = 'samos.dhe.duke.edu:/mnt/paros_DB/Projects/'
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
        elif ':' in inpath:
            inpath = inpath.split(":")[1]
        outpath=inpath
    if project == 'AD_Decode' or project == 'Chavez':
        outpath = os.path.join(outpath, project_rename[project], 'Analysis')
        inpath = os.path.join(inpath, project_rename[project], 'Analysis')
    else:
        outpath = os.path.join(outpath, project_rename[project])
        inpath = os.path.join(inpath, project_rename[project])

    return inpath, outpath, atlas_folder, sftp

def get_sftp(remote, username=None, password=None):
    computer_name = socket.gethostname()
    server='samos'
    if remote and not 'samos' in computer_name:
        if username is None or password is None:
            username = input("Username:")
            password = getpass.getpass("Password for " + username + ":")
        ssh = paramiko.SSHClient()
        ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
        ssh.load_host_keys(os.path.expanduser(os.path.join("~", ".ssh", "known_hosts")))
        ssh.connect(server, username=username, password=password)
        sftp = ssh.open_sftp()
    else:
        sftp=None

    return sftp

def get_atlas(atlas_folder, atlas_type):
    if atlas_type == 'IIT':
        index_path = os.path.join(atlas_folder,'IITmean_RPI','IITmean_RPI_index.xlsx')
    elif atlas_type == 'CHASSSYMM3':
        index_path = os.path.join(atlas_folder,'CHASSSYMM3AtlasLegends.xlsx')
    else:
        raise Exception('unknown atlas')
    return index_path


def make_temppath(path):
    return f'{os.path.join(os.path.expanduser("~"), os.path.basename(path))}'


def load_nifti_remote(niipath, sftp=None):
    if sftp is None:
        img = nib.load(niipath)
        data = img.get_data()
        vox_size = img.header.get_zooms()[:3]
        affine = img.affine
        header = img.header
        ref_info = get_reference_info(niipath)
    else:
        temp_path = f'{os.path.join(os.path.expanduser("~"), os.path.basename(niipath))}'
        sftp.get(niipath, temp_path)
        try:
            img = nib.load(temp_path)
            data = img.get_data()
            vox_size = img.header.get_zooms()[:3]
            affine = img.affine
            header = img.header
            ref_info = get_reference_info(temp_path)
            os.remove(temp_path)
        except Exception as e:
            os.remove(temp_path)
            raise Exception(e)
    return data, affine, vox_size, header, ref_info


def save_nifti_remote(niiobject,niipath, sftp):

    if sftp is None:
        nib.save(niiobject, niipath)
    else:
        nib.save(niiobject, make_temppath(niipath))
        sftp.put(make_temppath(niipath),niipath)
        os.remove(make_temppath(niipath))
    return


def remove_remote(path, sftp=None):
    if sftp is None:
        os.remove(path)
    else:
        sftp.remove(path)


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


def remote_pickle(picklepath, sftp=None):
    if sftp is not None:
        picklepath_tconnectome = make_temppath(picklepath)
        sftp.get(picklepath, picklepath_tconnectome)
    else:
        picklepath_tconnectome = picklepath
    with open(picklepath_tconnectome, 'rb') as f:
        M = pickle.load(f)
    if sftp is not None:
        os.remove(picklepath_tconnectome)
    return M


def load_trk_remote(trkpath,reference,sftp=None):
    #from dipy.io.streamline import load_trk
    from streamline_nocheck import load_trk as load_trk_spe
    if sftp is not None:
        temp_path = f'{os.path.join(os.path.expanduser("~"), os.path.basename(trkpath))}'
        sftp.get(trkpath, temp_path)
        try:
            trkdata = load_trk_spe(temp_path, reference)
            #trkdata = load_trk(temp_path, reference)
            os.remove(temp_path)
        except Exception as e:
            if os.path.exists(temp_path):
                os.remove(temp_path)
            raise Exception(e)
    else:
        trkdata = load_trk_spe(trkpath, reference)
    return trkdata

def loadmat_remote(matpath, sftp):
    from scipy.io import loadmat
    temp_path = f'{os.path.join(os.path.expanduser("~"), os.path.basename(matpath))}'
    sftp.get(matpath, temp_path)
    try:
        mymat = loadmat(temp_path)
        os.remove(temp_path)
    except Exception as e:
        os.remove(temp_path)
        raise Exception(e)
    return mymat

def glob_remote(path, sftp):
    match_files = []
    if sftp is not None:
        if '.' not in path:
            allfiles = sftp.listdir(path)
            for filepath in allfiles:
                match_files.append(os.path.join(path, filepath))
            return match_files
        else:
            dirpath = os.path.dirname(path)
            try:
                sftp.stat(dirpath)
            except:
                return match_files
            allfiles = sftp.listdir(dirpath)
            #if '*' in path:
            #    for filepath in allfiles:
            #            match_files.append(os.path.join(dirpath,filepath))
            #else:
            for filepath in allfiles:
                if fnmatch.fnmatch(os.path.basename(filepath), os.path.basename(path)):
                    match_files.append(os.path.join(dirpath, filepath))
    else:
        if '.' not in path:
            match_files = glob.glob(path)
        else:
            dirpath = os.path.dirname(path)
            if not os.path.exists(dirpath):
                return(match_files)
            else:
                allfiles = glob.glob(os.path.join(dirpath,'*'))
                for filepath in allfiles:
                    if fnmatch.fnmatch(os.path.basename(filepath), os.path.basename(path)):
                        match_files.append(os.path.join(dirpath, filepath))
    return(match_files)

def pickledump_remote(var,path,sftp=None):
    if sftp is None:
        pickle.dump(var, open(path, "wb"))
    else:
        temp_path = make_temppath(path)
        pickle.dump(var, open(temp_path, "wb"))
        sftp.put(temp_path, path)
        os.remove(temp_path)

def checkfile_exists_remote(path, sftp=None):
    match_files = glob_remote(path,sftp)
    if np.size(match_files)>0:
        return True
    else:
        return False