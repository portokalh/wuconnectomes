
"""
Created by Jacques Stout
Small useful tools for files and data management

"""

import os
import glob
import numpy as np
from pathlib import Path
import warnings
import subprocess, pipes
from scipy.io import loadmat

def mkcdir(folderpaths):
    #creates new folder only if it doesnt already exists
    if np.size(folderpaths) == 1:
        if not os.path.exists(folderpaths):
            os.mkdir(folderpaths)
    else:
        for folderpath in folderpaths:
            if not os.path.exists(folderpath):
                os.mkdir(folderpath)

def check_for_duplicates(checklist,checklist2):
    lenlist = np.size(checklist)
    lenlist2 = np.size(checklist)
    samelist=False
    if np.all(checklist==checklist2):
        samelist=True
    for i in np.arange(lenlist):
        for j in np.arange(lenlist2):
            if checklist[i]==checklist2[j] and (i!=j or not samelist):
                if samelist:
                    print(f"Found the value {checklist[i]} at {i} and {j}")
                else:
                    print(f"Found the value {checklist[i]} in first list at {i} and second list at {j}")

def get_single_var_matlab(path):
    mat_struct = loadmat(path)
    var_name = list(mat_struct.keys())[0]
    if np.size(var_name) >2:
        raise Exception('Multiple values in struct, load it properly with scipy.io.loadmat')
    myvar = mat_struct[var_name]
    return myvar

def exists_remote(host, path):
    """Test if a file exists at path on a host accessible with SSH."""
    status = subprocess.call(
        ['ssh', host, 'test -f {}'.format(pipes.quote(path))])
    if status == 0:
        return True
    if status == 1:
        return False
    raise Exception('SSH failed')

def file_rename(folder, initstring, finalstring, identifier_string="*", anti_identifier_string='the answer is obv 42'):
    #Renames files in a folder by replacing init string with final string, identifier_string being an optional
    #identifier (* format)
    files = glob.glob(os.path.join(folder, identifier_string))
    for myfile in files:
        filename = os.path.basename(myfile)
        newfilename = filename.replace(initstring, finalstring)
        newfilepath = os.path.join(folder, newfilename)
        if newfilepath!=myfile and anti_identifier_string not in filename:
            os.rename(myfile, newfilepath)
            #print(myfile, newfilepath)

def check_files(files):
    exists=[]
    newfiles = []
    for file in files:
        if "*" in file:
            testfile = glob.glob(file)
            if np.size(testfile) == 1:
                exists.append(1)
                newfiles.append(testfile[0])
            elif np.size(testfile) == 0:
                print(f"{file} does not exist")
                exists.append(0)
                newfiles.append("")
            else:
                raise warnings.warn('Too many files of similar names, will take first one but might cause problems!!!')
                exists.append(np.size())
                newfiles.append(testfile[0])
        else:
            if not os.path.exists(file):
                print(f"{file} does not exist")
                exists.append(0)
                newfiles.append("")
            else:
                exists.append(1)
                newfiles.append(file)
    return newfiles, exists


def largerfile(path, identifier=""):
    max_size=0
    max_file=None
    if os.path.isdir(path):
        for folder, subfolders, files in os.walk(path):

            # checking the size of each file
            for file in files:
                if identifier in file:
                    size = os.stat(os.path.join(folder, file)).st_size

                    # updating maximum size
                    if size > max_size:
                        max_size = size
                        max_file = os.path.join(folder, file)
    else:
        if identifier != "":
            files = glob.glob(os.path.join(path,"*"+identifier+"*"))
        else:
            files = glob.glob(path)
        for file in files:
            size = os.stat(file).st_size

            # updating maximum size
            if size > max_size:
                max_size = size
                max_file = file

    return max_file


def getrelativepath(destination, origin):
    if not os.path.isdir(origin):
        origin=os.path.dirname(origin)
    origin      = Path(origin).resolve()
    destination = Path(destination).resolve()
    return(os.path.relpath(destination, start=origin))


def buildlink(real_file, linked_file):
    if os.path.islink(linked_file) and not os.path.exists(os.readlink(linked_file)):
        os.unlink(linked_file)
    if not os.path.islink(linked_file) and os.path.isfile(real_file):
        if os.path.islink(real_file):
            real_file = os.readlink(real_file)
        relpath = getrelativepath(real_file, linked_file)
        link_cmd=f"ln -s ./{relpath} {linked_file}"
        os.system(link_cmd)

def buildlink_folder(folder, initstring, finalstring, identifier_string="*"):
    files = glob.glob(os.path.join(folder, identifier_string))
    for myfile in files:
        filename = os.path.basename(myfile)
        newfilename = filename.replace(initstring, finalstring)
        newfilepath = os.path.join(folder, newfilename)
        if newfilepath != myfile:
            buildlink(myfile, newfilepath)
            #print(myfile, newfilepath)

def getext(file):
    filesplit=file.split('.')
    ext=""
    for i in np.arange(1,np.size(filesplit)):
        ext=ext+"."+str(filesplit[i])
    return ext


def splitpath(filepath):
    dirname, filename = os.path.split(os.path.abspath(filepath))
    ext = getext(filename)
    filename = filename.split('.')[0]
    return dirname, filename, ext
