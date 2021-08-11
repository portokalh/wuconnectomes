
"""
Created by Jacques Stout
Small useful tools for files and data management

"""

import os
import glob
import numpy as np
from pathlib import Path

def mkcdir(folderpaths):
    #creates new folder only if it doesnt already exists
    if np.size(folderpaths) == 1:
        if not os.path.exists(folderpaths):
            os.mkdir(folderpaths)
    else:
        for folderpath in folderpaths:
            if not os.path.exists(folderpath):
                os.mkdir(folderpath)


def file_rename(folder, initstring, finalstring, identifier_string="*"):
    #Renames files in a folder by replacing init string with final string, identifier_string being an optional
    #identifier (* format)
    files = glob.glob(os.path.join(folder, identifier_string))
    for file in files:
        newfile = file.replace(initstring, finalstring)
        if newfile!=file:
            os.rename(file, newfile)


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


def buildlink(linked_file, real_file):
    if os.path.islink(linked_file) and not os.path.exists(os.readlink(linked_file)):
        os.unlink(linked_file)
    if not os.path.islink(linked_file) and os.path.isfile(real_file):
        relpath = getrelativepath(real_file, linked_file)
        link_cmd=f"ln -s ./{relpath} {linked_file}"
        os.system(link_cmd)


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