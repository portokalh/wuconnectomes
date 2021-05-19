
"""
Created by Jacques Stout
Small useful tools for files and data management

"""

import os
import glob
import numpy as np

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
        os.rename(file, newfile)


def largerfile(path):
    max_size=0
    max_file=None
    if os.path.isdir(path):
        for folder, subfolders, files in os.walk(path):

            # checking the size of each file
            for file in files:
                size = os.stat(os.path.join(folder, file)).st_size

                # updating maximum size
                if size > max_size:
                    max_size = size
                    max_file = os.path.join(folder, file)
    else:
        files = glob.glob(path)
        for file in files:
            size = os.stat(file).st_size

            # updating maximum size
            if size > max_size:
                max_size = size
                max_file = file

    return max_file

def getext(file):
    filesplit=file.split('.')
    ext=""
    for i in np.arange(1,np.size(filesplit)):
        ext=ext+"."+str(filesplit[i])
    return ext