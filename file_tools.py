import os
import glob
import numpy as np

def mkcdir(folderpaths):
    if np.size(folderpaths) == 1:
        if not os.path.exists(folderpaths):
            os.mkdir(folderpaths)
    else:
        for folderpath in folderpaths:
            if not os.path.exists(folderpath):
                os.mkdir(folderpath)


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