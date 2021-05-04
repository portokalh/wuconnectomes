
import pandas as pd

import numpy as np
import matplotlib.pyplot as plt

def get_matrix_connectomes(filepath):

    #dataframe = pd.read_excel("/Volumes/Data/Badea/Lab/mouse/C57_JS/VBM_whiston_Figs_inclusive/H29060_stepsize_2_all_wholebrain_connectomes.xlsx")
    dataframe = pd.read_excel(filepath)
    ROIlist = dataframe['Unnamed: 0']
    dataframe.pop('Unnamed: 0')
    matrix=dataframe.values[:]
    return matrix,ROIlist

dataframepath = "/Volumes/Data/Badea/Lab/mouse/C57_JS/VBM_whiston_Figs_inclusive/H29060_stepsize_2_all_wholebrain_connectomes.xlsx"
fig, ax = plt.subplots()
matrix,ROIlist=get_matrix_connectomes(dataframepath)

ax.matshow(matrix, cmap=plt.cm.Blues)

print('hi')