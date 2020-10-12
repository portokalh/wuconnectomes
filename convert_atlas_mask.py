
import pandas as pd
import numpy as np
from dipy.io.image import load_nifti, save_nifti


atlas_legends = "/Users/alex/jacques/connectomes_testing/atlases/CHASSSYMM3AtlasLegends.xlsx"
df = pd.read_excel(atlas_legends, sheet_name='Sheet1')
df['Structure'] = df['Structure'].str.lower()
index1=df['index']
index2=df['index2']

l = ['N57442', 'N57446', 'N57447','N57449','N57451','N57496','N57498','N57500','N57502','N57504','N57513','N57515','N57518', 'N57520','N57522','N57546','N57447','N57449','N57451','N57496','N57498','N57500','N57502','N57504','N57513','N57515','N57518','N57520','N57522','N57546','N57548', 'N57550', 'N57552', 'N57554', 'N57559', 'N57580', 'N57582', 'N57584', 'N57587', 'N57590', 'N57692', 'N57694', 'N57700', 'N57702', 'N57709']


atlas_folder = "/Volumes/Data/Badea/Lab/mouse/C57_JS/DWI_RAS_40subj/"

for subject in l:
    labels, affine_labels = load_nifti(atlas_folder + subject + "_chass_symmetric3_labels_RAS.nii.gz")
    rslt_whitem = df.loc[df['Subdivisions_7'] == "7_whitematter"]

    labels_leftright = np.copy(labels)
    labels_combinedlr = np.copy(labels)
    #for i in range(int(np.shape(labels)[0]/2-2), int(np.shape(labels)[0])):
    for i in range(np.shape(labels)[0]):
        for j in range(np.shape(labels)[1]):
            for k in range(np.shape(labels)[2]):
                if labels[i,j,k]>1000:
                    labels_leftright[i, j, k] = labels[i, j, k] - 834
                    labels_combinedlr[i, j, k] = labels[i, j, k] - 1000

    save_nifti(atlas_folder + subject + "_chass_symmetric3_labels_RAS_combined.nii.gz", labels_combinedlr, affine_labels)
    save_nifti(atlas_folder + subject + "_chass_symmetric3_labels_RAS_lrordered.nii.gz", labels_leftright, affine_labels)

    print("done")

print("hi")

"""
atlas_legends = "/Users/alex/jacques/connectomes_testing/atlases/CHASSSYMM3AtlasLegends.xlsx"
df = pd.read_excel(atlas_legends, sheet_name='Sheet1')
df['Structure'] = df['Structure'].str.lower()
index1=df['index']
index2=df['index2']

labels, affine_labels = load_nifti("/Users/alex/jacques/connectomes_testing/DWI_RAS_labelsbackup/N57434_chass_symmetric3_labels_RAS.nii.gz")
rslt_whitem = df.loc[df['Subdivisions_7'] == "7_whitematter"]

labels_leftright = np.copy(labels)
labels_combinedlr = np.copy(labels)
#for i in range(int(np.shape(labels)[0]/2-2), int(np.shape(labels)[0])):
for i in range(np.shape(labels)[0]):
    for j in range(np.shape(labels)[1]):
        for k in range(np.shape(labels)[2]):
            if labels[i,j,k]>1000:
                labels_leftright[i, j, k] = labels[i, j, k] - 834
                labels_combinedlr[i, j, k] = labels[i, j, k] - 1000

save_nifti("/Users/alex/jacques/connectomes_testing/DWI_RAS/N57434_chass_symmetric3_labels_RAS_combined.nii.gz", labels_combinedlr, affine_labels)
save_nifti("/Users/alex/jacques/connectomes_testing/DWI_RAS_labelsbackup/N57434_chass_symmetric3_labels_RAS_lrordered.nii.gz", labels_leftright, affine_labels)

print("done")

atlas_legends = "/Users/alex/jacques/connectomes_testing/atlases/CHASSSYMM3AtlasLegends.xlsx"
df = pd.read_excel(atlas_legends, sheet_name='Sheet1')
df['Structure'] = df['Structure'].str.lower()
index1=df['index']
index2=df['index2']

labels, affine_labels = load_nifti("/Users/alex/jacques/connectomes_testing/DWI_RAS_labelsbackup/N57435_chass_symmetric3_labels_RAS.nii.gz")
rslt_whitem = df.loc[df['Subdivisions_7'] == "7_whitematter"]

labels_leftright = np.copy(labels)
labels_combinedlr = np.copy(labels)
#for i in range(int(np.shape(labels)[0]/2-2), int(np.shape(labels)[0])):
for i in range(np.shape(labels)[0]):
    for j in range(np.shape(labels)[1]):
        for k in range(np.shape(labels)[2]):
            if labels[i,j,k]>1000:
                labels_leftright[i, j, k] = labels[i, j, k] - 834
                labels_combinedlr[i, j, k] = labels[i, j, k] - 1000

save_nifti("/Users/alex/jacques/connectomes_testing/DWI_RAS/N57435_chass_symmetric3_labels_RAS_combined.nii.gz", labels_combinedlr, affine_labels)
save_nifti("/Users/alex/jacques/connectomes_testing/DWI_RAS_labelsbackup/N57435_chass_symmetric3_labels_RAS_lrordered.nii.gz", labels_leftright, affine_labels)

print("done")

atlas_legends = "/Users/alex/jacques/connectomes_testing/atlases/CHASSSYMM3AtlasLegends.xlsx"
df = pd.read_excel(atlas_legends, sheet_name='Sheet1')
df['Structure'] = df['Structure'].str.lower()
index1=df['index']
index2=df['index2']

labels, affine_labels = load_nifti("/Users/alex/jacques/connectomes_testing/DWI_RAS_labelsbackup/N57436_chass_symmetric3_labels_RAS.nii.gz")
rslt_whitem = df.loc[df['Subdivisions_7'] == "7_whitematter"]

labels_leftright = np.copy(labels)
labels_combinedlr = np.copy(labels)
#for i in range(int(np.shape(labels)[0]/2-2), int(np.shape(labels)[0])):
for i in range(np.shape(labels)[0]):
    for j in range(np.shape(labels)[1]):
        for k in range(np.shape(labels)[2]):
            if labels[i,j,k]>1000:
                labels_leftright[i, j, k] = labels[i, j, k] - 834
                labels_combinedlr[i, j, k] = labels[i, j, k] - 1000

save_nifti("/Users/alex/jacques/connectomes_testing/DWI_RAS/N57436_chass_symmetric3_labels_RAS_combined.nii.gz", labels_combinedlr, affine_labels)
save_nifti("/Users/alex/jacques/connectomes_testing/DWI_RAS_labelsbackup/N57436_chass_symmetric3_labels_RAS_lrordered.nii.gz", labels_leftright, affine_labels)

print("done")


import pandas as pd
import numpy as np
from dipy.io.image import load_nifti, save_nifti


atlas_legends = "/Users/alex/jacques/connectomes_testing/atlases/CHASSSYMM3AtlasLegends.xlsx"
df = pd.read_excel(atlas_legends, sheet_name='Sheet1')
df['Structure'] = df['Structure'].str.lower()
index1=df['index']
index2=df['index2']

labels, affine_labels = load_nifti("/Users/alex/jacques/connectomes_testing/DWI_RAS_labelsbackup/N57437_chass_symmetric3_labels_RAS.nii.gz")
rslt_whitem = df.loc[df['Subdivisions_7'] == "7_whitematter"]

labels_leftright = np.copy(labels)
labels_combinedlr = np.copy(labels)
for i in range(int(np.shape(labels)[0]/2-2), int(np.shape(labels)[0])):
    for j in range(np.shape(labels)[1]):
        for k in range(np.shape(labels)[2]):
            if labels[i,j,k]>1000:
                labels_leftright[i, j, k] = labels[i, j, k] - 834
                labels_combinedlr[i, j, k] = labels[i, j, k] - 1000

save_nifti("/Users/alex/jacques/connectomes_testing/DWI_RAS/N57437_chass_symmetric3_labels_RAS_combined.nii.gz", labels_combinedlr, affine_labels)
save_nifti("/Users/alex/jacques/connectomes_testing/DWI_RAS_labelsbackup/N57437_chass_symmetric3_labels_RAS_lrordered.nii.gz", labels_leftright, affine_labels)

print("done")


import pandas as pd
import numpy as np
from dipy.io.image import load_nifti, save_nifti


atlas_legends = "/Users/alex/jacques/connectomes_testing/atlases/CHASSSYMM3AtlasLegends.xlsx"
df = pd.read_excel(atlas_legends, sheet_name='Sheet1')
df['Structure'] = df['Structure'].str.lower()
index1=df['index']
index2=df['index2']

labels, affine_labels = load_nifti("/Users/alex/jacques/connectomes_testing/DWI_RAS_labelsbackup/N57440_chass_symmetric3_labels_RAS.nii.gz")
rslt_whitem = df.loc[df['Subdivisions_7'] == "7_whitematter"]

labels_leftright = np.copy(labels)
labels_combinedlr = np.copy(labels)
for i in range(int(np.shape(labels)[0]/2-2), int(np.shape(labels)[0])):
    for j in range(np.shape(labels)[1]):
        for k in range(np.shape(labels)[2]):
            if labels[i,j,k]>1000:
                labels_leftright[i, j, k] = labels[i, j, k] - 834
                labels_combinedlr[i, j, k] = labels[i, j, k] - 1000

save_nifti("/Users/alex/jacques/connectomes_testing/DWI_RAS/N57440_chass_symmetric3_labels_RAS_combined.nii.gz", labels_combinedlr, affine_labels)
save_nifti("/Users/alex/jacques/connectomes_testing/DWI_RAS_labelsbackup/N57440_chass_symmetric3_labels_RAS_lrordered.nii.gz", labels_leftright, affine_labels)

print("done")

"""