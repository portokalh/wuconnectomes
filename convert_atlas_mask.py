
import pandas as pd
import numpy as np
from dipy.io.image import load_nifti, save_nifti
from computer_nav import load_nifti_remote, make_temppath

def atlas_converter(ROI_excel):

    df = pd.read_excel(ROI_excel, sheet_name='Sheet1')
    df['Structure'] = df['Structure'].str.lower()
    index1=df['index']
    index2=df['index2']
    index3=df['index3']
    structures=df['Structure']
    hemispheres=df['Hemisphere']
    hemispheres_new = []
    for i in np.arange(np.size(hemispheres)):
        if hemispheres[i] == "Left":
            hemispheres_new.append('_left')
        if hemispheres[i] == "Right":
            hemispheres_new.append('_right')
    converter_lr = {}
    converter_comb = {}
    index_to_struct_lr = {}
    index_to_struct_comb = {}
    for i in np.arange(np.size(hemispheres_new)):
        if hemispheres_new[i] in ['_left','_right']:
            converter_lr[index2[i]] = index1[i]
            converter_comb[index2[i]] = index3[i]
            index_to_struct_lr[index1[i]] = structures[i] + hemispheres_new[i]
            index_to_struct_comb[index3[i]] = structures[i]
    if 0 not in converter_lr:
        converter_lr[0] = 0
    if 0 not in converter_comb:
        converter_comb[0] = 0
    return converter_lr, converter_comb, index_to_struct_lr, index_to_struct_comb

def IIT_converter(ROI_excel):

    df = pd.read_excel(ROI_excel, sheet_name='Sheet1')
    df['Structure'] = df['Structure'].str.lower()
    index1=df['index']
    index2=df['index2']
    structures=df['Structure']
    hemispheres=df['Hemisphere']
    hemispheres_new = []
    for i in np.arange(np.size(hemispheres)-1):
        if hemispheres[i] == "Left":
            hemispheres_new.append('_left')
        if hemispheres[i] == "Right":
            hemispheres_new.append('_right')
    converter_lr = {}
    converter_comb = {}
    index_to_struct_lr = {}
    index_to_struct_comb = {}
    for i in np.arange(np.size(index1)-1):
        converter_lr[index2[i]] = index1[i]
        converter_comb[index2[i]] = index3[i]
        index_to_struct_lr[index1[i]] = structures[i] + hemispheres_new[i]
        index_to_struct_comb[index3[i]] = structures[i]
    return converter_lr, converter_comb, index_to_struct_lr, index_to_struct_comb

def convert_labelmask(atlas, converter, atlas_outpath = None, affine_labels=None, sftp=None):

    if isinstance(atlas, str):
        labels, affine_labels, _, _, _ = load_nifti_remote(atlas, sftp)
    else:
        if affine_labels is None:
            raise TypeError('Need to add the affine labels if directly including label array')
        else:
            labels = atlas

    labels_new = np.copy(labels)

    for i in range(np.shape(labels)[0]):
        for j in range(np.shape(labels)[1]):
            for k in range(np.shape(labels)[2]):
                try:
                    labels_new[i, j, k] = converter[labels[i, j, k]]
                except:
                    print('hi')

    if sftp is not None:
        save_nifti(make_temppath(atlas_outpath), labels_new, affine_labels)
        sftp.put(make_temppath(atlas_outpath), atlas_outpath)
    else:
        save_nifti(atlas_outpath, labels_new, affine_labels)

    return(labels_new)

def run_onall():

    atlas_legends = "/Users/alex/jacques/connectomes_testing/atlases/CHASSSYMM3AtlasLegends.xlsx"
    df = pd.read_excel(atlas_legends, sheet_name='Sheet1')
    df['Structure'] = df['Structure'].str.lower()
    index1=df['index']
    index2=df['index2']

    l = ['N57442', 'N57446', 'N57447','N57449','N57451','N57496','N57498','N57500','N57502','N57504','N57513','N57515','N57518', 'N57520','N57522','N57546','N57447','N57449','N57451','N57496','N57498','N57500','N57502','N57504','N57513','N57515','N57518','N57520','N57522','N57546','N57548', 'N57550', 'N57552', 'N57554', 'N57559', 'N57580', 'N57582', 'N57584', 'N57587', 'N57590', 'N57692', 'N57694', 'N57700', 'N57702', 'N57709']

    converter_lr, converter_comb = chassym3_converter(atlas_legends)

    atlas_folder = "/Volumes/Data/Badea/Lab/mouse/C57_JS/DWI_RAS_40subj/"

    for subject in l:
        labels, affine_labels = load_nifti(atlas_folder + subject + "_chass_symmetric3_labels_RAS.nii.gz")
        rslt_whitem = df.loc[df['Subdivisions_7'] == "7_whitematter"]

        labels_leftright = np.copy(labels)
        labels_combinedlr = np.copy(labels)
        # for i in range(int(np.shape(labels)[0]/2-2), int(np.shape(labels)[0])):
        for i in range(np.shape(labels)[0]):
            for j in range(np.shape(labels)[1]):
                for k in range(np.shape(labels)[2]):
                    if labels[i, j, k] > 1000:
                        labels_leftright[i, j, k] = converter_lr[labels[i,j,k]]
                        labels_combinedlr[i, j, k] = converter_comb[labels[i,j,k]]

        save_nifti(atlas_folder + subject + "_chass_symmetric3_labels_RAS_combined.nii.gz", labels_combinedlr,
                   affine_labels)
        save_nifti(atlas_folder + subject + "_chass_symmetric3_labels_RAS_lrordered.nii.gz", labels_leftright,
                   affine_labels)

        print("done")

    """
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
