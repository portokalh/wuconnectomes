from nifti_handler import getfa, getdiffdata_all, getdiffdata, getdiffpath, getgtab, getlabelmask, move_bvals, getmask, getb0s
import numpy as np
from convert_atlas_mask import convert_labelmask, atlas_converter
import os

subjects = ["N58214", "N58215", "N58216", "N58217", "N58218", "N58219", "N58221", "N58222", "N58223", "N58224",
        "N58225", "N58226", "N58228",
        "N58229", "N58230", "N58231", "N58232", "N58633", "N58634", "N58635", "N58636", "N58649", "N58650",
        "N58651", "N58653", "N58654",
        'N58408', 'N58398', 'N58714', 'N58740', 'N58477', 'N58734', 'N58309', 'N58792', 'N58302',
        'N58784', 'N58706', 'N58361', 'N58355', 'N58712', 'N58790', 'N58606', 'N58350', 'N58608',
        'N58779', 'N58500', 'N58604', 'N58749', 'N58510', 'N58394', 'N58346', 'N58344', 'N58788', 'N58305',
        'N58514', 'N58794', 'N58733', 'N58655', 'N58735', 'N58310', 'N58400', 'N58708', 'N58780', 'N58512',
        'N58747', 'N58303', 'N58404', 'N58751', 'N58611', 'N58745', 'N58406', 'N58359', 'N58742', 'N58396',
        'N58613', 'N58732', 'N58516', 'N58402']
subjects = ['N57437', 'N57442', 'N57446', 'N57447', 'N57449', 'N57451', 'N57496', 'N57498', 'N57500', 'N57502', 'N57504', 'N57513', 'N57515', 'N57518', 'N57520', 'N57522', 'N57546', 'N57548', 'N57550', 'N57552', 'N57554', 'N57559', 'N57580', 'N57582', 'N57584', 'N57587', 'N57590', 'N57692', 'N57694', 'N57700', 'N57702', 'N57709']
#subjects = ["N59120"]
diffpath = '/mnt/paros_MRI/jacques/APOE/DWI_allsubj_RAS'
#diffpath = '/Users/jas/jacques/APOE_testing/labels'
verbose = True
ROI_excel = "/mnt/paros_MRI/jacques/atlases/CHASSSYMM3AtlasLegends.xlsx"
#ROI_excel = '/Volumes/Data/Badea/Lab/atlases/CHASSSYMM3AtlasLegends.xlsx'
overwrite=True

for subject in subjects:
	labelmask, labelaffine, labelpath = getlabelmask(diffpath, subject, verbose)
	if np.size(np.shape(labelmask)) == 1:
	    labelmask = labelmask[0]
	if np.size(np.shape(labelmask)) == 4:
	    labelmask = labelmask[:, :, :, 0]
	print("Mask shape is " + str(np.shape(labelmask)))

	converter_lr, converter_comb, index_to_struct_lr, index_to_struct_comb = atlas_converter(ROI_excel)

	labeloutpath = labelpath.replace('_labels','_labels_lr_ordered')
	if not os.path.isfile(labeloutpath) or overwrite:
	    labelmask = convert_labelmask(labelmask, converter_lr, atlas_outpath=labeloutpath,
		                  affine_labels=labelaffine)
