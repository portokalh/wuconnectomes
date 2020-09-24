
from tract_manager import gettrkpath, tract_getroi, getdwidata
from dipy.io.utils import create_tractogram_header
import pandas as pd
import numpy as np
from dipy.io.streamline import load_trk



trkpath = "/Volumes/Data/Badea/Lab/mouse/C57_JS/TRK_RAS_40subj/N57442_small_fimbria_stepsize_2_pruned.trk"

trkpath = "/mnt/munin6/Badea/Lab/mouse/C57_JS/TRK_RAS_40subj"
dwipath = "/mnt/munin6/Badea/Lab/mouse/C57_JS/DWI_RAS_40subj/"
str_identifier = "_small_fimbria_stepsize_2_pruned"
l = ['N57496','N57554','N57580','N57692','N57437']
verbose = True


atlas_legends = "/Users/alex/jacques/connectomes_testing//atlases/CHASSSYMM3AtlasLegends.xlsx"
atlas_legends = "/mnt/munin6/Badea/Lab/atlases/CHASSSYMM3AtlasLegends.xlsx"
df = pd.read_excel(atlas_legends, sheet_name='Sheet1')
df['Structure'] = df['Structure'].str.lower()

targetrois = ["hippocampus"]
labelslist = []



for targetroi in targetrois:

    rslt_df = df.loc[df['Structure'] == targetroi.lower()]
    if targetroi.lower() == "wholebrain" or targetroi.lower() == "brain":
        labelslist = None
    else:
        labelslist = np.concatenate((labelslist, np.array(rslt_df.index2)))

    for subject in l:
        trkfilepath = gettrkpath(trkpath, subject, str_identifier, verbose)
        trkdata = load_trk(trkfilepath, "same")
        if verbose:
            print("loaded ")
        trkdata.to_vox()
        if hasattr(trkdata, 'space_attribute'):
            header = trkdata.space_attribute
        elif hasattr(trkdata, 'space_attributes'):
            header = trkdata.space_attributes
        affine = trkdata._affine
        trkstreamlines = trkdata.streamlines

        trkroipath = trkpath + '/' + subject + str_identifier + '_' + targetroi + '.trk'
        myheader = create_tractogram_header(trkroipath, *header)
        _, _, _, labelmask, _, _, _, _ = getdwidata(dwipath, subject)
        tract_getroi(trkstreamlines, affine, myheader, labelslist, labelmask, trkroipath), verbose

    labelslist = []

