
#Just something to extract the pickle information that was created and convert it to an excel. Useful
#if there was a major problem with the generation of connectomes

import pickle
from tract_manager import connectomes_to_excel

BIGGUS_DISKUS = "/Volumes/dusom_dibs_ad_decode/all_staff/munin3badea/mouse"

outpath = BIGGUS_DISKUS + "/C57_JS/Figures_RAS_40subj/"

l = ['N57442', 'N57446', 'N57447','N57449','N57451','N57496','N57498','N57500','N57502','N57504','N57513','N57515','N57518', 'N57520','N57522','N57546','N57447','N57449','N57451','N57496','N57498','N57500','N57502','N57504','N57513','N57515','N57518','N57520','N57522','N57546','N57548', 'N57550', 'N57552', 'N57554', 'N57559', 'N57580', 'N57582', 'N57584', 'N57587', 'N57590', 'N57692', 'N57694', 'N57700', 'N57702', 'N57709']
tractsize="all"

atlas_legends = BIGGUS_DISKUS + "/atlases/CHASSSYMM3AtlasLegends.xlsx"
ROI_excel = "/Users/alex/jacques/connectomes_testing//atlases/CHASSSYMM3AtlasLegends.xlsx"

for subject in l:
    excel_path = outpath + subject + "_" + tractsize + "_connectomes_v3ssh.xlsx"
    picklepath_connect = outpath + subject + "_" + tractsize + '_connectomes_v2.p'
    M = pickle.load(open(picklepath_connect, "rb"))
    connectomes_to_excel(M, ROI_excel, excel_path)
    print('did subject '+subject)
