
from tract_manager import reducetractnumber

#oldtrkfile='/Volumes/dusom_dibs_ad_decode/all_staff/munin3badea/mouse/C57_JS/TRK_RAS_40subj/N57449_wholebrain_all_stepsize_2.trk'
#newfile = '/Volumes/dusom_dibs_ad_decode/all_staff/munin3badea/mouse/C57_JS/TRK_RAS_40subj/N57449_wholebrain_little_stepsize_2.trk'
oldtrkfile = '/Volumes/dusom_dibs_ad_decode/all_staff/VBM_whiston_QA/H29056_stepsize_2_ratio_100_wholebrain_pruned.trk'
newfile = '/Volumes/dusom_dibs_ad_decode/all_staff/VBM_whiston_QA/H29056_stepsize_2_ratio_1000_wholebrain_pruned.trk'

reducetractnumber(oldtrkfile, newfile, False, 10, False)