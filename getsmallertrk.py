
from tract_manager import reducetractnumber

oldtrkfile='/Volumes/dusom_dibs_ad_decode/all_staff/munin3badea/mouse/C57_JS/TRK_RAS_40subj/N57449_wholebrain_all_stepsize_2.trk'
newfile = '/Volumes/dusom_dibs_ad_decode/all_staff/munin3badea/mouse/C57_JS/TRK_RAS_40subj/N57449_wholebrain_little_stepsize_2.trk'

reducetractnumber(oldtrkfile, newfile, False, 10, True)