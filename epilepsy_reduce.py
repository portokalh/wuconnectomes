from tract_handler import target, prune_streamlines, get_trk_params, get_tract_params, gettrkpath, reducetractnumber, reducetractnumber_all
import os

folder = "/Volumes/Data/Badea/Lab/human/Sinha_epilepsy/TRK_basedondenoise/"
folder = "/mnt/munin6/Badea/Lab/human/Sinha_epilepsy/TRK_basedondenoise/"
for filename in os.listdir(folder):
    if filename.find("pruned"):
        newfile = filename.replace("pruned", "pruned_100")
        reducetractnumber(os.path.join(folder + filename), os.path.join(folder, newfile), getdata=True, ratio=100, return_affine= False, verbose=False)