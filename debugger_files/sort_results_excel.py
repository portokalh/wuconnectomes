import pandas as pd
import shutil
import os
import math

myexcel = "/Volumes/Data/Badea/Lab/mouse/C57_JS/AMD_KEY/Whitson_AMD_groups_of_runnos.xlsx"
input = "/Volumes/Data/Badea/Lab/mouse/C57_JS/VBM_whiston_Figs/"
output = "/Volumes/Data/Badea/Lab/mouse/C57_JS/Whiston_reordered/"
df = pd.read_excel(myexcel, sheet_name='Sheet1')
#df['Structure'] = df['Structure'].str.lower()
#data_top = df.head()
for column in df.columns:
    for val in df[column]:
        if isinstance(val, str):
            print(str(val))
            oldfilename = input + val + "_stepsize_2_all_wholebrain_connectomes.xlsx"
            newfilename = output + column + "/" + val + "_connectomes.xlsx"
            if not os.path.exists(newfilename):
                shutil.copyfile(oldfilename, newfilename)
