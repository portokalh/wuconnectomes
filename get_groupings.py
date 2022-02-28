import pandas as pd

file_project = '/Users/alex/AD_Decode_hydra/ADDECODE2021/AD_DECODE_data_020722_JS.xlsx'

dataframe = pd.read_excel(file_project)
listM = dataframe.loc[dataframe['sex'] == 'M']['MRI_Exam']
listF = dataframe.loc[dataframe['sex'] == 'F']['MRI_Exam']

listM_str = []
listF_str = []

for subj in listM:
    listM_str.append('S0'+str(subj))
for subj in listF:
    listF_str.append('S0'+str(subj))

print(listM_str)

print(listF_str)