import xlsxwriter
import pandas as pd
import numpy as np

df_path = '/Users/alex/Downloads/18ABB11_readable101921_fixed.xlsx'
df_APOE = pd.read_excel(df_path)
#df_control = df_APOE.query('Age_Months < 16 and Genotype == APOE33')
df_control=df_APOE[(df_APOE["Genotype"] == 'APOE33') & (df_APOE['Diet'] == 'Control')]
df_else = df_APOE[(df_APOE["Genotype"] != 'APOE33') | (df_APOE['Diet'] != 'Control')]
print(df_control['DWI'])
dwis = []
for dwi in df_else['DWI']:
    dwis.append(dwi)

dwis_str = ','.join(dwis)
print(dwis_str, np.shape(dwis))