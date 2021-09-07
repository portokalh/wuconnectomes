import xlsxwriter
import pandas as pd
import numpy as np

QCLab_path = "/Users/alex/Downloads/QCLAB_AD_mice062921_clean.csv"
Concatweight_path = "/Users/alex/Downloads/MasterSheet_Experiments2021-4_9_3_21.xlsx"
df_QC = pd.read_csv(QCLab_path)
df_ConWeight = pd.read_excel(Concatweight_path)
#df['Structure'] = df['Structure'].str.lower()

#rd stands for human readable for humans, R is the one used as R input. Both are built here from same principles
all_csv_rd = "/Users/alex/Downloads/Mouse_experiments_readable.csv"
all_csv_R = "/Users/alex/Downloads/Mouse_experiments_R.csv"

#df_ConWeight = df_ConWeight[df_ConWeight.Success == True]
glucose_weeks = [0, 6, 12, 20, 30]
glucose_weeks_list=[]
glucose_weeks_dir = {}
for glucose_week in glucose_weeks:
    oldformat = f'Glucose_week{glucose_week}'
    newformat = f'glucose_w{glucose_week}'
    glucose_weeks_dir[oldformat]=newformat
    glucose_weeks_list.append(newformat)

rename_cols_dir={'Animal_ID': 'animalid','Genotype':'genotype', 'Sex':'sex'}
#df_ConWeight.rename(columns=rename_cols_dir, inplace=True)

rename_cols_dir.update(glucose_weeks_dir)
weight_weeks_list=[]
weight_weeks_dir = {}

number_weeks_weight = 28
for i in np.arange(number_weeks_weight+1):
    oldformat = f'Weight_week{i}'
    newformat = f'weight_w{i}'
    weight_weeks_dir[oldformat]=newformat
    weight_weeks_list.append(newformat)

if number_weeks_weight>np.max(glucose_weeks):
    number_weeks = np.max(glucose_weeks)
else:
    number_weeks = np.max(glucose_weeks)

rename_cols_dir.update(weight_weeks_dir)

df_ConWeight.rename(columns=rename_cols_dir, inplace=True)

allcolumns = ["animalid", "genotype","sex"] + glucose_weeks_list +weight_weeks_list
df_allcsv_human = df_ConWeight[allcolumns]

df_QC = df_QC[df_QC.DWI != "Blank"]
df_QC = df_QC[df_QC['DWI'].notna()]

#df_allcsv_human['DWI'] = ''
#df_allcsv_human['GRE'] = ''

#df_QC_2 = df_QC.copy()

#df_QC_2.replace(to_replace=r':0$', value='', regex=True)
#df_QC_2.replace(to_replace=r':1$', value='', regex=True)

#Duplicates are 191205-1, 191205-2, 191205-3, 191205-4
df_QC = df_QC.drop_duplicates('Animal')
#duplicate = df_QC_2[df_QC_2.duplicated()]

df_QC = df_QC[['Animal', 'DWI', 'GRE']]
df_QC.rename(columns={'Animal':'animalid'}, inplace=True)

df_allcsv_human = df_allcsv_human.join(df_QC.set_index('animalid'), on='animalid')

df_allcsv_human.to_csv(all_csv_rd)

df_allcsv_human = df_allcsv_human[df_allcsv_human['animalid'].notna()]
if not pd.Series(df_allcsv_human["animalid"]).is_unique:
    raise Exception("The subject names in this dataframe are not unique, cannot proceed with next step for R")

subject_names = df_allcsv_human['animalid'].unique().tolist()
R_columns = [['week','animalid','genotype','sex','weight','glucose','DWI',"GRE"]]
df_allcsv_R = pd.DataFrame(columns = list(R_columns))

index=0
for week in np.arange(number_weeks + 1):
    for subject_name in subject_names :
        if set([f'weight_w{week}']).issubset(df_allcsv_human.columns):
            try:
                weight_val = df_allcsv_human[df_allcsv_human['animalid'] == subject_name][f'weight_w{week}'].values[0]
            except:
                raise Exception('this is for debugging, but clearly it crashed xD')
        else:
            weight_val = np.nan
        if set([f'glucose_w{week}']).issubset(df_allcsv_human.columns):
            glucose_val = df_allcsv_human[df_allcsv_human['animalid'] == subject_name][f'glucose_w{week}'].values[0]
        else:
            glucose_val = np.nan
        try:
            df_allcsv_R.loc[index]= [week,subject_name,df_allcsv_human[df_allcsv_human['animalid']==subject_name]['genotype'].values[0],
                        df_allcsv_human[df_allcsv_human['animalid'] == subject_name]['sex'].values[0], weight_val, glucose_val,
                        df_allcsv_human[df_allcsv_human['animalid'] == subject_name]['DWI'].values[0],
                        df_allcsv_human[df_allcsv_human['animalid'] == subject_name]['GRE'].values[0]]
        except:
            raise Exception('this is for debugging, but clearly it crashed xD')
        index += 1
df_allcsv_R.to_csv(all_csv_R)

"""
#df_dict = {name: df_allcsv_human.loc[df_allcsv_human['customer name'] == name] for name in customerNames}

#df_dict['Name3']

animals = df_QC.AnimalID.unique()
for animal in animals:
    newrow = df_ConWeight['animal-id'] == animal
    newrow =
    mydataframe = mydataframe.append(new_row, ignore_index=True)
    
glucose_weeks_dir = {'Glucose_week0':'glucose_week0', 'Glucose_week6':'glucose_week6',
                     'Glucose_week12':'glucose_week12', 'Glucose_week20':'glucose_week20',
                     'Glucose_Week30':'glucose_week30'}
"""
"""
df_QC = df_QC[df_QC.DWI != "Blank"]
df_QC = df_QC[df_QC['DWI'].notna()]

genotypes = df_QC.Genotype.unique()
treatments = df_QC.Treatment.unique()
for genotype in genotypes:
    for treatment in treatments:
        number_gen = len(df_QC[(df_QC['Genotype'] == genotype) & (df_QC['Treatment'] == treatment)])
        #if number_gen!=0:
        print(f"There are {number_gen} with genotype {genotype} and treatment {treatment}")
"""