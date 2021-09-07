import xlsxwriter
import pandas as pd

QCLab_path = "/Users/alex/Downloads/QCLAB_AD_mice062921.csv"

df = pd.read_csv(QCLab_path)
#df['Structure'] = df['Structure'].str.lower()


df = df[df.DWI != "Blank"]
df = df[df['DWI'].notna()]

genotypes = df.Genotype.unique()
treatments = df.Treatment.unique()
for genotype in genotypes:
    for treatment in treatments:
        number_gen = len(df[(df['Genotype'] == genotype) & (df['Treatment'] == treatment)])
        #if number_gen!=0:
        print(f"There are {number_gen} with genotype {genotype} and treatment {treatment}")


"""
genotype = df['index']
index2 = df['index2']
structures = df['Structure']
hemispheres = df['Hemisphere']
"""