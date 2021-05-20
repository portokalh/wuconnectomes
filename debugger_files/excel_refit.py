
l = ['N57437', 'N57442', 'N57446', 'N57447','N57449','N57451','N57496','N57498','N57500','N57502','N57504','N57513','N57515','N57518','N57520','N57522','N57546','N57548', 'N57550', 'N57552', 'N57554', 'N57559', 'N57580', 'N57582', 'N57584', 'N57587', 'N57590', 'N57692', 'N57694', 'N57700', 'N57702', 'N57709']

l = [ 'N57498']

base_path = "/Users/alex/jacques/connectomes_testing/connectomes_all/"
out_path = "/Users/alex/jacques/connectomes_testing/connectomes_v2/"

import xlsxwriter
import numpy as np
import pandas as pd

filepaths=[]
for subject in l:

    print("Begin process of subject " + subject)

    filepath = base_path + subject + "_all_connectomes.xlsx"
    fileoutpath = out_path + subject + "_all_connectomes.xlsx"

    df = pd.read_excel(filepath, sheet_name='Sheet1')
    data = df.values
    workbook = xlsxwriter.Workbook(fileoutpath)
    worksheet = workbook.add_worksheet()

    worksheet.write_row(0,1, data[:166,0])
    for i in np.arange(166):
        worksheet.write_row(i+1, 0, data[i, :167])

    workbook.close()

    print("end of process of subject " + subject)
