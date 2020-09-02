
import xlsxwriter
import pandas as pd
import numpy as np

def round_array(array,rounder):

    newarray = np.zeros(np.shape(array))
    if np.size(np.shape(array)) == 2:
        for i in np.arange(np.shape(array)[0]):
            for j in np.arange(np.shape(array)[1]):
                newarray[i,j] = round(array[i,j],rounder)

    return(newarray)

def excel_average(files, file_outpath):

    workbook = xlsxwriter.Workbook(file_outpath)
    worksheet = workbook.add_worksheet()

    for file in files:

        df = pd.read_excel(file, sheet_name='Sheet1')
        data = df.values
        row = 0
        if 'datasum' not in locals():
            worksheet.write_row(0, 1, data[:166, 0])
            worksheet.write_column(1, 0, data[:166, 0])
            datasum = data[:, 1:]
        else:
            datasum = datasum + data[:, 1:]

    numfiles = np.size(files)
    dataavg = datasum / numfiles
    dataavg = round_array(dataavg, 1)
    for i in np.arange(np.shape(dataavg)[0]):
        worksheet.write_row(i+1, 1, dataavg[i,:])

    workbook.close()
    #for col, data in enumerate(connectome):
    #    worksheet.write(row + 1, col + 1, data)
