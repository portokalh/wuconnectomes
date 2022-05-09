
import xlsxwriter
import pandas as pd
import numpy as np
from os import remove as file_remove
from BIAC_tools import send_mail, isempty
from computer_nav import make_temppath, glob_remote, checkfile_exists_remote, remove_remote
"""
Created by Jacques Stout

Set of helper functions
Currently has:
make every float an int
Calculate the average in each cell of values from multiple excel files (found in same folder)
set up the excel files for connectomes and grouping
"""


def connectomes_to_excel(connectome, index_to_struct, output_path, overwrite=False, verbose=False, sftp=None):

    #df = pd.read_excel(ROI_excel, sheet_name='Sheet1')
    #structure = df['Structure']

    if checkfile_exists_remote(output_path, sftp=sftp):
        if overwrite:
            remove_remote(output_path, sftp=sftp)
        else:
            if verbose:
                print(f'Connectome already saved at {output_path}')
            return

    if sftp is None:
        output_tpath = output_path
    else:
        output_tpath = make_temppath(output_path)

    workbook = xlsxwriter.Workbook(output_tpath)
    worksheet = workbook.add_worksheet()

    for num in np.arange(1, np.shape(connectome)[0]+1):
        worksheet.write(0, num, index_to_struct[num])
        worksheet.write(num, 0, index_to_struct[num])

    row=0
    for col, data in enumerate(connectome):
        worksheet.write_column(row+1, col+1, data)

    workbook.close()

    if sftp is not None:
        sftp.put(output_tpath, output_path)
        file_remove(output_tpath)

    if verbose:
        print(f'Saved connectome at {output_path}')

    return


def grouping_to_excel(grouping, index_to_struct, output_path, overwrite=False, verbose=False, sftp=None):

    if checkfile_exists_remote(output_path):
        if overwrite:
            remove_remote(output_path, sftp)
        else:
            if verbose:
                print(f'Grouping already saved at {output_path}')
            return

    if sftp is None:
        output_tpath = output_path
    else:
        output_tpath = make_temppath(output_path)

    workbook = xlsxwriter.Workbook(output_tpath)
    worksheet = workbook.add_worksheet()

    for num in np.arange(1, np.shape(grouping)[0]):
        worksheet.write(0, num, index_to_struct[num])
        worksheet.write(num, 0, index_to_struct[num])

    row = 0
    for i in np.arange(np.shape(grouping)[0]):
        for j in np.arange(np.shape(grouping)[1]):
            worksheet.write(i + 1, j + 1, str(grouping[i, j]))

    workbook.close()

    if sftp is not None:
        sftp.put(output_tpath, output_path)
        file_remove(output_tpath)

    if verbose:
        print(f'Saved grouping at {output_path}')

    return

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

def M_grouping_excel_save(M,grouping,M_path, grouping_path, index_to_struct, verbose=False, sftp=None):

    matrix_sl = np.empty(np.shape(M), dtype=object)
    for i in np.arange(np.shape(matrix_sl)[0]):
        for j in np.arange(np.shape(matrix_sl)[1]):
            matrix_sl[i, j] = []
    for key in grouping.keys():
        matrix_sl[key] = grouping[key]
        matrix_sl[tuple(np.flip(key))] = grouping[key]

    M = np.delete(M, 0, 0)
    M = np.delete(M, 0, 1)

    matrix_sl = np.delete(matrix_sl, 0, 0)
    matrix_sl = np.delete(matrix_sl, 0, 1)

    if sftp is None:
        M_tpath = M_path
        grouping_tpath = grouping_path
    else:
        M_tpath = make_temppath(M_path)
        grouping_tpath = make_temppath(grouping_path)

    connectomes_to_excel(M, index_to_struct, M_tpath)
    grouping_to_excel(matrix_sl, index_to_struct, grouping_tpath)

    if sftp is not None:
        sftp.put(M_tpath, M_path)
        sftp.put(grouping_tpath, grouping_path)
        file_remove(M_tpath)
        file_remove(grouping_tpath)

    if verbose:
        txt = (f"The excelfile for connectome and grouping were saved at {M_path} and {grouping_path}")
        send_mail(txt, subject="Excel save")
        print(txt)


def extract_grouping(grouping_path, index_to_struct, shape=None, verbose=False, sftp=None):

    #grouping_array = np.empty(shape, dtype=object)
    grouping_newdic = {}
    if sftp is not None:
        grouping_paths = glob_remote(grouping_path, sftp)
        if np.size(grouping_paths == 1):
            grouping_path = grouping_paths[0]
            sftp.get(grouping_path, make_temppath(grouping_path))
            grouping_frame = pd.read_excel(make_temppath(grouping_path))
            file_remove(make_temppath(grouping_path))
        elif np.size(grouping_paths > 1):
            raise Exception('too many grouping files???')
        else:
            raise Exception('could not find grouping_path')
    else:
        grouping_frame = pd.read_excel(grouping_path)
    if shape is None:
        shape = list(grouping_frame.shape)
        shape[0] = shape[0]+1
        shape = tuple(shape)
    for i in np.arange(shape[0]-1):
        for j in np.arange(shape[1]-1):
            liststring = grouping_frame.iloc[i, j+1]
            liststring = liststring.replace('[','')
            liststring = liststring.replace(']','')
            liststring = (liststring.split(','))
            #print(liststring[-1])
            #liststring.pop(-1)
            if liststring[0] != '' and liststring[0] != ' ':
                try:
                    int(liststring[-1])
                except ValueError:
                    liststring = liststring[:-1]
                liststring = [int(i) for i in liststring]
            else:
                liststring = []
            grouping_newdic[i + 1, j + 1] = liststring
    return grouping_newdic

    #for key in grouping.keys():
    #    matrix_sl[key] = grouping[key]
    #    matrix_sl[tuple(np.flip(key))] = grouping[key]

