

from tract_manager import connectomes_to_excel
import pandas as pd
import pickle
import xlsxwriter

connectome_path = "/Users/alex/jacques/connectomes_testing/Fig_RAS_connectomes/N57433vsmall_grouping.p"
ROI_excel = "/Users/alex/jacques/connectomes_testing//atlases/CHASSSYMM3AtlasLegends.xlsx"
output_path = "/Volumes/Data/Badea/Lab/mouse/C57_JS/Figures_RAS_diff/N57433_vsmall_connectomes.xlsx"

df = pd.read_excel(ROI_excel, sheet_name='Sheet1')
structure = df['Structure']
connectome = pickle.load(open(connectome_path, "rb"))

workbook = xlsxwriter.Workbook(output_path)

worksheet = workbook.add_worksheet()

num = 1
for struct in structure:
    worksheet.write(0, num, struct)
    worksheet.write(num, 0, struct)
    num += 1

row=0
for col, data in enumerate(connectome):
    worksheet.write_column(row+1, col+1, data)

workbook.close()
