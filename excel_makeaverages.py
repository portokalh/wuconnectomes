
from excel_management import excel_average

lAPOE22 = ["N57496", "N57500", "N57694", "N57702", "N57709", "N57498", "N57502", "N57513", "N57552", "N57692", "N57700"]
lAPOE33 = ["N57582", "N57584", "N57590", "N57546", "N57548", "N57550", "N57580", "N57587"]
lAPOE44 = ["N57442", "N57447", "N57451", "N57518", "N57522", "N57437", "N57446", "N57449", "N57515", "N57520"]
HN = ["N57554", "N57559"]

base_path = "/Users/alex/jacques/connectomes_testing/connectomes_v2/"
out_folder = base_path

filepaths=[]
for subject in HN:

    filepath = base_path + subject + "_all_connectomes.xlsx"
    filepaths.append(filepath)
file_outpath = out_folder + "HN_connectomes_average.xlsx"
excel_average(filepaths,file_outpath)