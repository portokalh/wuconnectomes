
import os
import sys

attribute_file = sys.argv[1]

if not os.path.exists(attribute_file):
    print("The file could not be found at "+attribute_file+", end of scrypt run")

with open(attribute_file, 'rb') as source:
    if verbose: print('INFO    : Extracting pipeline parameters')
    header_size = source.read(4)
    header_size = struct.unpack('I', header_size)
    if verbose: print('INFO    : Header size = ', int(header_size[0]))
    i = 0
    stopsign = 200
    for line in source:

        # pattern1 = '<ParamLong."BaseResolution">  {*'
        # rx1 = re.compile(pattern1, re.IGNORECASE|re.MULTILINE|re.DOTALL)

        # pattern1 = 'z_Agilent_bvalue_m00='
        pattern1 = 'Diffusion folder'
        rx1 = re.compile(pattern1, re.IGNORECASE | re.MULTILINE | re.DOTALL)
        pattern2 = 'Tract folder'
        rx2 = re.compile(pattern2, re.IGNORECASE | re.MULTILINE | re.DOTALL)
        pattern3 = 'Figures folder'
        rx3 = re.compile(pattern3, re.IGNORECASE | re.MULTILINE | re.DOTALL)
        pattern4 = 'Atlas file'
        rx4 = re.compile(pattern4, re.IGNORECASE | re.MULTILINE | re.DOTALL)
        pattern5 = 'Subject list'
        rx5 = re.compile(pattern5, re.IGNORECASE | re.MULTILINE | re.DOTALL)

        i += 1
        if i == stopsign:
            print("hi")
        for a in rx1.findall(str(line)):
            bvals_all = str(line).split(',')[1]
            bvals = bvals_all.split('\\')[0]
        for a in rx2.findall(str(line)):
            dsl_all = str(line).split(',')[1]
            dsl = dsl_all.split('\\')[0]
            # dsl=([float(s) for s in dsl_all.split() if s.isnumeric()])
        for a in rx3.findall(str(line)):
            dpe_all = str(line).split(',')[1]
            dpe = dpe_all.split('\\')[0]
        for a in rx4.findall(str(line)):
            dro_all = str(line).split(',')[1]
            dro = dro_all.split('\\')[0]