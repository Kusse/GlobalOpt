from pandas import DataFrame
import pandas as pd
import numpy as np
import xlwt
import re
myCNAfile="CNanalysis.log"
ovito = {
'200':[(2, 0, 0), (0.06, 0.29, 0.17, 0.12, 0.13, 0.09, 0.09, 0.04, 0.04, 0.04, 0.04),[]],
'211':[(2, 1, 1), (0.06, 0.08, 0.05, 0.05, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04),[]],
'300':[(3, 0, 0), (0.01, 0.01, 0.02, 0.02, 0.02, 0.01, 0.01, 0.00, 0.00, 0.00, 0.00),[]],
'311':[(3, 1, 1), (0.26, 0.22, 0.26, 0.27, 0.27, 0.26, 0.26, 0.23, 0.23, 0.21, 0.21),[]],
'322':[(3, 2, 2), (0.18, 0.17, 0.21, 0.19, 0.20, 0.24, 0.24, 0.24, 0.24, 0.25, 0.25),[]],
'421':[(4, 2, 1), (0.07, 0.00, 0.01, 0.03, 0.03, 0.03, 0.03, 0.04, 0.04, 0.03, 0.03),[]],
'422':[(4, 2, 2), (0.25, 0.10, 0.17, 0.19, 0.20, 0.20, 0.20, 0.19, 0.19, 0.16, 0.16),[]],
'433':[(4, 3, 3), (0.03, 0.03, 0.03, 0.03, 0.03, 0.05, 0.05, 0.08, 0.08, 0.08, 0.08),[]],
'544':[(5, 4, 4), (0.04, 0.01, 0.01, 0.03, 0.02, 0.04, 0.04, 0.08, 0.08, 0.10, 0.10),[]],
'555':[(5, 5, 5), (0.05, 0.03, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.05, 0.05),[]]
}
SGS = {}
RC = []
COL = []
counter=0
f=open(myCNAfile, "r")
white_spaces = [""," ","\t", "\t\t","\n","\n\n"]
while True:
    CNApairs = []
    line = f.readline()
    counter = counter + 1
    line = line.rstrip()
    line = re.split(r'\t+', line)
    if(line[0] not in white_spaces and line[0].split()[0]=="CutOff"):
        rc = round(float(line[-1]),2) 
        RC.append(rc)
    if (len(line)==4 and line[0]=="" and len(line[1])==9):
        print("line = ", line)
        pair = line[1] 			#;  print("pair = ", pair)
        relative_abundance = line[-1] 	#; print("rabd = ", relative_abundance)
        CNApairs.append([pair, relative_abundance])
        table_end = "*****************************************************"
        while line != table_end:
            line = f.readline()
            counter = counter + 1
            line = line.rstrip()
            line = re.split(r'\t+', line)
            if(len(line) != 4):
                break
            pair = line[1]; print("pair = ", pair)
            relative_abundance = line[-1] ; print("rabd = ", relative_abundance)
            CNApairs.append([pair, relative_abundance])
        SGS[rc] = CNApairs
        if(len(RC)>=10):
            break
    print("line counter = ", len(RC), counter, RC)
f.close()
for key in ovito.keys():
    print(key," : ", ovito[key][2])
my_table = {}
for key in SGS.keys():
    TEMP = SGS[key]
    values = []
    pairs = []
    for TMP in TEMP:
        pairs.append(TMP[0])
        values.append(np.round(float(TMP[1])*0.01,2))
    my_table[key] = pd.Series(np.array(values), index=pairs)
df = pd.DataFrame(my_table)
print(df)
writer = pd.ExcelWriter('Au55_CNA_SGS.xlsx', engine='xlsxwriter')
df.to_excel(writer, sheet_name='Sheet1')
writer.save()
