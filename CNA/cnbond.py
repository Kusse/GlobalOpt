
#************ Bond Indices Of EACH Atom ************
#  3 Atom:      1 (CNtag=6)
#  4         Nbor:     20    (211) Pair
#  5         Nbor:      6    (322) Pair
#  6         Nbor:     53    (322) Pair
#  7         Nbor:     27    (543) Pair
#  8         Nbor:     38    (211) Pair
#  9         Nbor:     49    (322) Pair
# 10 ***


import numpy as np
import pandas as pd
import collections
bondsfile = "CNbond.log"
natoms = 55
f=open(bondsfile, "r")
f.readline()
f.readline()
cnaPairs = {}
for i in range(natoms):
	line=f.readline().rstrip()
	line=line.split()
	if line[0] == "Atom:":
		nnbrs = int(line[2][7:-1])
		print("nnbrs = ", nnbrs)
		for j in range(nnbrs):
			line=f.readline().rstrip()
			line=line.split()
			pair = line[2][1:4]
			#print("{}: key = ", pair, "pair=", int(pair))
			if i==0 and j==0:
				key = int(pair)
				cnaPairs[key] = 1 
				continue
			isunique = True
			for key in cnaPairs.keys():
				print("key = ", pair,  key==int(pair), cnaPairs[key])
				if key == int(pair):
					isunique = False
					break
			if isunique:
				key=int(pair)
				cnaPairs[key] = 1 
			else:
				key=int(pair)
				cnaPairs[key] = cnaPairs[key] + 1	

			#print("{} : 	{}	{}".format(j+1, int(pair), cnaPairs[int(pair)]))
		
		f.readline().rstrip()
		f.readline().rstrip()
print(cnaPairs)
ordered_cnaPairs=collections.OrderedDict(sorted(cnaPairs.items()))
print(ordered_cnaPairs)
#table = pd.Series(np.array(ordered_cnaPairs.values()), index=ordered_cnaPairs.keys())
f.close()
pairs = []
values = []
for key, value in ordered_cnaPairs.items():
# 59     TEMP = SGS[key]
# 60     values = []
# 61     pairs = []
# 62     for TMP in TEMP:
         pairs.append(key)
         values.append(np.round(float(value),2))
table = pd.Series(np.round(np.array(values)/np.sum(values),2), index=pairs)
df = pd.DataFrame(table)
print(df)
# 68 writer = pd.ExcelWriter('Au55_CNA_SGS.xlsx', engine='xlsxwriter')
# 69 df.to_excel(writer, sheet_name='Sheet1')
# 70 writer.save()
#line = line.split()
#print(line) # = line.rstrip()
#print(line[0], line[1], line[2])
#print(line.split())
#print(line[0], line[1], line[2])
