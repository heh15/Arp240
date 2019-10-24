# import pandas as pd
import os
import sys
sys.path.insert(0,os.path.abspath('.'))
import numpy as np

# data=pd.read_csv('results_latex.tex', header=None, sep='&',skiprows=6, skipfooter=19)
text=[]
data=[]

filename='results_latex.tex'
with open(filename,'r') as input:
    for line in input.readlines():
        text.append(line.strip('\\\\\n'))
text_value=text[4:-2]

for line in text_value:
    list_tmp=line.split('&')[1:]
    list_tmp2=[float(i) for i in list_tmp]
    data.append(list_tmp2)

# Row
# 0 density
# 1 temperature
# 2 column density
# 3 filling factor
# 6 pressure
# 7 beam averaged column density

# Columns
# 4 Mode Mean
# 5 Mode Sigma

data=np.array(data)
data_table=np.copy(data[[0,1,2,3,6,7]][:, [4,5]])
error=np.copy(data_table[:,1])
data_table=np.column_stack((data_table,error))
logT=data_table[1][0]; logT_error=data_table[1][1]
data_table[1][1]=10**(logT+logT_error)-10**logT
data_table[1][2]=10**logT-10**(logT-logT_error)
data_table[1][0]=10**logT

logbf=data_table[3][0]; logbf_error=data_table[3][1]
data_table[3][1]=10**(logbf+logbf_error)-10**logbf
data_table[3][2]=10**logbf-10**(logbf-logbf_error)
data_table[3][0]=10**logbf
print(data_table)


# data_output=np.array2string(data_table, precision=3, separator='\t')
np.savetxt('results4table.txt', data_table, fmt='%3f', delimiter='\t', newline='\n')
