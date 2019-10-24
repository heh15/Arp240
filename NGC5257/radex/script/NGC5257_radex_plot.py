'''
Oct 4th, 2018

Read the data from csv file instead of text file. 
read_ratio: infile is csv file from ratio directory. 
read_value: filename is the result from radex file. 
get_file: given the column density and abundance ratio, get the filename for read_value. 

'''


import numpy as np
import pandas as pd
from shutil import copy
import os
from matplotlib.ticker import MaxNLocator
from matplotlib import rc
from matplotlib import rcParams
rcParams['mathtext.default']='regular'
import copy

def read_ratio(infile):
    ratio = []
    uncertainty=[]
    with open(infile,'r') as input:
        line=input.readline()
        line=input.readline()
        words=line.split()
        while (words != []):
            ratio.append(float(words[-2]))
            uncertainty.append(float(words[-1]))
            line=input.readline()
            words=line.split()
    return ratio,uncertainty

def read_value(filename,regions,ratio_table,mol):
    data=np.transpose(np.loadtxt(filename))
    temp_tmp=data[0]
    dens_tmp=data[1]
    ratio=data[2]
    if mol=='2110':
        ratio=np.power(ratio,-1)
    temp=dict.fromkeys(regions)
    dens=dict.fromkeys(regions)
    index=dict.fromkeys(regions)
    for region in regions:
        index[region]=np.where(np.abs(ratio-ratio_table[mol+'_ratio'][region]) < ratio_table[mol+'_uncertainty'][region])
        temp[region]=temp_tmp[index[region]]
        dens[region]=dens_tmp[index[region]]
    return index, temp, dens


def get_file(col,abu):
    co_12='co_1-0_2-1'
    co_13='12co_13co'
    co_12=co_12+'_'+col
    co_13=co_13+'_'+col
    co_13=co_13.replace(col,abu+'_'+col)
    co_12=co_12+'.dat'
    co_13=co_13+'.dat'
    picturename='radex_'+col+'_'+abu+'.png'
    return co_13,co_12,picturename


#parameters
col='166'
abu='50'

# basic setting
Dir='/home/heh15/workingspace/Arp240/radex/'
picDir=Dir+'picture/'
logDir=Dir+'log/'
scriptDir=Dir+'script/radex_run/'
co_12=get_file(col,abu)[0]
co_13=get_file(col,abu)[1]
picturename=picDir+'NGC5257_cal_'+get_file(col,abu)[2]

os.chdir(logDir)


# read the data from the ratio data file.
file_1213='../../ratio/log/NGC5257_1213_ratio_mod_cal.csv'
file_2110='../../ratio/log/NGC5257_2110_ratio_mod_cal.csv'

data_1213=pd.read_csv(file_1213,header=0,index_col=0)
data_1213=data_1213.drop(index='nonarm')

regions=list(data_1213.index)
values=['1213_ratio','1213_uncertainty','2110_ratio','2110_uncertainty']
variables=['temp','dens']

ratio_table=pd.DataFrame()
ratio_table[0]=data_1213['ratio']
ratio_table[1]=data_1213['uncertainty']

# read the data from file_2110.
data_2110=pd.read_csv(file_2110,header=0,index_col=0)
ratio=data_2110['ratio']['total']
uncertainty=data_2110['uncertainty']['total']
ratio_table[2]=ratio
ratio_table[3]=uncertainty

# add the header and index. 
ratio_table.index=regions
ratio_table.columns=values


os.chdir(scriptDir)


############################################################
# load data from radex file

filename=get_file(col,abu)[0]
mol='1213'
index_1213,temp_1213,dens_1213 = read_value(filename,regions,ratio_table,mol)

temp=np.transpose(np.loadtxt(filename))[0]
dens=np.transpose(np.loadtxt(filename))[1]

filename=get_file(col,abu)[1]
mol='2110'
index_2110,temp_2110,dens_2110 = read_value(filename,regions,ratio_table,mol)



############################################################
# find the cross region in parameter space. 

index=dict.fromkeys(regions)
temp_cross=dict.fromkeys(regions)
dens_cross=dict.fromkeys(regions)

for region in regions:
    index[region]=np.intersect1d(index_1213[region],index_2110[region])
    temp_cross[region]=temp[index[region]]
    dens_cross[region]=dens[index[region]]


############################################################
# plot the figure.    
size=20
legend=regions
a,b=legend.index('spiral'),legend.index('restdisk')
legend[a],legend[b]=legend[b],legend[a]

fig=plt.figure(figsize=(10,7.5))
ax=fig.add_subplot(111)
ax.set_xlim((1.0,4.0))
ax.set_ylim((0,120))
ax.set_xlabel('log(number density) ($cm^{-3}$)',size=size)
ax.set_ylabel('Temperature (K)',size=size)
for region in regions:
    ax.plot(dens_cross[region],temp_cross[region],'.')
    ax.legend(legend,prop={'size':28},loc='upper left',markerscale=4.)

# ax.plot(dens_5257,temp_5257,'.')
ax.tick_params(labelsize=size)
ax.yaxis.set_major_locator(MaxNLocator(prune='both'))
ax.xaxis.set_major_locator(MaxNLocator(prune='both'))
plt.tight_layout()
fig.savefig(picturename)

#    plt.legend(loc="lower left", markerscale=2., scatterpoints=1, fontsize=10)
