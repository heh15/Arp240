'''
May 5th, 2018
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
    data=np.loadtxt(filename)
    cdmol_tmp=data[:,0]
    dens_tmp=data[:,1]
    ratio=data[:,2]
    if mol=='2110':
        ratio=np.power(ratio,-1)
    cdmol=dict.fromkeys(regions)
    dens=dict.fromkeys(regions)
    index=dict.fromkeys(regions)
    for region in regions:
        index[region]=np.where(np.abs(ratio-ratio_table[mol+'_ratio'][region]) < ratio_table[mol+'_uncertainty'][region])
        cdmol[region]=cdmol_tmp[index[region]]
        dens[region]=dens_tmp[index[region]]
    return index, cdmol, dens


def get_file(col,abu):
    co_12='co_1-0_2-1'
    co_13='12co_13co'
    co_12=co_12+'_'+col
    co_13=co_13+'_'+col
    co_13=co_13.replace(col,abu+'_'+col)
    co_12=co_12+'.dat'
    co_13=co_13+'.dat'
    picturename='radex_'+col+'_'+abu+'.png'
    return co_12,co_13,picturename


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
picturename=picDir+'NGC5258_uncal_'+get_file(col,abu)[2]

os.chdir(logDir)


# read the data from the ratio data file.
file_1213='../log/NGC5258_1213_ratio_quadrature.txt'
file_2110='../log/NGC5258_2110_ratio_quadrature.txt'

regions=['northarm','southarm','center','ring']
values=['1213_ratio','1213_uncertainty','2110_ratio','2110_uncertainty']
variables=['cdmol','dens']

ratio_table=pd.DataFrame()
ratio_table[0]=read_ratio(file_1213)[0]
ratio_table[1]=read_ratio(file_1213)[1]
ratio_table[2]=read_ratio(file_2110)[0]
ratio_table[3]=read_ratio(file_2110)[1]
ratio_table.index=regions
ratio_table.columns=values


os.chdir(scriptDir)


############################################################
# load data

filename=co_13
mol='1213'
index_1213,cdmol_1213,dens_1213 = read_value(filename,regions,ratio_table,mol)


filename='co_20K.dat'
mol='2110'
index_2110,cdmol_2110,dens_2110 = read_value(filename,regions,ratio_table,mol)

region='northarm'
fig=plt.figure(figsize=(10,7.5))
plt.plot(cdmol_2110[region],dens_2110[region],linestyle='None',marker='o')
plt.show()

'''
############################################################
# find the cross region in parameter space. 

index=dict.fromkeys(regions)
cdmol_cross=dict.fromkeys(regions)
dens_cross=dict.fromkeys(regions)

for region in regions:
    index[region]=np.intersect1d(index_1213[region],index_2110[region])
    cdmol_cross[region]=cdmol[index[region]]
    dens_cross[region]=dens[index[region]]


############################################################
# plot the figure.    
size=32
legend=['northarm','southarm','center','ring']

fig=plt.figure(figsize=(10,7.5))
ax=fig.add_subplot(111)
ax.set_xlim((1.0,4.0))
ax.set_ylim((0,120))
ax.set_xlabel('log(number density) ($cm^{-3}$)',size=size)
ax.set_ylabel('Cdmolerature (K)',size=size)
for region in regions:
    ax.plot(dens_cross[region],cdmol_cross[region],'.')
    ax.legend(legend,prop={'size':28},loc='upper left')
# ax.plot(dens_5257,cdmol_5257,'.')
ax.tick_params(labelsize=size)
ax.yaxis.set_major_locator(MaxNLocator(prune='both'))
ax.xaxis.set_major_locator(MaxNLocator(prune='both'))
plt.tight_layout()
fig.savefig(picturename)
'''
