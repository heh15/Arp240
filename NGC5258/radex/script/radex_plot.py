'''
May 5th, 2018
'''


import numpy as np
from shutil import copy
import os
from matplotlib.ticker import MaxNLocator
from matplotlib import rc
from matplotlib import rcParams
rcParams['mathtext.default']='regular'

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
abu='30'

# basic setting
Dir='/home/heh15/workingspace/Arp240/radex/'
picDir=Dir+'picture/'
logDir=Dir+'log/'
scriptDir=Dir+'script/radex_run/'
co_12=get_file(col,abu)[0]
co_13=get_file(col,abu)[1]
picturename=picDir+'NGC5258_cal_'+get_file(col,abu)[2]

# copy the ratio data file into the current directory
#filename='/home/heh15/workingspace/Arp240/ratio/script/NGC5258_2110_ratio_quadrature.txt'
#copy(filename,logDir)
#filename='/home/heh15/workingspace/Arp240/ratio/script/NGC5258_1213_ratio_quadrature.txt'
#copy(filename,logDir)


# read the data from the ratio data file.
dir='../../ratio/log/'
ratio_file={}
ratio_file['1213']=dir+'NGC5258_1213_ratio_error_cal.txt'
ratio_file['2110']=dir+'NGC5258_2110_ratio_error_cal.txt'

regions=['north','south','center','ring']
values=['flux','uncertainty']

ratio_1213=read_ratio(ratio_file['1213'])[0]               
uncertainty_1213=read_ratio(ratio_file['1213'])[1]

ratio_2110=read_ratio(ratio_file['2110'])[0]
uncertainty_2110=read_ratio(ratio_file['2110'])[1]

os.chdir(scriptDir)
                                                    
############################################################
# load the 1-0_2-1 data from radex file

# load radex output file
filename=co_12
data=np.loadtxt(filename)
temp=data[:,0]
dens=data[:,1]
ratio_radex_2110=data[:,2]
ratio_radex_2110=np.power(ratio_radex_2110,-1)


# draw the region. 
temp_2110=dict.fromkeys(regions)
dens_2110=dict.fromkeys(regions)


for i in range(len(regions)):
    index=np.where(np.abs(ratio_radex_2110-ratio_2110[i])< uncertainty_2110[i])
    temp_2110[regions[i]]=temp[index]
    dens_2110[regions[i]]=dens[index]

#fig=plt.figure()
#for i in range(len(regions)):
#    plt.plot(dens_2110[regions[i]],temp_2110[regions[i]],'r.')



############################################################
# load the 12co-13co data

filename=co_13

data=np.loadtxt(filename)
temp=data[:,0]
dens=data[:,1]
ratio_radex_1213=data[:,2]

# draw the region. 
temp_1213=dict.fromkeys(regions)
dens_1213=dict.fromkeys(regions)


for i in range(len(regions)):
    index=np.where(np.abs(ratio_radex_1213-ratio_1213[i])< uncertainty_1213[i])
    temp_1213[regions[i]]=temp[index]
    dens_1213[regions[i]]=dens[index]

#for i in range(len(regions)):
#    plt.plot(dens_1213[regions[i]],temp_1213[regions[i]],'b.')


############################################################
# find the cross region in parameter space. 

# draw the region. 
temp_op=dict.fromkeys(regions)
dens_op=dict.fromkeys(regions)

for i in range(len(regions)):
    mask=((np.abs(ratio_radex_2110-ratio_2110[i])<uncertainty_2110[i]) & \
          (np.abs(ratio_radex_1213-ratio_1213[i])< uncertainty_1213[i]))
    mask=np.logical_not(mask)
    temp_op[regions[i]]=np.ma.masked_where(mask,temp)
    dens_op[regions[i]]=np.ma.masked_where(mask,dens)

############################################################
# do the same thing for NGC 5257

ratio_1213=11.19
error_1213=1.11
ratio_2110=0.71
error_2110=0.05

mask=((np.abs(ratio_radex_2110-ratio_2110)<error_2110) & \
      (np.abs(ratio_radex_1213-ratio_1213)< error_1213))
mask=np.logical_not(mask)
temp_5257=np.ma.masked_where(mask,temp)
dens_5257=np.ma.masked_where(mask,dens)


############################################################
# plot the figure.    
size=32
legend=['northarm','southarm','center','ring']

fig=plt.figure(figsize=(10,7.5))
ax=fig.add_subplot(111)
ax.set_xlim((1.0,4.0))
ax.set_ylim((0,120))
ax.set_xlabel('log(number density) ($cm^{-3}$)',size=size)
ax.set_ylabel('Temperature (K)',size=size)
for i in range(len(regions)):
    ax.plot(dens_op[regions[i]],temp_op[regions[i]],'o')
    ax.legend(legend,prop={'size':28},loc='upper left',markerscale=4.0)
# ax.plot(dens_5257,temp_5257,'.')
ax.tick_params(labelsize=size)
ax.yaxis.set_major_locator(MaxNLocator(prune='both'))
ax.xaxis.set_major_locator(MaxNLocator(prune='both'))
plt.tight_layout()
fig.savefig(picturename)

os.chdir('..')
