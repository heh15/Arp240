import time
import matplotlib.pyplot as plt
import numpy as np
import math
from matplotlib import rcParams
rcParams['mathtext.default']='regular'

############################################################
# Basic settings

############################################################ 
# function 

# read the diskfile output
def readvels(infile):
    radius=[];velocity=[];error=[]
    line=infile.readline()
    words=line.split()
    marker=words[0]
    while (marker != "Fitted"):
        line=infile.readline()
        words=line.split()
        if words==[]:
            continue
        marker=words[0]
    line=infile.readline()
    line=infile.readline()
    line=infile.readline()
    words=line.split()
    while(words != []):
        radius.append(words[0])
        velocity.append(words[2])
        error.append(words[3])
        line=infile.readline()
        words=line.split()
    radius=0.3*np.array([float(l) for l in radius])
    velocity=np.array([float(l) for l in velocity])
    error=np.array([float(l) for l in error])
    return radius, velocity,error

        
############################################################
# main program

Dir='/1/home/heh15/workingspace/Arp240/NGC5257/RotationCurve/'
logDir=Dir+'log/'
picDir=Dir+'picture/'

i=0
radius=[];velocity=[]
filename=Dir+'DiskFit/NGC5257/output_radcut/NGC5257_12CO10_vel.out'
with open (filename,'r') as Input:
    radius,velocity,error=readvels(Input)
radius_arcsec=radius
vel=velocity
error_CO10=error

i=0
radius=[];velocity=[]
filename=Dir+'DiskFit/NGC5257/output_radcut/NGC5257_12CO21_vel.out'
with open (filename,'r') as Input:
    radius,velocity,error=readvels(Input)

radius_2d_CO21=radius/3
vel_CO21=velocity
error_CO21=error
radius21=radius_2d_CO21


filename=Dir+'Bbarolo/12CO21/test5/run3/ringlog1.txt'
with open (filename,'r') as Input:
    header=Input.readline()
    header=header.split()
    data=np.loadtxt(Input,skiprows=0)
data=np.transpose(data)
radius_3dB=data[1]
vel_3dB=data[2]
error_low_3dB=-data[13]
error_upp_3dB=data[14]

# filename='Bbarolo/NGC5257_12CO21/ringlog1.txt'
# with open (filename,'r') as Input:
#     header=Input.readline()
#     header=header.split()
#     data=np.loadtxt(Input,skiprows=0)
# data=np.transpose(data)
# radius_3dB_CO21=data[1]

vel_3dB_CO21=data[2]


# R_ba=np.arange(1,15,2)
# V_rot=np.array([162,162,219,261,294,315,348])

# R_ba21=np.arange(0.5,14,1)
# V_rot21=np.array([156,150,135,171,210,228,234,252,261,285,285,336,363,279])

# read the data from the literature
data=np.loadtxt(logDir+'rotationcurve_app.txt')
data=data.transpose()

data2=np.loadtxt(logDir+'rotationcurve_rec.txt')
data2=data2.transpose()

 # draw the figure. 
fig=plt.figure()
CO10_2d=plt.errorbar(radius_arcsec,vel,error_CO10,color='red',marker='o',linestyle='None',label='DiskFit 12CO10')
plt.show()
CO21_2d=plt.errorbar(radius_2d_CO21,vel_CO21,error_CO21,color='blue',marker='o',linestyle='None',label='DiskFit 12CO21')
#plt.legend([CO10_2d,CO21_2d],['CO10 Diskfit','CO21 Diskfit'],loc='upper left')
# CO10_3d=plt.scatter(R_ba,V_rot,marker='o',color='orange',label='Bbarolo 12CO10')
CO21_3d=plt.errorbar(radius_3dB,vel_3dB,[error_low_3dB,error_upp_3dB],marker='o',color='green',label='Bbarolo 12CO21')
# Ha_rec=plt.errorbar(data[0],data[1],data[2],color='black',marker='h',linestyle='none',label='Ha approaching side')
Ha_rec=plt.errorbar(data[0],data[1],color='black',marker='*',linestyle='none',label='Ha approaching side')
# Ha_app=plt.errorbar(data2[0],data2[1],data2[2],color='black',marker='*',linestyle='none',label='Ha approaching side')
Ha_app=plt.plot(data2[0],data2[1],color='black',marker='h',linestyle='none',label='Ha receding side')
plt.legend(loc='lower right')
plt.xlabel('radius(arcsec)')
plt.ylabel('velocity(km/s)')
plt.show()
plt.savefig(picDir+'NGC5257_12CO_rot.png')

# output the velocity result
output=np.transpose(np.vstack((radius_arcsec,vel,error_CO10)))
fmt= '%10.3e %10.3e %10.3e \n'
filename=logDir+'vel_12CO10.txt'
np.savetxt(filename,output,fmt=fmt,delimiter='',newline='\n')

output=np.transpose(np.vstack((radius_2d_CO21,vel_CO21,error_CO21)))
fmt= '%10.3e %10.3e %10.3e \n'
filename='log/vel_12CO21.txt'
np.savetxt(filename,output,fmt=fmt,delimiter='',newline='\n')

