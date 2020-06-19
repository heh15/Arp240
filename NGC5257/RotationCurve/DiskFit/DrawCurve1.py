import time
import matplotlib.pyplot as plt
import numpy as np
import math

def readvels(infile):
    radius=[];velocity=[];velocity2=[]
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
        velocity2.append(words[6])
        line=infile.readline()
        words=line.split()
    radius=[float(l) for l in radius]
    velocity=[float(l) for l in velocity]
    velocity2=[float(l) for l in velocity2]
    return radius, velocity,velocity2

def getvel(velocity,velocity2,theta,phi_b):
    velocity=np.array(velocity)
    velocity2=np.array(velocity2)
    vel=np.arange(len(velocity),dtype='float')
    for i in range(len(velocity)):
        vel[i]=velocity[i]+velocity2[i]*math.cos((2*theta+180-2*phi_b)*math.pi/180)
    return vel

start=time.time()
i=0
radius=[];velocity=[]
filename='NGC5257_12CO10_vel.out'
with open (filename,'r') as Input:
    radius,velocity,velocity2=readvels(Input)

phi_b=65.72;theta=0;theta2=90
radius_arcsec=[0.3*l for l in radius]
vel1=getvel(velocity, velocity2,theta,phi_b)
vel2=getvel(velocity, velocity2,theta2,phi_b)


#plt.scatter(radius_arcsec,vel1)
plt.scatter(radius_arcsec,vel2)

stop=time.time()
t=stop-start

filename='NGC5257_12CO21_vel.out'
with open (filename,'r') as Input:
    radius21,velocity21,velocity21_2=readvels(Input)
radius_arcsec21=[0.1*l for l in radius21]
vel21=getvel(velocity21,velocity21_2,theta,phi_b)

# plt.scatter(radius_arcsec21,vel21)

plt.show()
