#! /usr/bin/python
#
import math
import os
import time
import numpy as np
from globalvar import col 
from globalvar import abundance

# parameter
cdmol=col/abundance
type='12co_13co'


# Grid boundaries
#
tmin = 5.0  # minimum kinetic temperature (K)
tmax = 105.0 # maximum kinetic temperature (K)
nmin = 1e1   # minimum H2 density (cm^-3)
nmax = 1e4   # maximum H2 density (cm^-3)
#
# Parameters to keep constant
#
tbg   = 2.73 # background radiation temperature
dv    = 1.0 # line width (km/s)
#
# Numerical parameters
#
ntemp = 100  # number of temperature points
ndens = 100  # number of density points
bw    = 0.001 # "bandwidth": free spectral range around line
#
# No user changes needed below this point.
#
def write_input(infile,tkin,nh2):
    infile.write(mole+'.dat\n')
    infile.write('radex_13co.out\n')
    infile.write(str(flow*(1-bw))+' '+str(fupp/(1-bw))+'\n')
    infile.write(str(tkin)+'\n')
    infile.write('1\n')
    infile.write('H2\n')
    infile.write(str(nh2)+'\n')
    infile.write(str(tbg)+'\n')
    infile.write(str(cdmol)+'\n')
    infile.write(str(dv)+'\n')

def read_radex(outfile):
    line  = outfile.readline()
    words = line.split()
    while (words[1] != "T(kin)"):
        line  = outfile.readline()
        words = line.split()
    temp  = float(words[-1])
    line  = outfile.readline()
    words = line.split()
    dens  = float(words[-1])
    while (words[-1] != "FLUX"):
        line  = outfile.readline()
        words = line.split()
    line  = outfile.readline()
    line  = outfile.readline()
    words = line.split()
    ftmp  = float(words[4])
    while ((ftmp < flow*(1-bw)) or (ftmp > flow/(1-bw))):
        line  = outfile.readline()
        words = line.split()
        ftmp  = float(words[4])
    flux   = float(words[-2])
    return temp,dens,flux
 
def output_filename(col):
    colsub=str(int(round(math.log(col,10)*10,0)))
    filename=type+'_'+colsub
    filename=filename.replace(colsub,str(abundance)+'_'+colsub)
    filename=filename+'.dat'
    return filename

# Begin of main program


start = time.time()

mole = '13co'
act=[110.20,110.20]

flow = act[0]
fupp = act[1]
gfil=output_filename(col)

infile = open('radex.inp','w')
print "Starting ",gfil

for itemp in range(ntemp+1):
    for idens in range(ndens+1):

        temp = tmin*((tmax/tmin)**(float(itemp)/ntemp))
        dens = nmin*((nmax/nmin)**(float(idens)/ndens))
        write_input(infile,temp,dens)
        if (itemp == ntemp and idens == ndens):
            infile.write('0\n')
            infile.close()
        else:
            infile.write('1\n')

os.system('radex < radex.inp > /dev/null')


stop = time.time()
dure = stop - start
print "Run time = ",dure, "seconds"


# write down the 13co temperature, density and flux
T=[]
n=[]
Flux_13=[]

outfile  = open('radex_13co.out')

for itemp in range(ntemp+1):
    for idens in range(ndens+1):

        temp = tmin*((tmax/tmin)**(float(itemp)/ntemp))
        dens = nmin*((nmax/nmin)**(float(idens)/ndens))
        temp,dens,flux = read_radex(outfile)
        T.append(temp)
        n.append(dens)
        Flux_13.append(flux)

outfile.close()


# write down the 12co temperature, density and flux
act=[115.27,115.27]

flow = act[0]
fupp = act[1]

T=[]
n=[]
Flux_12=[]

outfile  = open('radex.out')

for itemp in range(ntemp+1):
    for idens in range(ndens+1):

        temp = tmin*((tmax/tmin)**(float(itemp)/ntemp))
        dens = nmin*((nmax/nmin)**(float(idens)/ndens))
        temp,dens,flux = read_radex(outfile)
        T.append(temp)
        n.append(dens)
        Flux_12.append(flux)

outfile.close()


# write down the temperature, density and ratio.
T=np.array(T)
n=np.array(n);n=np.log10(n)
ratio=np.array(Flux_12)/np.array(Flux_13)

fmt  = '%10.3e %10.3e %10.3e \n'

output=np.transpose(np.vstack([T,n,ratio]))

np.savetxt(gfil,output,fmt=fmt,delimiter='',newline='\n')
