#! /usr/bin/python
#
import math
import os
import time
import temp_header


# modified by Hao He at Mcmaster University
'''
############################################################
# set the parameters.

tkin=20

cdmin = 1e15  
cdmax = 1e18
nmin = 1e1   
nmax = 1e4   

tbg   = 2.73
dv    = 1.0 

ncdmol  = 100
ndens = 100  
bw    = 0.001

type='2110'


############################################################
# define function

def write_input(infile,cdmol,nh2):
    infile.write(mole+'.dat\n')
    infile.write('radex.out\n')
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
    while (words[3] != "H2"):
        line  = outfile.readline()
        words = line.split()
    dens  = float(words[-1])
    while (words[1] != "Column"):
        line = outfile.readline()
        words= line.split()
    cdmol= float(words[-1])
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
    low   = float(words[-2])
    line  = outfile.readline()
    words = line.split()
    if fupp == flow:
        upp = low
    else:
        ftmp  = float(words[4])
        while ((ftmp < fupp*(1-bw)) or (ftmp > fupp/(1-bw))):
            line  = outfile.readline()
            words = line.split()
            ftmp  = float(words[4])
        upp   = float(words[-2])
    ratio = low / upp
    return cdmol,dens,upp,low

def namefile(tkin,type):
    datafile='co'+type+str(tkin)+'K'+'.dat'
'''
############################################################ 
# Begin of main program

start = time.time()

gfil=namefile(tkin,type)

mole = 'co'
act=[115.27,230.54]

flow = act[0]
fupp = act[1]

write_infile(infile,cdmin_12,cdmax_12,flow,fupp,mole)

infile = open('radex.inp','w')
print "Starting ",gfil

for icdmol in range(ncdmol+1):
    for idens in range(ndens+1):
        cdmol = cdmolmin*((cdmolmax/cdmolmin)**(float(icdmol)/ncdmol))
        dens = nmin*((nmax/nmin)**(float(idens)/ndens))
        write_input(infile,cdmol,dens)
        if (icdmol == ncdmol and idens == ndens):
            infile.write('0\n')
            infile.close()
        else:
            infile.write('1\n')

os.system('radex < radex.inp > /dev/null')
grid = open(gfil,'w')
fmt  = '%10.3e %10.3e %10.3e \n'

outfile  = open('radex.out')

rmin = 100
rmax = 0.1

for icdmol in range(ncdmol+1):
    for idens in range(ndens+1):
        cdmol = cdmolmin*((cdmolmax/cdmolmin)**(float(icdmol)/ncdmol))
        dens = nmin*((nmax/nmin)**(float(idens)/ndens))
        cdmol,dens,ratio = read_radex(outfile)
        if (ratio > 0.0):
            if (ratio < rmin):
                rmin = ratio
            if (ratio > rmax):
                rmax = ratio
        grid.write(fmt %(math.log10(cdmol), math.log10(dens), ratio))

grid.close()
outfile.close()

print "Min, max:", rmin, rmax

stop = time.time()
dure = stop - start
print "Run time = ",dure, "seconds"


start = time.time()

freq_12=115.27
freq_13=110.20
file_12='radex.out'
file_13='radex_13co.out'

mole = '13co'
with open('radex.inp','w') as infile:
    write_infile(infile,cdmin_13,cdmax_13,freq_13,freq_13,mole)

print "Starting ",gfil


os.system('radex < radex.inp > /dev/null')

stop = time.time()
dure = stop - start
print "Run time = ",dure, "seconds"

result = write_dat_1213(file_12,file_13,freq_12,freq_13)
fmt  = '%10.3e %10.3e %10.3e \n'
np.savetxt(gfil,output,fmt=fmt,delimiter='',newline='\n')
