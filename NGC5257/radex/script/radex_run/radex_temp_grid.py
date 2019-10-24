#! /usr/bin/python
#
import math
import os
import time
import numpy as np

############################################################
# set the parameter
abundance= 50 
tkin=20

cdmin = 1e15  
cdmax = 1e18

cdmin_13=cdmin_12/abundance
cdmax_13=cdmax_12/abundance

nmin = 1e1   
nmax = 1e4   

tbg   = 2.73
dv    = 1.0 

ncdmol  = 100
ndens = 100  
bw    = 0.001

gfil='co13_20K_50.dat'


############################################################
# define function

def write_input(infile,cdmol,nh2,flow,fupp,mole):
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

def write_infile(infile,cdmin,cdmax,flow,fupp,mole):
    for icdmol in range(ncdmol+1):
        for idens in range(ndens+1):
            cdmol = cdmin*((cdmax/cdmin)**(float(icdmol)/ncdmol))
            dens = nmin*((nmax/nmin)**(float(idens)/ndens))
            write_input(infile,cdmol,dens,flow,fupp,mole)
            if (icdmol == ncdmol and idens == ndens):
                infile.write('0\n')
                infile.close()
            else:
                infile.write('1\n')

def read_radex(outfile,fupp,flow):
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
    low = float(words[-2])
    while ((ftmp < fupp*(1-bw)) or (ftmp > fupp/(1-bw))):
        line  = outfile.readline()
        words = line.split()
        ftmp  = float(words[4])
    upp = float(words[-2])
    return cdmol,dens,upp,low


def write_dat(file,flow,fupp):
    cd=[];n=[];Ratio=[];rmin=100;rmax=0.1
    for icdmol in range(ncdmol):
        for idens in range(ndens):
            cdmol = cdmin*((cdmax/cdmin)**(float(icdmol)/ncdmol))
            dens = nmin*((nmax/nmin)**(float(idens)/ndens))
            output = read_radex(file,fupp,flow)
            cdmol=output[0];dens=output[1];ratio=output[2]
            cd.append(cdmol);n.append(dens);Ratio.append(ratio)
            if (ratio > 0.0):
                if (ratio < rmin):
                    rmin = ratio
                if (ratio > rmax):
                    rmax = ratio
    cd=np.array(cd)
    n=np.array(n);n=np.log10(n)
    Ratio=np.array(Ratio)
    result =np.transpose(np.vstack([cd,n,Ratio]))
    return result            

def write_dat_1213(file_12,file_13,freq_12,freq_13):
    Flux_12=[];Flux_13=[];cd=[];n=[]
    with open(file_12,'r') as outfile:
        flow=freq_12;fupp=freq_12
        for icdmol in range(ncdmol+1):
            for idens in range(ndens+1):
                cdmol = cdmin*((cdmax/cdmin)**(float(icdmol)/ncdmol))
                dens = nmin*((nmax/nmin)**(float(idens)/ndens))
                output=read_radex(outfile,fupp,flow)
                cdmol = output[0];dens=output[1];flux=output[2]
                cd.append(cdmol)
                n.append(dens)
                Flux_12.append(flux)
    with open(file_13,'r') as outfile:
        for icdmol in range(ncdmol+1):
            for idens in range(ndens+1):
                cdmol_13 = cdmin_13*((cdmax_13/cdmin_13)**(float(icdmol)/ncdmol))
                dens = nmin*((nmax/nmin)**(float(idens)/ndens))
                output = read_radex(outfile,freq_13,freq_13)
                flux=output[3]
                Flux_13.append(flux)
        cd=np.array(cd)
        n=np.array(n);n=np.log10(n)
        ratio=np.array(Flux_12)/np.array(Flux_13)
        result =np.transpose(np.vstack([cd,n,ratio]))
    return result

    
############################################################
# run the main program

# 12CO10 

flow=115.27;fupp=230.54;mole='co';file='co12_20K.dat'

with open('radex.inp','w') as infile:
    write_infile(infile,cdmin_12,cdmax_12,flow,fupp,mole)

os.system('radex < radex.inp > /dev/null')

outfile=open('radex.out','r')
result_2110=write_dat(outfile,flow,fupp)
outfile.close()

# 13CO10 

freq_12=115.27;freq_13=110.2;mole='13co';file='co12-13_20K.dat'

with open('radex.inp','w') as infile:
    write_infile(infile,cdmin_13,cdmax_13,freq_13,freq_13,mole)
os.system('radex < radex.inp > /dev/null')

file_12='radex.out';file_13='radex_13co.out'
result_1213=write_dat_1213(file_12,file_13,freq_12,freq_13)

