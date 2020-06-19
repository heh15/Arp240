import time
import os
from subprocess import Popen, PIPE, STDOUT
from subprocess import call
import glob
import numpy as np
import copy

def getval(input):
    line=input.readline()
    words=line.split()
    marker=words[0]
    while(marker != "Best"):
        line=input.readline()
        words=line.split()
        if words==[]:
            continue
        marker=words[0]
    line=input.readline()
    line=input.readline()
    words=line.split()
    PA=words[4]
    line=input.readline()
    words=line.split()
    eps=words[2]
    line=input.readline()
    words=line.split()
    incl=words[3]
    line=input.readline()
    words=line.split()
    xcenter=words[4]
    ycenter=words[-3]
    line=input.readline()
    line=input.readline()
    words=line.split()
    vsys=words[2]
    return PA, eps, incl, xcenter, ycenter, vsys

def modify(lines,PA, eps, xcenter, ycenter, vsys,i):
    newlines=copy.copy(lines)
    words=newlines[6].split()
    words[1]=PA
    words[2]=eps
    newlines[6]='  '.join(words)+'\n'
    newlines[7]=outfile[i]+'\n'
    words=newlines[9].split()
    words[0]=PA
    words[1]=eps
    newlines[9]='  '.join(words)+'\n'
    words=newlines[10].split()
    words[0]=xcenter
    words[1]=ycenter
    newlines[10]='  '.join(words)+'\n'
    words=newlines[13].split()
    words[1]=vsys
    newlines[13]='  '.join(words)+'\n'
    return newlines

        
i=0;n=10
infile=[];outfile=[]
for i in range(n):
    infile.append('velsf_12CO21_'+str(i)+'.inp')
    outfile.append('NGC5257_12CO21_vel_'+str(i)+'.out')

origfile='velsf_12CO21_orig.inp'
with open(origfile,'r') as input:
    content=input.readlines()
content=''.join(content)
with open(infile[0],'w') as output:
    output.write(content)

incls=[]
epses=[]

for i in range(n):
    if i !=0:
        newlines=modify(lines,PA, eps, xcenter, ycenter, vsys,i)
        content=''.join(newlines)
        with open(infile[i],'w') as output:
            output.write(content)
        call(['find', '.','-name','NGC5257_12CO21*','-and','!','-name','NGC5257_12CO21*.out','-delete'])
# run the file
    filename=infile[i].encode('utf-8')
    p = Popen(['./DiskFit'], stdout=PIPE, stdin=PIPE, stderr=STDOUT)  
    p.communicate(input=filename)
# read the input file and output file
    with open(infile[i],'r') as input:
        lines=input.readlines()
    with open(outfile[i],'r') as input:
        PA, eps, incl, xcenter, ycenter, vsys=getval(input)
    incls.append(incl)
    epses.append(eps)

    
    
