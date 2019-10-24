import time
import os
from subprocess import Popen, PIPE, STDOUT
from subprocess import call
import glob
import numpy as np
import copy

############################################################
# directory
Dir='/home/heh15/workingspace/Arp240/NGC5258/RotationCurve/'
diskfitDir=Dir+'diskfit/'
C12O_10=diskfitDir+'12CO10/'
C12O_21=diskfitDir+'12CO21/'
logDir=Dir+'log/'

############################################################
# main program

os.chdir(diskfitDir)

#### 12CO10 ####

## run the file
infile=diskfitDir+'velsf_5258_12CO10.inp'
filename=infile.encode('utf-8')
p = Popen(['./DiskFit.bin'], stdout=PIPE, stdin=PIPE, stderr=STDOUT)  
p.communicate(input=filename)

## copy the file to the 
call(['cp',C12O_10+'NGC5258_12CO10_vel_radcut.out',logDir])

#### 12CO2-1 ####
infile=diskfitDir+'velsf_5258_12CO21.inp'
filename=infile.encode('utf-8')
p = Popen(['./DiskFit.bin'], stdout=PIPE, stdin=PIPE, stderr=STDOUT)  
p.communicate(input=filename)

## copy the file to the 
call(['cp',C12O_21+'NGC5258_12CO21_vel_radcut.out',logDir])

    

    
