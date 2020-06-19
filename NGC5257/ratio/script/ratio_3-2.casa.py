'''
Nov.23rd, 2018

smooth the image to the same beam size
'''
import glob, os
import re
from distutils.dir_util import copy_tree

Dir='/home/heh15/workingspace/Arp240/NGC5257/ratio/'
workDir=Dir+'CO3-2/'
scriptDir=Dir+'script'
picDir=Dir+'picture/'
logDir=Dir+'log/'

############################################################
# basic settings

############################################################
# function

############################################################
# program

os.chdir(workDir)

files=glob.glob('NGC5257_*')

for imagename in files:
    imsmooth(imagename=imagename,
             kernel='gauss',
             major='3.789arcsec',
             minor='2.989arcsec',
             pa='-17.357deg',
             targetres=True,
             outfile=imagename.replace('.image','_smooth.image'))

images_10km=glob.glob('NGC5257_12CO*_smooth.image')

for imagename in images_10km:
    immoments(imagename=imagename,
              moments=0,
              chans='10~60',
              outfile=imagename.replace('.image','.mom0'))

image_20km=glob.glob('NGC5257_CO13*_smooth.image')
for imagename in image_20km:
   immoments(imagename=imagename,
          moments=0,
          chans='5~30',
          outfile=imagename.replace('.image','.mom0'))    
