'''Jan 31 th
test casa version is 5.1.1

Nov. 4th, 2018
switch to casa version 5.4.0
rms=0.00301

'''

import numpy as np
import os
import scipy.ndimage as sni
import sys

vis = '/home/heh15/workingspace/Arp240/12CO21/'\
      'calibrated/NGC5258_CO12_21_combine.ms/'
cleanDir = '/home/heh15/workingspace/Arp240/'\
             '12CO21/NGC5258/casa5.4/test/'
preName = cleanDir + 'NGC5258_12CO21_combine_noise45'
rmsCRTF = '/home/heh15/workingspace/Arp240/12CO21/'\
          '/NGC5258/emission_free_54.crtf'
rmsChan = '60~69'
pblimit=0.2
phasecenter='J2000 13h39m57.675 0d49m51.5'
spw='0~3'

# Additional parameter for multithresh
sidelobethreshold=4.0
noisethreshold=4.5
lownoisethreshold=1.5
negativethreshold=0.0
minbeamfrac=0.3
niter=100000

# default
mode = 'velocity'
restfreq='225.46GHz'
width = '10km/s' 
nchan = 70 
start = '-300km/s' 
cell='0.1arcsec'  
imsize = 960
weighting = 'briggs'
robust = 0.5
imagermode = 'mosaic'
cyclefactor = 1.0  # default value
specmode='cube'
outframe='BARY'
deconvolver='hogbom' 
gridder='mosaic'

# CLEAN output names
myImage = preName + '.image'
myFlux = preName + '.flux'
myMask = preName + '.mask'
myResidual = preName + '.residual'
 
# change directory to workaround makemask failing
# when specifying paths
origDir = os.getcwd()
os.chdir(cleanDir)

############################################################
# Make dirty map

delmod(vis=vis)

tclean(vis=vis,
      imagename=preName,
      phasecenter=phasecenter,
      spw=spw,
      specmode='cube',
      outframe='BARY',
      restfreq=restfreq,
      width=width,
      nchan=nchan,
      start=start,
      cell=cell,
      imsize=imsize,
      weighting=weighting,
      robust=robust,
      deconvolver='hogbom',
      gridder='mosaic',
      niter=None,
      threshold='0.0Jy',
      cyclefactor=cyclefactor,
      pblimit=0.2,
      interactive=False)

# measure noise region

if rmsCRTF == '':
    done = False
    rmsCRTF = raw_input('Go draw and save a CASA'
                        +'viewer region for '
                        +'estimating the RMS in the '
                        +'dirty image \n'+myImage
                        +'.\nEnter the full path to'
                        +'that CRTF and press '
                        +'return.\n')
    while not done:
        if not os.path.exists(rmsCRTF):
            rmsCRTF = raw_input(rmsCRTF+' does not'
                                +'exist, try again.'
                                +'\n')
        else:
            print rmsCRTF, 'successfully found.', \
                  'Continuing...'
        done = True
rms = imstat(imagename=preName+'.image',
             region=rmsCRTF, chans=rmsChan)['rms'][0]

# grab some details from the dirty image
peak = imstat(imagename=myImage)['max'][0]
thresh = 2*rms

delmod(vis=vis)

# set the environment
# os.environ['SAVE_ALL_AUTOMASKS']="true"

tclean(vis=vis,
       imagename=preName,
       spw=spw,
       phasecenter=phasecenter,
       specmode=specmode,
       outframe=outframe,
       restfreq=restfreq,
       width=width,
       nchan=nchan,
       start=start,
       deconvolver=deconvolver,
       cell=cell,
       imsize=imsize,
       weighting=weighting,
       robust=robust,
       gridder=gridder,
       usemask='auto-multithresh',
       niter=niter,
       threshold=str(thresh)+'Jy/beam',
       sidelobethreshold=sidelobethreshold,
       noisethreshold=noisethreshold,
       lownoisethreshold=lownoisethreshold,
       negativethreshold=negativethreshold,
       minbeamfrac=minbeamfrac,
       restoringbeam=restoringbeam,
       cyclefactor=cyclefactor,
       pblimit=pblimit,
       interactive=False)

# Make moment map

immoments(imagename=myImage,moments=0,includepix=[2*rms,100],chans='10~60',outfile=myImage+'.mom0')

# go back to where we started
os.chdir(origDir)







