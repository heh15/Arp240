'''Jan 31 th
test casa version is 5.1.1

Nov 4th, 2018 
switch to casa version 5.4.0
the phase center is changed to coordinates. 
rms=0.000644
'''

import numpy as np
import os
import scipy.ndimage as sni
import sys

vis = '/home/heh15/workingspace/Arp240/' \
      '13CO10/calibrated/NGC5258_CO13_12m_contsub.ms'
cleanDir = '/home/heh15/workingspace/Arp240/'\
             '13CO10/NGC5258/casa5.4/'
preName = cleanDir + 'NGC5258_CO13_10_12m_mthresh_split_contsub_noise45'
field = '0'
phasecenter='J2000 13h39m57.675 0d49m51.5'
mode = 'velocity'
restfreq='107.78GHz'
width = '20km/s' 
nchan = 35 
start = '-300km/s' 
cell='0.3arcsec'  
imsize = [320,320]
weighting = 'briggs'
robust = 0.5
imagermode = 'mosaic'
cyclefactor = 1.0  # default value
# stuff related to auto-masking
rmsCRTF = '/1/home/heh15/workingspace/Arp240/13CO10/' \
          'NGC5258/emission_free_54.crtf'
rmsChan = '1~5'

# Additional parameter for tclean
specmode='cube'
outframe='BARY'
deconvolver='hogbom' 
gridder='mosaic' 
pblimit=0.2


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
      field=field,
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

# Additional parameter for multithresh
sidelobethreshold=3.0
noisethreshold=4.2
lownoisethreshold=1.5
negativethreshold=0.0
minbeamfrac=0.2

# set the environment
# os.environ['SAVE_ALL_AUTOMASKS']="true"

tclean(vis=vis,
       imagename=preName,
       field=field,
       antenna=antenna,
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
       niter=100000,
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

# immoments(imagename=myImage,moments=0,includepix=[2*rms,100],chans='5~30',outfile=myImage+'.mom0')

# go back to where we started
os.chdir(origDir)







