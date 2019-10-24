'''Jan 31th, 2018
test casa version is 5.1.1

Nov 4th, 2018

switch to casa 5.4.0,
rms=1.6e-3


'''

import numpy as np
import os
import scipy.ndimage as sni
import sys

vis = '/home/heh15/workingspace/Arp240/' \
      '12CO10/calibrated/NGC5257_combine_CO.ms.contsub'
cleanDir = '/home/heh15/workingspace/Arp240/'\
             '12CO10/NGC5257/casa5.4/'
preName = cleanDir + 'NGC5257_12CO10_combine_contsub'
field = '0'
phasecenter='J2000 13h39m52.922 0d50m24.1'
mode = 'velocity'
restfreq='112.73GHz'
# restfreq='115.27GHz'
width = '10km/s' 
nchan = 70 
start = '-300km/s' 
cell='0.3arcsec'  
imsize = [320,320]
weighting = 'briggs'
robust = 0.5
imagermode = 'mosaic'
cyclefactor = 1.0  # default value
# stuff related to auto-masking
rmsCRTF = '/1/home/heh15/workingspace/Arp240/' \
          '12CO10/region/NGC5257_emission_free_54.crtf'
rmsChan = '60~69'

# Additional parameter for tclean
specmode='cube'
outframe='BARY'
deconvolver='hogbom' 
gridder='mosaic' 
pblimit=0.2

# Additional parameter for multithresh
sidelobethreshold=3.0
noisethreshold=4.2
lownoisethreshold=1.5
negativethreshold=0.0

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
       restoringbeam=restoringbeam,
       cyclefactor=cyclefactor,
       pblimit=pblimit,
       interactive=False)

# Make moment map

# immoments(imagename=myImage,moments=0,chans='40~90',includepix=[2*rms,100],outfile=myImage+'.mom0')

# go back to where we started
os.chdir(origDir)



