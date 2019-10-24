'''Apr.10th, 2019
Try to image the CN1-0 line. 
rms=6.8e-4


'''

import numpy as np
import os
import scipy.ndimage as sni
import sys

vis = '/home/heh15/workingspace/Arp240/NGC5257/' \
      '12CO10/calibrated/NGC5257_combine.ms'
cleanDir = '/home/heh15/workingspace/Arp240/NGC5257/'\
             'CN10/cleanDir/'
preName = cleanDir + 'NGC5257_CN10_combine'
field = '0'
phasecenter='J2000 13h39m52.922 0d50m24.1'
mode = 'velocity'
restfreq='110.97GHz'
# restfreq='115.27GHz'
width = '40km/s' 
nchan = 30
start = '-600km/s' 
cell='0.3arcsec'  
imsize = [320,320]
weighting = 'briggs'
robust = 0.5
imagermode = 'mosaic'
cyclefactor = 1.0  # default value
spw='1,3'

# stuff related to auto-masking
rmsCRTF = '/1/home/heh15/workingspace/Arp240/NGC5257/' \
          '12CO10/region/NGC5257_emission_free_54.crtf'
rmsChan = '0~10'

# Additional parameter for tclean
specmode='cube'
outframe='BARY'
deconvolver='hogbom' 
gridder='mosaic' 
pblimit=0.2

# Additional parameter for multithresh
sidelobethreshold=3.0
noisethreshold=5.0
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
      interactive=False,
       spw=spw)

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
       interactive=False,
       spw=spw)

# Make moment map

# immoments(imagename=myImage,moments=0,chans='40~90',includepix=[2*rms,100],outfile=myImage+'.mom0')

# go back to where we started
os.chdir(origDir)




