import numpy as np
import os
import scipy.ndimage as sni
import sys

############################################################
# directory 
Dir='/home/heh15/workingspace/Arp240/NGC5258/uvtaper/'
regionDir=Dir+'region/'
imageDir=Dir+'image/'


vis='/1/home/heh15/workingspace/Arp240/NGC5258/13CO10/calibrated/NGC5258_CO13_12m_contsub.ms'

# # select one channel so that it saves the time. 
# plotms(vis=vis,xaxis="uvwave", spw='0:480')

### start the tclean 
cleanDir = '/home/heh15/workingspace/Arp240/NGC5258/'\
             'uvtaper/image/'
preName = cleanDir + 'NGC5258_13CO10_12m_uvrange'
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

specmode='cube'
outframe='BARY'
deconvolver='hogbom' 
# gridder='mosaic' 
pblimit=0.2

rmsCRTF = regionDir+'emissionFree_12CO10.crtf'
rmsChan = '29~34'

delmod(vis=vis)
tclean(vis=vis,
      imagename=preName,
      phasecenter=phasecenter,
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
      uvrange='>6klambda')

rms = imstat(imagename=preName+'.image',
             region=rmsCRTF, chans=rmsChan)['rms'][0]
thresh=2*rms

sidelobethreshold=3.0
noisethreshold=4.2
lownoisethreshold=1.5
negativethreshold=0.0
minbeamfrac=0.3
niter=100000

delmod(vis=vis)


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
       interactive=False,
       uvrange='>6klambda')

impbcor(imagename=preName+'.image',
        pbimage=preName+'.pb',
        outfile=preName+'_pbcor.image')
