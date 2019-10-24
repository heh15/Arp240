'''Apr 2nd, 2018
test casa version is 5.1.1

Aug. 22nd, 2019
rms_12CO21=0.0038
'''

import numpy as np
import os
import scipy.ndimage as sni
import sys


# # change directory to workaround makemask failing
# # when specifying paths
# origDir = os.getcwd()
# os.chdir(cleanDir)


# #check the uvrange of 12CO21
# plotms(vis=vis_12CO21,xaxis="uvwave",spw="0:300,1:300,2:300,3:300")

# #contiuous from 6k to 300k, up to 500k.

# # check the uvrange of 12CO10
# plotms(vis=vis_12CO10,xaxis="uvwave",spw="0:600,1:600")
# # continuous from 3k to 140, up to 160k.

# clean the 12CO10 data
vis = '/home/heh15/workingspace/Arp240/NGC5257/12CO10/'\
      'calibrated/NGC5257_combine_CO.ms.contsub'
cleanDir = '/home/heh15/workingspace/Arp240/'\
             'NGC5257/uvtaper/test/'
preName = cleanDir + 'NGC5257_12CO10_combine_contsub_uvrange'
phasecenter='J2000 13h39m52.922 0d50m24.1'
mode = 'velocity'
restfreq='112.73GHz'
width = '10km/s' 
nchan = 70 
start = '-300km/s' 
cell='0.3arcsec'  
imsize = 320
weighting = 'briggs'
robust = 0.5
imagermode = 'mosaic'
outframe='BARY'
cyclefactor = 1.0  # default value
# stuff related to auto-masking
rmsCRTF = '/1/home/heh15/workingspace/Arp240/NGC5257/' \
          '12CO10/region/NGC5257_emission_free_54.crtf'
rmsChan = '60~69'


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

rms = imstat(imagename=preName+'.image',region=rmsCRTF, chans=rmsChan)['rms'][0]
thresh=2*rms

sidelobethreshold=3.0
noisethreshold=4.5
lownoisethreshold=1.5
negativethreshold=0.0
minbeamfrac=0.3
niter=100000

delmod(vis=vis)
tclean(vis=vis,
       imagename=preName,
       phasecenter=phasecenter,
       specmode='cube',
       outframe='Bary',
       restfreq=restfreq,
       width=width,
       nchan=nchan,
       start=start,
       deconvolver='hogbom',
       cell=cell,
       imsize=imsize,
       weighting=weighting,
       robust=robust,
       gridder='mosaic',
       usemask='auto-multithresh',
       niter=100000,
       threshold=str(thresh)+'Jy/beam',
       sidelobethreshold=sidelobethreshold,
       noisethreshold=noisethreshold,
       lownoisethreshold=lownoisethreshold,
       negativethreshold=negativethreshold,       
       interactive=False,
       uvrange='>6klambda')

# get the taper parameters

target_beam = ['1.93arcsec', '1.545arcsec', '-69.762deg']
orig_beam = ['0.9984899163246155arcsec', '0.5498111248016357arcsec', '-64.42718505859375deg']
taper_beam = ia.beamforconvolvedsize(source=orig_beam,convolved=target_beam)
taper_beam = ['{:.7f}arcsec'.format(taper_beam['major']['value']),
               '{:.7f}arcsec'.format(taper_beam['minor']['value']),
               '{:.4f}deg'.format(taper_beam['pa']['value'])]



# clean the 12CO21 data
vis_12CO21 = '/home/heh15/workingspace/Arp240/NGC5257/'\
      '12CO21/calibrated/NGC5257_CO12_21_combine.ms.contsub/'
cleanDir = '/home/heh15/workingspace/Arp240/NGC5257/'\
             'uvtaper/test/'
preName_12CO21 = cleanDir + 'NGC5257_12CO21_combine_contsub_uvtaper'
rmsCRTF = '/home/heh15/workingspace/Arp240/NGC5257/12CO21/'\
          'emission_free.crtf'
rmsChan = '60~69'
pblimit=0.2
phasecenter='J2000 13h39m52.922 0d50m24.1'

# default
mode = 'velocity'
restfreq='225.46GHz'
width = '10km/s' 
nchan = 70 
start = '-300km/s' 
cell='0.3arcsec'  
imsize = 320
weighting = 'briggs'
robust = 0.5
imagermode = 'mosaic'
cyclefactor = 1.0  # default value
specmode='cube'
outframe='BARY'
deconvolver='hogbom' 
gridder='mosaic'

delmod(vis=vis_12CO21)

tclean(vis=vis_12CO21,
      imagename=preName_12CO21,
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
      uvrange='>6klambda',
      uvtaper=taper_beam)

rms = imstat(imagename=preName_12CO21+'.image',
             region=rmsCRTF, chans=rmsChan)['rms'][0]
thresh=2*rms

delmod(vis=vis_12CO21)

# automultithresh parameters: use the default settings in finalImage
sidelobethreshold=3.0
noisethreshold=5.0
lownoisethreshold=1.5
negativethreshold=0.0
minbeamfrac=0.3
niter=100000

tclean(vis=vis_12CO21,
       imagename=preName_12CO21,
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
       interactive=False,
       uvrange='>6klambda',
       uvtaper=taper_beam)






