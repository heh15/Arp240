'''
Oct 19th, 2018

to split the 13CO10 data to spw 0,1 and spw 2,3
'''

############################################################
# basic setting

Dir='/home/heh15/workingspace/Arp240/continuum/'
imageDir=Dir+'image/'

vis='/home/heh15/Data/Arp240/arp240-110GHz.ms/'
imagename='NGC5257_13CO10_cont_spw01.image'

tclean(vis=vis,
      imagename=imagename,
      field='0',
      specmode='mfs',
      spw='0:0~400;700~900,1',
      outframe='BARY',
      cell='0.3arcsec',
      imsize=[320,320],
      weighting='briggs',
      robust=0.5,
      deconvolver='hogbom',
      gridder='mosaic',
      niter=100000,
      cyclefactor=1.0,
      pblimit=0.2,
      interactive=True,
      chanchunks=-1,
      threshold='0.1mJy/beam',
      mask=imageDir+'NGC5257_13CO.mask')

# spw=2,3
imagename='NGC5257_13CO10_cont_spw23.image'
tclean(vis=vis,
      imagename=imagename,
      field='0',
      specmode='mfs',
      spw='2,3',
      outframe='BARY',
      cell='0.3arcsec',
      imsize=[320,320],
      weighting='briggs',
      robust=0.5,
      deconvolver='hogbom',
      gridder='mosaic',
      niter=100000,
      cyclefactor=1.0,
      pblimit=0.2,
      interactive=True,
      chanchunks=-1,
      threshold='0.1mJy/beam',
      mask=imageDir+'NGC5257_13CO.mask')
