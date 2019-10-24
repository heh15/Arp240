vis = '/home/heh15/practice/Arp240/' \
      'workingspace/CO13_10/calibrated/NGC5257_CO13_10_12m.ms.cont/'
cleanDir = '/home/heh15/practice/Arp240/'\
             'workingspace/CO13_10/continuum/'
maskDir = '/home/heh15/practice/Arp240/'\
          'workingspace/continuum/'
preName = cleanDir + 'NGC5257_12m_13CO_spw0'

delmod(vis=vis)

tclean(vis=vis,
      imagename=preName,
      specmode='mfs',
      spw='0:0~400;600~900',
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
      mask=maskDir+'NGC5257_13CO.mask/',
      threshold='0.047mJy/beam')
