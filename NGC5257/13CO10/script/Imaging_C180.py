vis='/home/heh15/Data/Arp240/arp240-110GHz.ms.contsub'
cleanDir = '/home/heh15/workingspace/Arp240/'\
             'NGC5257/13CO10/test/'
preName = cleanDir + 'NGC5257_C18O'

field = '0'
phasecenter='J2000 13h39m52.922 0d50m24.1'
mode = 'velocity'
restfreq='107.345GHz'
width = '40km/s' 
nchan = 20 
start = '-400km/s' 
cell='0.3arcsec'  
imsize = [320,320]
weighting = 'briggs'
robust = 0.5
imagermode = 'mosaic'
cyclefactor = 1.0  # default value
spw='0'

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
       spw='0')
