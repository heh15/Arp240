'''
Oct.18th, 2018

Smooth the VLA image to have the same beam size as the ALMA image

'''

Dir='/1/home/heh15/workingspace/Arp240/continuum/'
workDir=Dir+'image/'

origDir=os.getcwd()

os.chdir(workDir)

imsmooth(imagename='ngc5257_Ka_c_r0.5_ms.fits',
         kernel='gauss',
         major='2.186arcsec',
         minor='1.896arcsec',
         pa='-87.314deg',
         targetres=True,
         outfile='NGC5257_33GHz_smooth.image')

imsmooth(imagename='NGC5257_12CO10_cont_12m.image/',
         kernel='gauss',
         major='2.186arcsec',
         minor='1.896arcsec',
         pa='-87.314deg',
         targetres=True,
         outfile='NGC5257_12CO10_cont_12m_smooth.image')

imsmooth(imagename='NGC5257_12m_cont_12m_13CO.image/',
         kernel='gauss',
         major='2.186arcsec',
         minor='1.896arcsec',
         pa='-87.314deg',
         targetres=True,
         outfile='NGC5257_13CO10_cont_12m_smooth.image')
