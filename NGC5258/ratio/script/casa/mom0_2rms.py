import os
from subprocess import call

Dir='/1/home/heh15/workingspace/Arp240/NGC5258/ratio/'
imageDir=Dir+'image/'
ratioDir=imageDir+'ratio/2110/'

origDir=os.getcwd()

if os.path.isdir('temp')==False:
    os.mkdir('temp')

os.chdir('temp')
rmtables('*')

### 12CO 10 image. 
imsmooth(imagename=imageDir+'12CO10/NGC5258_12CO10_combine_contsub_uvrange.image',
         major='2.186arcsec',
         minor='1.896arcsec',
         pa='-87.31deg',
         targetres=True, 
         outfile='NGC5258_12CO10_combine_contsub_uvrange_smooth.image')

immoments(imagename='NGC5258_12CO10_combine_contsub_uvrange_smooth.image', 
          includepix=[2*1.65e-3, 100], 
          chans='10~60',
          outfile='NGC5258_12CO10_combine_contsub_uvrange_smooth.mom0')

impbcor(imagename='NGC5258_12CO10_combine_contsub_uvrange_smooth.mom0', 
        pbimage=imageDir+'12CO10/NGC5258_12CO10_combine_contsub_uvrange_mom0.pb',
        outfile='NGC5258_12CO10_combine_contsub_uvrange_smooth_pbcor.mom0')

imagename='NGC5258_12CO10_combine_contsub_uvrange_smooth_pbcor.mom0'
rmtables(imageDir+'12CO10/'+imagename)
call(['mv', imagename, imageDir+'12CO10/'])

### 12CO 21 image. 

imsmooth(imagename=imageDir+'12CO21/NGC5258_12CO21_combine_uvtaper.image',
         major='2.186arcsec',
         minor='1.896arcsec',
         pa='-87.31deg',
         targetres=True,
         outfile='NGC5258_12CO21_combine_contsub_uvtaper_smooth.image')

immoments(imagename='NGC5258_12CO21_combine_contsub_uvtaper_smooth.image/', 
          chans='10~60',
          includepix=[2*5.4e-3, 100], 
          outfile='NGC5258_12CO21_combine_contsub_uvtaper_smooth.mom0/')
impbcor(imagename='NGC5258_12CO21_combine_contsub_uvtaper_smooth.mom0/',
        pbimage=imageDir+'12CO21/NGC5258_12CO21_combine_uvtaper_mom0.pb',
        outfile='NGC5258_12CO21_combine_contsub_uvtaper_smooth_pbcor.mom0')

imagename='NGC5258_12CO21_combine_contsub_uvtaper_smooth_pbcor.mom0'
rmtables(imageDir+'12CO21/'+imagename)
call(['mv', imagename, imageDir+'12CO21/'])

### 13CO10 image 
imsmooth(imagename=imageDir+'13CO10/NGC5258_13CO10_12m_uvrange.image',
         major='2.186arcsec',
         minor='1.896arcsec',
         pa='-87.31deg',
         targetres=True, 
         outfile='NGC5258_13CO10_12m_uvrange_smooth.image')

immoments(imagename='NGC5258_13CO10_12m_uvrange_smooth.image', 
          chans='5~30',
          includepix=[2*6.4e-4, 100],
          outfile='NGC5258_13CO10_12m_uvrange_smooth.mom0')

impbcor(imagename='NGC5258_13CO10_12m_uvrange_smooth.mom0',
        pbimage=imageDir+'13CO10/NGC5258_13CO10_12m_uvrange_mom0.pb',
        outfile='NGC5258_13CO10_12m_uvrange_smooth_pbcor.mom0')

imagename='NGC5258_13CO10_12m_uvrange_smooth_pbcor.mom0'
rmtables(imageDir+'13CO10/'+imagename)
call(['mv', imagename, imageDir+'13CO10/'])

os.chdir(origDir)
call(['rm', '-rf', 'temp'])
