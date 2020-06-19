Dir='/1/home/heh15/workingspace/Arp240/NGC5258/ratio/'
regionDir=Dir+'region/'

############################################################

rmtables('NGC5258_33GHz_smooth.image')
imsmooth(imagename='ngc5258_Ka_c_r0.5_ms.fits',
         major='2.186arcsec',
         minor='1.896arcsec',
         pa='-87.6deg',
         outfile='NGC5258_33GHz_smooth.image',
         targetres=True)

imregrid(imagename='NGC5258_33GHz_smooth.image',
         template='../12CO10/NGC5258_12CO10_combine_contsub_pbcor.mom0',
         output='NGC5258_33GHz_smooth_regrid.image')

exportfits(imagename='NGC5258_33GHz_smooth_regrid.image',fitsimage='NGC5258_33GHz_smooth_regrid.fits',overwrite=True)

imsmooth(imagename='ngc5258_Ka_c_r0.5_ms.fits',
         major='3.79arcsec',
         minor='2.99arcsec',
         pa='-17.4deg',
         outfile='NGC5258_33GHz_smooth_co32.image',
         targetres=True)

imregrid(imagename='NGC5258_33GHz_smooth_co32.image',
         template='../12CO10/NGC5258_12CO10_combine_contsub_pbcor.mom0',
         output='NGC5258_33GHz_smooth_co32_regrid.image')


############################################################
# 12CO10 line. 
imsmooth(imagename='12CO10/NGC5258_12CO10_combine_contsub_pbcor.mom0/',
         outfile='12CO10/NGC5258_12CO10_combine_contsub_pbcor_smooth.mom0',
         major='2.186arcsec',
         minor='1.896arcsec',
         pa='-87.6deg',
         targetres=True)

exportfits(imagename='12CO10/NGC5258_12CO10_combine_contsub_pbcor_smooth.mom0/',fitsimage='12CO10/NGC5258_12CO10_combine_contsub_pbcor_smooth_mom0.fits')

imsmooth(imagename='12CO10/NGC5258_12CO10_combine_contsub_pbcor.mom0/',
         major='3.79arcsec',
         minor='2.99arcsec',
         pa='-17.4deg',
         outfile='12CO10/NGC5258_12CO10_combine_contsub_pbcor_co32.mom0',
         targetres=True)

exportfits(imagename='12CO10/NGC5258_12CO10_combine_contsub_pbcor_co32.mom0',fitsimage='12CO10/NGC5258_12CO10_combine_contsub_pbcor_co32_mom0.fits')

# make the moment 2 maps. 
immoments(imagename='NGC5258_12CO10_combine_contsub_pbcor_smooth.image',moments=2, includepix=[4*1.6e-3, 100], outfile='NGC5258_12CO10_combine_contsub_pbcor_smooth.mom2')
exportfits(imagename='NGC5258_12CO10_combine_contsub_pbcor_smooth.mom2', fitsimage='NGC5258_12CO10_combine_contsub_pbcor_smooth_mom2.fits')

## 12CO10 uvrange map smooth to the 12CO 3-2 resolution and 2rms cut for moment 0 map. 
rmtables('NGC5258_12CO10_combine_contsub_uvrange_smooth_co32.image')
imsmooth(imagename='NGC5258_12CO10_combine_contsub_uvrange.image',
         major='3.789arcsec',
         minor='2.989arcsec',
         pa='-17.357deg',
         targetres=True,
         outfile='NGC5258_12CO10_combine_contsub_uvrange_smooth_co32.image')

imstat(imagename='NGC5258_12CO10_combine_contsub_uvrange_smooth_co32.image', region=regionDir+'emission_free.crtf')['rms'][0]

rmtables('NGC5258_12CO10_combine_contsub_uvrange_smooth_co32.mom0')
immoments(imagename='NGC5258_12CO10_combine_contsub_uvrange_smooth_co32.image', includepix=[2*2.3e-3, 100], moments=0, outfile='NGC5258_12CO10_combine_contsub_uvrange_smooth_co32.mom0')

rmtables('NGC5258_12CO10_combine_contsub_uvrange_smooth_co32_pbcor.mom0')
impbcor(imagename='NGC5258_12CO10_combine_contsub_uvrange_smooth_co32.mom0', pbimage='NGC5258_12CO10_combine_contsub_uvrange_mom0.pb', outfile='NGC5258_12CO10_combine_contsub_uvrange_smooth_co32_pbcor.mom0')

exportfits(imagename='NGC5258_12CO10_combine_contsub_uvrange_smooth_co32_pbcor.mom0', fitsimage='NGC5258_12CO10_combine_contsub_uvrange_smooth_co32_pbcor_mom0.fits', overwrite=True)

### make the 12CO10 moment 8 map
imsmooth(imagename='NGC5258_12CO10_combine_contsub_uvrange.image/',
         major='2.186arcsec',
         minor='1.896arcsec',
         pa='-87.6deg',
         targetres=True, 
         outfile='NGC5258_12CO10_combine_contsub_uvrange_smooth.image')
immoments(imagename='NGC5258_12CO10_combine_contsub_uvrange_smooth.image/', moments=8, includepix=[2*1.6e-3, 100], outfile='NGC5258_12CO10_combine_contsub_uvrange_smooth.mom8')
impbcor(imagename='NGC5258_12CO10_combine_contsub_uvrange_smooth.mom8', pbimage='NGC5258_12CO10_combine_contsub_uvrange_mom0.pb/', outfile='NGC5258_12CO10_combine_contsub_uvrange_smooth_pbcor.mom8')
exportfits(imagename='NGC5258_12CO10_combine_contsub_uvrange_smooth_pbcor.mom8', fitsimage='NGC5258_12CO10_combine_contsub_uvrange_smooth_pbcor_mom8.fits')

############################################################
# 13CO10 line. 
imsmooth(imagename='13CO10/NGC5258_13CO10_contsub_pbcor.mom0/',
         outfile='13CO10/NGC5258_13CO10_contsub_pbcor_smooth.mom0',
         major='2.186arcsec',
         minor='1.896arcsec',
         pa='-87.6deg',
         targetres=True)

exportfits(imagename='13CO10/NGC5258_13CO10_contsub_pbcor_smooth.mom0',fitsimage='13CO10/NGC5258_13CO10_contsub_pbcor_smooth_mom0.fits',overwrite=True)

imsmooth(imagename='13CO10/NGC5258_13CO10_contsub_pbcor.mom0/',
         major='3.79arcsec',
         minor='2.99arcsec',
         pa='-17.4deg',
         outfile='13CO10/NGC5258_13CO10_contsub_pbcor_smooth_co32.mom0/',
         targetres=True)

exportfits(imagename='13CO10/NGC5258_13CO10_contsub_pbcor_smooth_co32.mom0/',fitsimage='13CO10/NGC5258_13CO10_contsub_pbcor_smooth_co32_mom0.fits',overwrite=True)

## 13CO10 uvrange map smooth to the 13CO 3-2 resolution and 2rms cut for moment 0 map. 
rmtables('NGC5258_13CO10_12m_uvrange_smooth_co32.image')
imsmooth(imagename='NGC5258_13CO10_12m_uvrange.image',
         major='3.789arcsec',
         minor='2.989arcsec',
         pa='-17.357deg',
         targetres=True,
         outfile='NGC5258_13CO10_12m_uvrange_smooth_co32.image')

imstat(imagename='NGC5258_13CO10_12m_uvrange_smooth_co32.image', region=regionDir+'emission_free.crtf')['rms'][0]

rmtables('NGC5258_13CO10_12m_uvrange_smooth_co32.mom0')
immoments(imagename='NGC5258_13CO10_12m_uvrange_smooth_co32.image', includepix=[2*7.3e-4, 100], chans='5~30', moments=0, outfile='NGC5258_13CO10_12m_uvrange_smooth_co32.mom0')

rmtables('NGC5258_13CO10_12m_uvrange_smooth_co32_pbcor.mom0')
impbcor(imagename='NGC5258_13CO10_12m_uvrange_smooth_co32.mom0', pbimage='NGC5258_13CO10_12m_uvrange_mom0.pb', outfile='NGC5258_13CO10_12m_uvrange_smooth_co32_pbcor.mom0')

exportfits(imagename='NGC5258_13CO10_12m_uvrange_smooth_co32_pbcor.mom0', fitsimage='NGC5258_13CO10_12m_uvrange_smooth_co32_pbcor_mom0.fits', overwrite=True)

############################################################
# 12CO21 line.
imsmooth(imagename='12CO21/NGC5258_12CO21_combine_uvtaper_pbcor.mom0/',
         outfile='12CO21/NGC5258_12CO21_combine_uvtaper_pbcor_smooth.mom0',
         major='2.186arcsec',
         minor='1.896arcsec',
         pa='-87.6deg',
         targetres=True)

exportfits(imagename='12CO21/NGC5258_12CO21_combine_uvtaper_pbcor_smooth.mom0',fitsimage='12CO21/NGC5258_12CO21_combine_uvtaper_pbcor_smooth_mom0.fits')

imsmooth(imagename='12CO21/NGC5258_12CO21_combine_uvtaper_pbcor.mom0/',
         major='3.79arcsec',
         minor='2.99arcsec',
         pa='-17.4deg',
         outfile='12CO21/NGC5258_12CO21_combine_uvtaper_pbcor_smooth_co32.mom0',
         targetres=True)

exportfits(imagename='12CO21/NGC5258_12CO21_combine_uvtaper_pbcor_smooth_co32.mom0',fitsimage='12CO21/NGC5258_12CO21_combine_uvtaper_pbcor_smooth_co32_mom0.fits')

# exportfits of the 12CO21 cube. 
exportfits(imagename='12CO21/NGC5258_12CO21_combine_uvtaper.image',fitsimage='12CO21/NGC5258_12CO21_combine_uvtaper.fits')

## smooth the 12CO 2-1 cube to 12CO 3-2 resolution and make the moment 0 maps. 

rmtables('NGC5258_12CO21_combine_uvtaper_smooth_co32.image')
imsmooth(imagename='NGC5258_12CO21_combine_uvtaper.image',
         major='3.789arcsec',
         minor='2.989arcsec',
         pa='-17.357deg',
         targetres=True,
         outfile='NGC5258_12CO21_combine_uvtaper_smooth_co32.image')
imstat(imagename='NGC5258_12CO21_combine_uvtaper_smooth_co32.image', region=regionDir+'emission_free.crtf')['rms'][0]

rmtables('NGC5258_12CO21_combine_uvtaper_smooth_co32.mom0')
immoments(imagename='NGC5258_12CO21_combine_uvtaper_smooth_co32.image', includepix=[2*9.8e-3, 100], chans='10~60', moments=0, outfile='NGC5258_12CO21_combine_uvtaper_smooth_co32.mom0')

# make the moment 2 map for smoothed cube

immoments(imagename='NGC5258_12CO21_combine_uvtaper_smooth_co32.image', includepix=[4*9.8e-3, 100], chans='10~60', moments=2, outfile='NGC5258_12CO21_combine_uvtaper_smooth_co32.mom2')

rmtables('NGC5258_12CO21_combine_uvtaper_smooth_co32_pbcor.mom0')
impbcor(imagename='NGC5258_12CO21_combine_uvtaper_smooth_co32.mom0', pbimage='NGC5258_12CO21_combine_uvtaper_mom0.pb', outfile='NGC5258_12CO21_combine_uvtaper_smooth_co32_pbcor.mom0')

exportfits(imagename='NGC5258_12CO21_combine_uvtaper_smooth_co32_pbcor.mom0', fitsimage='NGC5258_12CO21_combine_uvtaper_smooth_co32_pbcor_mom0.fits', overwrite=True)


## exportfits the 12CO21 uvtaper smoothed image cube
exportfits(imagename='NGC5258_12CO21_combine_uvtaper_smooth_co32.image', fitsimage='NGC5258_12CO21_combine_uvtaper_smooth_co32.fits')

############################################################
# 12CO3-2 line. 
importfits(fitsimage='12CO32/NGC5258co32_all.map40r.mom0.fits',imagename='12CO32/NGC5258co32_all_regrid.map40r.mom0',overwrite=True)

imregrid(imagename='12CO32/NGC5258co32_all.map40r.mom0.fits',template='12CO10/NGC5258_12CO10_combine_contsub_pbcor.mom0/',output='12CO32/NGC5258co32_all_regrid.map40r.mom0',overwrite=True)

exportfits(imagename='12CO32/NGC5258co32_all_regrid.map40r.mom0',fitsimage='12CO32/NGC5258co32_all_regrid.map40r.mom0.fits',overwrite=True)

# regrid the cube
imregrid(imagename='12CO32/NGC5258co32_all_map40r.fits',template='12CO10/NGC5258_12CO10_combine_contsub_pbcor.mom0/',output='12CO32/NGC5258co32_all_regrid.map40r.image',overwrite=True)

### make the moment 0 map from the cube. 
rmtables('NGC5258co32_all_map40r_2rms.mom0')
imstat(imagename='NGC5258co32_all_map40r.fits', region=regionDir+'emission_free_center.crtf', chans='1~10')['rms'][0]
immoments(imagename='NGC5258co32_all_map40r.fits', moments=0, includepix=[2*0.1, 100],chans='16~21',  outfile='NGC5258co32_all_map40r_2rms.mom0')

exportfits(imagename='NGC5258co32_all_map40r_2rms.mom0', fitsimage='NGC5258co32_all_map40r_2rms_mom0.fits', overwrite=True)

# import image from co32 fits cube
importfits(fitsimage='NGC5258co32_all_map40r.fits', imagename='NGC5258co32_all_map40r.image')
