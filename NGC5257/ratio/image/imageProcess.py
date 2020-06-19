Dir='/1/home/heh15/workingspace/Arp240/NGC5257/ratio/'
regionDir=Dir+'region/'
imageDir=Dir+'image/'

############################################################
# 12CO10 processing
# primary beam correction. 

os.chdir('12CO10')
impbcor(imagename='NGC5257_12CO10_combine_contsub_smooth_2rms.mom0/',
        pbimage='NGC5257_12CO10_combine_mom0.pb/',
        outfile='NGC5257_12CO10_combine_contsub_smooth_2rms_pbcor.mom0')
exportfits(imagename='NGC5257_12CO10_combine_contsub_smooth_2rms_pbcor.mom0',
           fitsimage='NGC5257_12CO10_combine_contsub_smooth_2rms_pbcor_mom0.fits')

# try the rms to be 11.45*13CO10rms

immoments(imagename='NGC5257_12CO10_combine_contsub_smooth.image/',
          moments=0,
          includepix=[11.45*0.00062,100],
          outfile='NGC5257_12CO10_combine_contsub_smooth_13COcut.mom0')

impbcor(imagename='NGC5257_12CO10_combine_contsub_smooth_13COcut.mom0',
        pbimage='NGC5257_12CO10_combine_mom0.pb/',
        outfile='NGC5257_12CO10_combine_contsub_smooth_13COcut_pbcor.mom0')

exportfits(imagename='NGC5257_12CO10_combine_contsub_smooth_13COcut_pbcor.mom0',
           fitsimage='NGC5257_12CO10_combine_contsub_smooth_13COcut_pbcor_mom0.fits')

# smooth the 12CO10 to be CO 3-2 resolution.
rm -rf NGC5257_12CO10_combine_contsub_smooth_13COcut_pbcor_co32.mom0

imsmooth(imagename='NGC5257_12CO10_combine_contsub_smooth_13COcut_pbcor.mom0',
         major='3.789arcsec',
         minor='2.989arcsec',
         pa='-17.357deg',
         targetres=True,
         outfile='NGC5257_12CO10_combine_contsub_smooth_13COcut_pbcor_co32.mom0')

rm -rf NGC5257_12CO10_combine_contsub_smooth_13COcut_pbcor_co32_mom0.fits
exportfits(imagename='NGC5257_12CO10_combine_contsub_smooth_13COcut_pbcor_co32.mom0',fitsimage='NGC5257_12CO10_combine_contsub_smooth_13COcut_pbcor_co32_mom0.fits')

# 12CO10 uvrange moment 0 map to fitsfile. 
exportfits(imagename='NGC5257_12CO10_combine_contsub_uvrange_smooth_pbcor.mom0', fitsimage='NGC5257_12CO10_combine_contsub_uvrange_smooth_pbcor_mom0.fits')

## 12CO10 uvrange map smooth to the 12CO 3-2 resolution and 2rms cut for moment 0 map. 
rmtables('NGC5257_12CO10_combine_contsub_uvrange_smooth_co32.image')
imsmooth(imagename='NGC5257_12CO10_combine_contsub_uvrange.image',
         major='3.789arcsec',
         minor='2.989arcsec',
         pa='-17.357deg',
         targetres=True,
         outfile='NGC5257_12CO10_combine_contsub_uvrange_smooth_co32.image')

imstat(imagename='NGC5257_12CO10_combine_contsub_uvrange_smooth_co32.image', region=regionDir+'emission_free.crtf')['rms'][0]

rmtables('NGC5257_12CO10_combine_contsub_uvrange_smooth_co32.mom0')
immoments(imagename='NGC5257_12CO10_combine_contsub_uvrange_smooth_co32.image', includepix=[2*2.3e-3, 100], moments=0, outfile='NGC5257_12CO10_combine_contsub_uvrange_smooth_co32.mom0')

rmtables('NGC5257_12CO10_combine_contsub_uvrange_smooth_co32_pbcor.mom0')
impbcor(imagename='NGC5257_12CO10_combine_contsub_uvrange_smooth_co32.mom0', pbimage='NGC5257_12CO10_combine_uvrange_mom0.pb', outfile='NGC5257_12CO10_combine_contsub_uvrange_smooth_co32_pbcor.mom0')

exportfits(imagename='NGC5257_12CO10_combine_contsub_uvrange_smooth_co32_pbcor.mom0', fitsimage='NGC5257_12CO10_combine_contsub_uvrange_smooth_co32_pbcor_mom0.fits', overwrite=True)

# make the moment 8 map. 
immoments(imagename='NGC5257_12CO10_combine_contsub_smooth.image', moments=8, includepix=[2*1.6e-3, 100], outfile='NGC5257_12CO10_combine_contsub_smooth.mom8')
impbcor(imagename='NGC5257_12CO10_combine_contsub_smooth.mom8', pbimage='NGC5257_12CO10_combine_mom0.pb/', outfile='NGC5257_12CO10_combine_contsub_smooth_pbcor.mom8')
exportfits(imagename='NGC5257_12CO10_combine_contsub_smooth_pbcor.mom8', fitsimage='NGC5257_12CO10_combine_contsub_smooth_pbcor_mom8.fits')

############################################################
# 13CO10 processing

os.chdir('..')
os.chdir('13CO10')
impbcor(imagename='NGC5257_13CO10_12m_contsub_smooth.mom0/',
        pbimage='NGC5257_13CO10_combine_mom0.pb/',
        outfile='NGC5257_13CO10_12m_contsub_smooth_pbcor.mom0')

# 2rms*width*sqrt(chans)
execfile(scriptDir+'casa/CO13_2rms.py')

exportfits(imagename='NGC5257_13CO10_12m_contsub_smooth_pbcor.mom0',
           fitsimage='NGC5257_13CO10_12m_contsub_smooth_pbcor_mom0.fits')

# smooth the 13CO10 to be CO 3-2 resolution.
rm -rf NGC5257_13CO10_12m_contsub_smooth_pbcor_co32.mom0

imsmooth(imagename='NGC5257_13CO10_12m_contsub_smooth_pbcor.mom0',
         major='3.789arcsec',
         minor='2.989arcsec',
         pa='-17.357deg',
         targetres=True,
         outfile='NGC5257_13CO10_12m_contsub_smooth_pbcor_co32.mom0')

exportfits(imagename='NGC5257_13CO10_12m_contsub_smooth_pbcor_co32.mom0',fitsimage='NGC5257_13CO10_12m_contsub_smooth_pbcor_co32_mom0.fits',overwrite=True)

# first smooth and then make the moment 0 maps. 
imsmooth(imagename='NGC5257_13CO10_12m_contsub_pbcor.image',
         major='3.789arcsec',
         minor='2.989arcsec',
         pa='-17.357deg',
         targetres=True,
         outfile='NGC5257_13CO10_12m_contsub_pbcor_smooth_co32.image')
immoments(imagename='NGC5257_13CO10_12m_contsub_pbcor_smooth_co32.image',
          moments=0,
          includepix=[2*6.4e-4,100],
          chans='5~25',
          outfile='NGC5257_13CO10_12m_contsub_pbcor_smooth_co32_image.mom0')


# export the uvrange smoothed 13CO10 cube. 
exportfits(imagename='NGC5257_13CO10_12m_uvrange_smooth_pbcor.mom0', fitsimage='NGC5257_13CO10_12m_uvrange_smooth_pbcor_mom0.fits')

## 13CO10 uvrange map smooth to the 13CO 3-2 resolution and 2rms cut for moment 0 map. 
rmtables('NGC5257_13CO10_12m_uvrange_smooth_co32.image')
imsmooth(imagename='NGC5257_13CO10_12m_uvrange.image',
         major='3.789arcsec',
         minor='2.989arcsec',
         pa='-17.357deg',
         targetres=True,
         outfile='NGC5257_13CO10_12m_uvrange_smooth_co32.image')

imstat(imagename='NGC5257_13CO10_12m_uvrange_smooth_co32.image', region=regionDir+'emission_free.crtf')['rms'][0]

rmtables('NGC5257_13CO10_12m_uvrange_smooth_co32.mom0')
immoments(imagename='NGC5257_13CO10_12m_uvrange_smooth_co32.image', includepix=[2*7.3e-4, 100], chans='5~30', moments=0, outfile='NGC5257_13CO10_12m_uvrange_smooth_co32.mom0')

rmtables('NGC5257_13CO10_12m_uvrange_smooth_co32_pbcor.mom0')
impbcor(imagename='NGC5257_13CO10_12m_uvrange_smooth_co32.mom0', pbimage='NGC5257_13CO10_12m_uvrange_mom0.pb', outfile='NGC5257_13CO10_12m_uvrange_smooth_co32_pbcor.mom0')

exportfits(imagename='NGC5257_13CO10_12m_uvrange_smooth_co32_pbcor.mom0', fitsimage='NGC5257_13CO10_12m_uvrange_smooth_co32_pbcor_mom0.fits', overwrite=True)

############################################################
# 12CO21 processing 

os.chdir('..')
os.chdir('12CO21')

imsmooth(imagename='NGC5257_12CO21_combine_contsub_uvtaper.image', 
         major='2.186arcsec',
         minor='1.896arcsec',
         pa='-87.31deg',
         targetres=True, 
         outfile='NGC5257_12CO21_combine_contsub_uvtaper_smooth.image')



# try the rms to be 2.98*11.45*13CO10rms

immoments(imagename='NGC5257_12CO21_combine_contsub_uvtaper_smooth_regrid.image/',
          moments=0,
          includepix=[2.98*11.45*0.00062,100],
          outfile='NGC5257_12CO21_combine_contsub_uvtaper_smooth_regrid_13COcut.mom0')

impbcor(imagename='NGC5257_12CO21_combine_contsub_uvtaper_smooth_regrid_13COcut.mom0',
        pbimage='NGC5257_12CO21_combine_uvtaper_mom0_regrid.pb/',
        overwrite=True,
        outfile='NGC5257_12CO21_combine_contsub_uvtaper_smooth_regrid_13COcut_pbcor.mom0')

exportfits(imagename='NGC5257_12CO21_combine_contsub_uvtaper_smooth_regrid_13COcut_pbcor.mom0',
           overwrite=True,
           fitsimage='NGC5257_12CO21_combine_contsub_uvtaper_smooth_regrid_13COcut_pbcor_mom0.fits')

# smooth the 12CO10 to be CO 3-2 resolution.
rm -rf NGC5257_12CO21_combine_contsub_uvtaper_smooth_regrid_13COcut_pbcor_co32.mom0

imsmooth(imagename='NGC5257_12CO21_combine_contsub_uvtaper_smooth_regrid_13COcut_pbcor.mom0',
         major='3.789arcsec',
         minor='2.989arcsec',
         pa='-17.357deg',
         targetres=True,
         outfile='NGC5257_12CO21_combine_contsub_uvtaper_smooth_regrid_13COcut_pbcor_co32.mom0')

exportfits(imagename='NGC5257_12CO21_combine_contsub_uvtaper_smooth_regrid_13COcut_pbcor_co32.mom0',fitsimage='NGC5257_12CO21_combine_contsub_uvtaper_smooth_regrid_13COcut_pbcor_co32_mom0.fits',overwrite=True)

# export the 2rms cut uvtaper 12CO 2-1 moment0 pbcor map. 
exportfits(imagename='NGC5257_12CO21_combine_contsub_uvtaper_smooth_pbcor.mom0', fitsimage='NGC5257_12CO21_combine_contsub_uvtaper_smooth_pbcor_mom0.fits')

## smooth the 12CO 2-1 cube to 12CO 3-2 resolution and make the moment 0 maps. 

rmtables('NGC5257_12CO21_combine_contsub_uvtaper_smooth_co32.image')
imsmooth(imagename='NGC5257_12CO21_combine_contsub_uvtaper.image',
         major='3.789arcsec',
         minor='2.989arcsec',
         pa='-17.357deg',
         targetres=True,
         outfile='NGC5257_12CO21_combine_contsub_uvtaper_smooth_co32.image')
imstat(imagename='NGC5257_12CO21_combine_contsub_uvtaper_smooth_co32.image', region=regionDir+'emission_free.crtf')['rms'][0]

rmtables('NGC5257_12CO21_combine_contsub_uvtaper_smooth_co32.mom0')
immoments(imagename='NGC5257_12CO21_combine_contsub_uvtaper_smooth_co32.image', includepix=[2*9.8e-3, 100], chans='10~60', moments=0, outfile='NGC5257_12CO21_combine_contsub_uvtaper_smooth_co32.mom0')

rmtables('NGC5257_12CO21_combine_contsub_uvtaper_smooth_co32_pbcor.mom0')
impbcor(imagename='NGC5257_12CO21_combine_contsub_uvtaper_smooth_co32.mom0', pbimage='NGC5257_12CO21_combine_uvtaper_mom0.pb', outfile='NGC5257_12CO21_combine_contsub_uvtaper_smooth_co32_pbcor.mom0')

exportfits(imagename='NGC5257_12CO21_combine_contsub_uvtaper_smooth_co32_pbcor.mom0', fitsimage='NGC5257_12CO21_combine_contsub_uvtaper_smooth_co32_pbcor_mom0.fits', overwrite=True)

## make the moment 2 map of the 12CO 2-1 smoothed co32 cube. 
immoments(imagename='NGC5257_12CO21_combine_contsub_uvtaper_smooth_co32.image', includepix=[4*9.8e-3, 100], chans='10~60', moments=2, outfile='NGC5257_12CO21_combine_contsub_uvtaper_smooth_co32.mom2')

# exportfits the 12CO 2-1 uvtapered cube
exportfits(imagename='NGC5257_12CO21_combine_contsub_uvtaper.image', fitsimage='NGC5257_12CO21_combine_contsub_uvtaper.fits')

############################################################ 
# 33GHz image

imsmooth(imagename='ngc5257_Ka_c_r0.5_ms.pbcor.fits',
         outfile='NGC5257_33GHz_pbcor_smooth.image',
         major='2.186arcsec',
         minor='1.896arcsec',
         pa='-87.31deg',
         targetres=True)

imregrid(imagename='NGC5257_33GHz_pbcor_smooth.image',template='../12CO10/NGC5257_12CO10_combine_contsub.image/',output='NGC5257_33GHz_smooth_regrid.image')

imsmooth(imagename='ngc5257_Ka_c_r0.5_ms.fits',
         outfile='NGC5257_33GHz_smooth.image',
         major='2.186arcsec',
         minor='1.896arcsec',
         pa='-87.31deg',
         targetres=True)

imregrid(imagename='NGC5257_33GHz_smooth.image',template='../12CO10/NGC5257_12CO10_combine_contsub.image/',output='NGC5257_33GHz_smooth_regrid.image')


imsmooth(imagename='NGC5257_33GHz_smooth_regrid.image',
         outfile='NGC5257_33GHz_smooth_regrid_co32.image',
         major='3.789arcsec',
         minor='2.989arcsec',
         pa='-17.357deg',
         targetres=True,
         overwrite=True)

exportfits(imagename='NGC5257_33GHz_smooth_regrid_co32.image',fitsimage='NGC5257_33GHz_smooth_regrid_co32.fits')

############################################################
# 12CO 3-2 processing
os.chdir('12CO32')

# regrid the 12CO 3-2 image to the same as 12CO10 image
imregrid(imagename='NGC5257co32_all_map40r_shift.fits', template=imageDir+'12CO10/NGC5257_12CO10_combine_contsub_uvrange_smooth_co32.mom0/', output='NGC5257co32_all_map40r_shift_regrid.mom0')
exportfits(imagename='NGC5257co32_all_map40r_shift_regrid.mom0', fitsimage='NGC5257co32_all_map40r_shift_regrid_mom0.fits')

imstat(imagename='NGC5257co32_all_map40r.fits', region=regionDir+'emission_free.crtf')['rms'][0]

rmtables('NGC5257co32_all_map40r_2rms.mom0')
imstat(imagename='NGC5257co32_all_map40r.fits', region=regionDir+'emission_free_center.crtf', chans='0~8')['rms'][0]
immoments(imagename='NGC5257co32_all_map40r.fits', includepix=[2*0.09,100],chans='10~20', outfile='NGC5257co32_all_map40r_2rms.mom0')
exportfits(imagename='NGC5257co32_all_map40r_2rms.mom0', fitsimage='NGC5257co32_all_map40r_2rms_mom0.fits', overwrite=True)

# import the fitsimage of 12CO 3-2 cube. 
importfits(fitsimage='NGC5257co32_all_map40r.fits', imagename='NGC5257co32_all_map40r.image')
