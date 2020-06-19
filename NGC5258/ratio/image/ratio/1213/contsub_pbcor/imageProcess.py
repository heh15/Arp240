# exportfits the 13CO 1-0 image. 

exportfits(imagename='NGC5258_13CO10_12m_uvrange_smooth_masked_pbcor.image.mom0', fitsimage='NGC5258_13CO10_12m_uvrange_smooth_masked_pbcor_mom0.fits', overwrite=True)

exportfits(imagename='NGC5258_12CO10_combine_contsub_uvrange_smooth_masked_pbcor.image.mom0', fitsimage='NGC5258_12CO10_combine_contsub_uvrange_smooth_masked_pbcor_mom0.fits', overwrite=True)

immoments(imagename='NGC5258_12CO10_combine_contsub_uvrange_smooth.image', includepix=[4*1.6e-3, 100], chans='10~60', moments=2, outfile='NGC5258_12CO10_combine_contsub_uvrange_smooth.mom2')

exportfits(imagename='NGC5258_12CO10_combine_contsub_uvrange_smooth.mom2', fitsimage='NGC5258_12CO10_combine_contsub_uvrange_smooth_mom2.fits', overwrite=True)
