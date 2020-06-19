imcollapse(imagename='NGC5257_12CO10_combine_contsub.pb/',
           axes=3,
           function="mean",
           outfile='NGC5257_12CO10_combine_mom0.pb')

impbcor(imagename='NGC5257_12CO10_combine_contsub.mom0/',
        pbimage='NGC5257_12CO10_combine_mom0.pb',
        outfile='NGC5257_12CO10_combine_pbcor.mom0')

exportfits(imagename='NGC5257_12CO10_combine_pbcor.mom0/',fitsimage='NGC5257_12CO10_combine_contsub_pbcor_mom0.fits')

# make the moment 8 map. 

immoments(imagename='NGC5257_12CO10_combine_contsub.image',moments=8,outfile='NGC5257_12CO10_combine_contsub.mom8')

exportfits(imagename='NGC5257_12CO10_combine_contsub.mom8',fitsimage='NGC5257_12CO10_combine_contsub_mom8.fits')


# correct for the cube for the primary beam. 

impbcor(imagename='NGC5257_12CO10_combine_contsub.image',pbimage='NGC5257_12CO10_combine_contsub.pb',outfile='NGC5257_12CO10_combine_contsub_pbcor.image')
exportfits(imagename='NGC5257_12CO10_combine_contsub_pbcor.image',fitsimage='NGC5257_12CO10_combine_contsub_pbcor.fits')

# moments map with 2rms cutoff
immoments(imagename='NGC5257_12CO10_combine_contsub.image',moments=0,chans='10~50',includepix=[2*1.6e-3,100],outfile='NGC5257_12CO10_combine_contsub_2rms.mom0')

impbcor(imagename='NGC5257_12CO10_combine_contsub_2rms.mom0/',
        pbimage='NGC5257_12CO10_combine_mom0.pb',
        outfile='NGC5257_12CO10_combine_pbcor_2rms.mom0')

exportfits(imagename='NGC5257_12CO10_combine_pbcor_2rms.mom0', fitsimage='NGC5257_12CO10_combine_pbcor_2rms_mom0.fits')
