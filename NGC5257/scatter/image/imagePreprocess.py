## Bin the image

imrebin(imagename='NGC5257_12CO21_combine_noise40.image',factor=[5,5],outfile='NGC5257_12CO21_combine_5bin.image')

# rms=0.0025
imstat(imagename='NGC5257_12CO21_combine_5bin.image/',region='rmsMeasure.crtf')['rms']

immoments(imagename='NGC5257_12CO21_combine_5bin.image/',moments=0,chans='10~60',includepix=[0.0024*2,100],outfile='NGC5257_12CO21_combine_5bin.image.mom0')
exportfits(imagename='NGC5257_12CO21_combine_5bin.image.mom0/',fitsimage='NGC5257_12CO21_combine_5bin_mom0.fits')

immoments(imagename='NGC5257_12CO21_combine_5bin.image/',moments=2,chans='10~60',includepix=[0.0024*4,100],outfile='NGC5257_12CO21_combine_5bin.image.mom2')
exportfits(imagename='NGC5257_12CO21_combine_5bin.image.mom2/',fitsimage='NGC5257_12CO21_combine_5bin_mom2.fits')

