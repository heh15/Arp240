''' Create in Dec. 11th, 2019
'''
############################################################
# Directory

Dir='/1/home/heh15/workingspace/Arp240/scatter/'
imageDir=Dir+'image/'
logDir=Dir+'log/'
picDir=Dir+'picture/'
mapDir=Dir+'map/'

############################################################

### regrid the gas component fraction.
## NGC 5258 
imregrid(imagename='NGC5258_gas_vol_fraction.fits', template=imageDir+'NGC5258/NGC5258_12CO21_combine_noise45.image', output='NGC5258_gas_vol_fraction_regrid.image')
exportfits(imagename='NGC5258_gas_vol_fraction_regrid.image', fitsimage='NGC5258_gas_vol_fraction_regrid.fits')

## NGC 5257
imregrid(imagename='NGC5257_gas_vol_fraction.fits', template=imageDir+'NGC5257/12CO21/NGC5257_12CO21_combine_noise40.image.mom0/', output='NGC5257_gas_vol_fraction_regrid.image')
exportfits(imagename='NGC5257_gas_vol_fraction_regrid.image', fitsimage='NGC5257_gas_vol_fraction_regrid.fits')
