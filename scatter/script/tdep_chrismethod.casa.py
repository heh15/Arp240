import os

Dir='/home/heh15/workingspace/Arp240/scatter/'
imageDir=Dir+'image/'
NGC5257dir=imageDir+'NGC5257/'
crtDir=os.getcwd()

# NGC 5257 - script for new processing for KS analysis
#---------------------------------------------------
# May 15, 2019

# adding noise map calculations

# things to adjust if applying to other galaxies:
#    factor in imrebin: here 
#           [10,10]
#    region in imsubimage: here
#           region= 'box [ [ 215pix  , 218pix] , [745pix , 748pix ] ]')
#    noise values in first calculation of noise images: here
#           0.045 CO,  1.0E-5 continuum
#    Jy/K factor in derived quantitites: here
#           1/(0.0109*1.1*0.8*(225.46/115.271)**2)=27.25
#    frequency in continuum name: here
#           33
#    ppb factor: here (pixel per beam) 1.1331*1.1*0.8/0.01=99.72
#           99.72
#    this galaxy needs AGN removed from center ***

# NEW: needed for sigma_v uncertainty
# total line width
#    vchan = 10 km/s, Nchan = 50 : here
# --> 500 km/s
# rms=3.1e-3
# sig(mom0): here 
#    = rms*10*sqrt(50)= 0.22 Jy/beam km/s
# sig(mom2) = 0.22*500**2/17.68*mom0/mom2: here
#    = 3111/mom0/mom2               *** CALCULATE ***


# before running script, need to copy files over from previous analysis area
# os.mkdir('temp')
os.chdir('temp')

immoments(imagename=NGC5257dir+'12CO21/NGC5257_12CO21_combine_smooth.image',moments=0, outfile='NGC5257_12CO21_combine_smooth.mom0')
immoments(imagename=NGC5257dir+'12CO21/NGC5257_12CO21_combine_smooth.image',moments=2,includepix=[5*3.1e-3,100],outfile='NGC5257_12CO21_5sig.mom2')
immoments(imagename=NGC5257dir+'12CO21/NGC5257_12CO21_combine_smooth.image',moments=0,includepix=[3*3.1e-3,100],outfile='NGC5257_12CO21_3sig.mom0')
impbcor(imagename='NGC5257_12CO21_3sig.mom0',pbimage=NGC5257dir+'12CO21/NGC5257_12CO21_combine_mom0.pb',outfile='NGC5257_12CO21_3sig.mom0.pbcor')

immath(imagename=NGC5257dir+'33GHz/NGC5257_33GHz_regrid.image',expr='IM0/IM0*3.1e-3',outfile='NGC5257_33GHz.noise')

imtrans(imagename='NGC5257_33GHz.noise',outfile='NGC5257_33GHz_trans.noise',order='0132')
imtrans(imagename=NGC5257dir+'33GHz/NGC5257_33GHz_pbcor_regrid_smooth.image',outfile='NGC5257_33GHz_pbcor_regrid_smooth_trans.image',order='0132')

# cp -rp imageDir+

#cp -rp newAnalysis/ngc7469_CO10_3sig.mom0.pbcor ks_with_noise/
#cp -rp newAnalysis/ngc7469_93GHz_smooth.pbcor ks_with_noise/
#cp -rp newAnalysis/ngc7469_CO10_4sig.mom2 ks_with_noise/

# in casa, change to right directory

# cd ks_with_noise

image1 = 'NGC5257_12CO21_combine_smooth.mom0'
image2 = 'NGC5257_12CO21_3sig.mom0.pbcor'
image3 = 'NGC5257_33GHz_trans.noise'
image4 = 'NGC5257_33GHz_pbcor_regrid_smooth_trans.image'
image5 = 'NGC5257_12CO21_5sig.mom2'

image6 = 'NGC5257_12CO21_sigmasquared4sig_over_Sigmagas'
image7 = 'NGC5257_12CO21_SigmaSFR_over_SigmaGas'
image8 = 'NGC5257_12CO21_sigma4sig_over_Sigmagas'
image9 = 'NGC5257_12CO21_SigmaSFR_over_SigmaGasSquared_times_sigma'

image10 = 'NGC5257_12CO21_4sig.mom2.noise'
image11 = 'NGC5257_12CO21_4sig.mom2.noisesquared'


################################


# 0. first make noise image for velocity dispersion map (small pixels)

# formula for mom2 noise = sig(I)/I (Nchan x vchan)^2/17.68 / mom2
#  (note: missing sqrt(2) now included)

# so sig(mom2) = 3110/mom0/mom2  (see starting info)

immath(imagename=[image2,image5],
      mode='evalexpr',outfile=image10,
      expr=' 3110/IM0/IM1')

################################

# re-organizing script somewhat

# 1. calculating some quantities at original resolution

# now combine some images to get scale height
# NOTE: leaving scalings on physical quantities the same as before


# 3.2 = MW conversion from CO to H2
# other factor is K/Jy for the CO data - see notebook
# K/Jy does not include redshifted frequency factor

immath(imagename=[image5,image2,image4],  # scale height, could ignore for a while. 
      mode='evalexpr',outfile=image6,
      expr=' IM0*IM0*IM2/IM2/IM1/27.25/3.2')


# now, I also need to calculate t_ff, t_gas, epsilon_ff; 
# *** last 2 depend on
# SFR and so need the differential dust correction applied

# t_gas

immath(imagename=[image4,image2],
      mode='evalexpr',outfile=image7,
      expr=' IM0/IM1')

# now calculate t_ff = sigma/Sigma_gas


# 3.2 = MW conversion from CO to H2
# other factor is K/Jy for the CO data - see notebook
# K/Jy = 399 if include redshifted frequency factor

immath(imagename=[image5,image2,image4],
      mode='evalexpr',outfile=image8,
      expr=' IM0*IM2/IM2/IM1/27.25/1.088')


# now for epsilon_ff proxy

# 33GHz, CO, sigma_v
# IM0, IM1, IM2

immath(imagename=[image4,image2,image5],
      mode='evalexpr',outfile=image9,
      expr=' IM0*IM2/IM1/IM1/27.25/3.2')



# 2B. make subimage

# example (not for this galaxy/binning):
# playing around with CO mom0 image
# 131,136 is 93 GHz peak; peak of CO is same
# going to bin by 9 so want odd multiple of 9 e.g. 63
# 63-1/2 = 31

imsubimage(imagename=image2,
    outfile=image2+'.smallks',
           region= 'box [ [ 215pix  , 218pix] , [745pix , 748pix ] ]')

imsubimage(imagename=image4,
    outfile=image4+'.smallks',
           region= 'box [ [ 215pix  , 218pix] , [745pix , 748pix ] ]')

imsubimage(imagename=image6,
    outfile=image6+'.smallks',
           region= 'box [ [215pix, 218pix], [745pix, 748pix ] ]')

imsubimage(imagename=image5,
    outfile=image5+'.smallks',
           region= 'box [ [215pix, 218pix], [745pix, 748pix ] ]')

imsubimage(imagename=image7,
    outfile=image7+'.smallks',
           region= 'box [ [215pix, 218pix], [745pix, 748pix ] ]')

imsubimage(imagename=image8,
    outfile=image8+'.smallks',
           region= 'box [ [215pix, 218pix], [745pix, 748pix ] ]')

imsubimage(imagename=image9,
    outfile=image9+'.smallks',
           region= 'box [ [215pix, 218pix], [745pix, 748pix ] ]')

imsubimage(imagename=image10,
    outfile=image10+'.smallks',
           region= 'box [ [215pix, 218pix], [745pix, 748pix ] ]')

# *** maybe check small mom2 to make sure not too many blank pixels

# 2BB. REBIN MOM2 UNCERTAINTY IMAGE

# needs to be squared and unsquared


immath(imagename=[image10+'.smallks'],
      mode='evalexpr',outfile=image11+'.smallks',
      expr=' IM0*IM0')

imrebin(imagename=image11+'.smallks',
    outfile=image11+'.smallks.rebin',
        factor=[10,10])

# note: imrebin calculates sum divided by Npix, not just sum
# then square root leaves me a factor of 1/sqrt(Npix)
# so need to multiply by sqrt(ppb/Npix)

# here:
# sqrt(99.72/[10,10]) = 0.9972      *** CALCULATE ***


immath(imagename=[image11+'.smallks.rebin'],
    outfile=image10+'.smallks.rebin',
    expr=' sqrt(IM0)*0.9972')


# 2C. REBIN OTHER IMAGES

imrebin(imagename=image2+'.smallks',
    outfile=image2+'.smallks.rebin',
    factor=[10,10])

imrebin(imagename=image4+'.smallks',
    outfile=image4+'.smallks.rebin',
    factor=[10,10])

imrebin(imagename=image6+'.smallks',
    outfile=image6+'.smallks.rebin',
    factor=[10,10])

imrebin(imagename=image5+'.smallks',
    outfile=image5+'.smallks.rebin',
    factor=[10,10])

imrebin(imagename=image7+'.smallks',
    outfile=image7+'.smallks.rebin',
    factor=[10,10])

imrebin(imagename=image8+'.smallks',
    outfile=image8+'.smallks.rebin',
    factor=[10,10])

imrebin(imagename=image9+'.smallks',
    outfile=image9+'.smallks.rebin',
    factor=[10,10])

# make noise images at rebinned scale
# sqrt((5%*IMAGE)**2 + sig**2 (ppb/Npix))

# adjust Npix: here [10,10]

immath(imagename=[image2+'.smallks.rebin'],
      mode='evalexpr',outfile=image1+'.smallks.rebin',
      expr='sqrt(0.22*0.22*99.72/100)*IM0/IM0')



immath(imagename=[image4+'.smallks.rebin'],
      mode='evalexpr',outfile=image3+'.smallks.rebin',
      expr='sqrt(1.0E-5*1.0E-5*99.72/100)*IM0/IM0')

# to check, calculate S/N images

immath(imagename=[image2+'.smallks.rebin',image1+'.smallks.rebin'],
      mode='evalexpr',outfile=image2+'.smallks.rebin.S2N',
      expr=' IM0/IM1')


immath(imagename=[image4+'.smallks.rebin',image3+'.smallks.rebin'],
      mode='evalexpr',outfile=image4+'.smallks.rebin.S2N',
      expr=' IM0/IM1')

immath(imagename=[image5+'.smallks.rebin',image10+'.smallks.rebin'],
      mode='evalexpr',outfile=image5+'.smallks.rebin.S2N',
      expr=' IM0/IM1')

# I can make a mask from the continuum and sigma_v S/N images: 
# choose S/N > 4 for both SFR and sigma_v
# (damn I hate masks ...)

immath(imagename=[image5+'.smallks.rebin.S2N',image4+'.smallks.rebin.S2N'],
      mode='evalexpr',outfile='premask1',
       expr=' iif( IM0 >=4.0, IM0/IM1, 0.0)')

immath(imagename=[image5+'.smallks.rebin.S2N',image4+'.smallks.rebin.S2N','premask1'],
      mode='evalexpr',outfile='premask2',
       expr=' iif( IM1 >=3.0, IM2/IM0*IM1, 0.0)')

ia.open('premask2')
ia.calcmask(mask='premask2 > 0.0', name='mymask1')
ia.done()
ia.close()


immath(imagename='premask2',
      mode='evalexpr',outfile='maskS2Neq4',
      expr=' IM0')

# calculate noise images for derived quantities: H, t_dep, epsilon, t_ff
# easiest to calculate S/N first and then divide image by S/N to get noise

# reminder of variables:  

#image2 = 'NGC5257_12CO21_3sig.mom0.pbcor'
#image1 = 'NGC5257_12CO21_noise.mom0.pbcor'

#image4 = 'NGC5257_93GHz_smooth.pbcor'
#image3 = 'NGC5257_93GHz_noise.pbcor'

#image5 = 'NGC5257_12CO21_4sig.mom2'
#image10 = 'NGC5257_12CO21_4sig.mom2.noise'

#image6 = 'NGC5257_12CO21_sigmasquared4sig_over_Sigmagas'
#image7 = 'NGC5257_12CO21_SigmaSFR_over_SigmaGas'
#image8 = 'NGC5257_12CO21_sigma4sig_over_Sigmagas'
#image9 = 'NGC5257_12CO21_SigmaSFR_over_SigmaGasSquared_times_sigma'


# S/N images have names like  image2+'.smallks.rebin.S2N'

# image7 = 1/t_dep = SFR/Gas: image2, image4

immath(imagename=[image2+'.smallks.rebin.S2N',image4+'.smallks.rebin.S2N'],
      mode='evalexpr',outfile=image7+'.smallks.rebin.S2N',
      expr=' 1.0/sqrt(1.0/IM0/IM0+1.0/IM1/IM1)')

immath(imagename=[image7+'.smallks.rebin',image7+'.smallks.rebin.S2N'],
      mode='evalexpr',outfile=image7+'.noise.smallks.rebin',
      expr=' IM0/IM1')

 # image8 = t_ff = sigma_v / gas:  image5, image2

immath(imagename=[image2+'.smallks.rebin.S2N',image5+'.smallks.rebin.S2N'],
      mode='evalexpr',outfile=image8+'.smallks.rebin.S2N',
      expr=' 1.0/sqrt(1.0/IM0/IM0+1.0/IM1/IM1)')

immath(imagename=[image8+'.smallks.rebin',image8+'.smallks.rebin.S2N'],
      mode='evalexpr',outfile=image8+'.noise.smallks.rebin',
      expr=' IM0/IM1')


# image6 = H = sigma_v **2 / gas:  image5, image2
# uncertainties in sigma_v correlate - factor of 4

immath(imagename=[image2+'.smallks.rebin.S2N',image5+'.smallks.rebin.S2N'],
      mode='evalexpr',outfile=image6+'.smallks.rebin.S2N',
      expr=' 1.0/sqrt(1.0/IM0/IM0+4.0/IM1/IM1)')

immath(imagename=[image6+'.smallks.rebin',image6+'.smallks.rebin.S2N'],
      mode='evalexpr',outfile=image6+'.noise.smallks.rebin',
      expr=' IM0/IM1')

# image9 = epsilon = SFR * sigma_v / gas **2 : image4, image5, image2
# uncertainties in gas correlate - factor of 4

immath(imagename=[image4+'.smallks.rebin.S2N',image5+'.smallks.rebin.S2N',image2+'.smallks.rebin.S2N'],
      mode='evalexpr',outfile=image9+'.smallks.rebin.S2N',
       expr=' 1.0/sqrt(1.0/IM0/IM0+1.0/IM1/IM1+4.0/IM2/IM2)')

immath(imagename=[image9+'.smallks.rebin',image9+'.smallks.rebin.S2N'],
      mode='evalexpr',outfile=image9+'.noise.smallks.rebin',
      expr=' IM0/IM1')


# finally apply the mask to all images and noise maps

immath(imagename=[image2+'.smallks.rebin','maskS2Neq4'],
      mode='evalexpr',outfile=image2+'.smallks.rebin.masked',
      expr=' IM0*IM1')

immath(imagename=[image4+'.smallks.rebin','maskS2Neq4'],
      mode='evalexpr',outfile=image4+'.smallks.rebin.masked',
      expr=' IM0*IM1')

immath(imagename=[image1+'.smallks.rebin','maskS2Neq4'],
      mode='evalexpr',outfile=image1+'.smallks.rebin.masked',
      expr=' IM0*IM1')

immath(imagename=[image3+'.smallks.rebin','maskS2Neq4'],
      mode='evalexpr',outfile=image3+'.smallks.rebin.masked',
      expr=' IM0*IM1')

immath(imagename=[image5+'.smallks.rebin','maskS2Neq4'],
      mode='evalexpr',outfile=image5+'.smallks.rebin.masked',
      expr=' IM0*IM1')

immath(imagename=[image10+'.smallks.rebin','maskS2Neq4'],
      mode='evalexpr',outfile=image10+'.smallks.rebin.masked',
      expr=' IM0*IM1')

immath(imagename=[image6+'.smallks.rebin','maskS2Neq4'],
      mode='evalexpr',outfile=image6+'.smallks.rebin.masked',
      expr=' IM0*IM1')

immath(imagename=[image7+'.smallks.rebin','maskS2Neq4'],
      mode='evalexpr',outfile=image7+'.smallks.rebin.masked',
      expr=' IM0*IM1')

immath(imagename=[image8+'.smallks.rebin','maskS2Neq4'],
      mode='evalexpr',outfile=image8+'.smallks.rebin.masked',
      expr=' IM0*IM1')

immath(imagename=[image9+'.smallks.rebin','maskS2Neq4'],
      mode='evalexpr',outfile=image9+'.smallks.rebin.masked',
      expr=' IM0*IM1')

# and noise images

immath(imagename=[image6+'.noise.smallks.rebin','maskS2Neq4'],
      mode='evalexpr',outfile=image6+'.noise.smallks.rebin.masked',
      expr=' IM0*IM1')

immath(imagename=[image7+'.noise.smallks.rebin','maskS2Neq4'],
      mode='evalexpr',outfile=image7+'.noise.smallks.rebin.masked',
      expr=' IM0*IM1')

immath(imagename=[image8+'.noise.smallks.rebin','maskS2Neq4'],
      mode='evalexpr',outfile=image8+'.noise.smallks.rebin.masked',
      expr=' IM0*IM1')

immath(imagename=[image9+'.noise.smallks.rebin','maskS2Neq4'],
      mode='evalexpr',outfile=image9+'.noise.smallks.rebin.masked',
      expr=' IM0*IM1')

# 3. write to fits

exportfits(imagename=image2+'.smallks.rebin.masked',
   dropdeg=True,fitsimage=image2+'_ksrebin_masked.fits')

exportfits(imagename=image4+'.smallks.rebin.masked',
   dropdeg=True,fitsimage=image4+'_ksrebin_masked.fits')

exportfits(imagename=image1+'.smallks.rebin.masked',
   dropdeg=True,fitsimage=image1+'_ksrebin_masked.fits')

exportfits(imagename=image3+'.smallks.rebin.masked',
   dropdeg=True,fitsimage=image3+'_ksrebin_masked.fits')

exportfits(imagename=image5+'.smallks.rebin.masked',
   dropdeg=True,fitsimage=image5+'_ksrebin_masked.fits')

exportfits(imagename=image10+'.smallks.rebin.masked',
   dropdeg=True,fitsimage=image10+'_ksrebin_masked.fits')

exportfits(imagename=image6+'.smallks.rebin.masked',
   dropdeg=True,fitsimage=image6+'_ksrebin_masked.fits')

exportfits(imagename=image7+'.smallks.rebin.masked',
   dropdeg=True,fitsimage=image7+'_ksrebin_masked.fits')

exportfits(imagename=image8+'.smallks.rebin.masked',
   dropdeg=True,fitsimage=image8+'_ksrebin_masked.fits')

exportfits(imagename=image9+'.smallks.rebin.masked',
   dropdeg=True,fitsimage=image9+'_ksrebin_masked.fits')


exportfits(imagename=image6+'.noise.smallks.rebin.masked',
   dropdeg=True,fitsimage=image6+'_noise_ksrebin_masked.fits')

exportfits(imagename=image7+'.noise.smallks.rebin.masked',
   dropdeg=True,fitsimage=image7+'_noise_ksrebin_masked.fits')

exportfits(imagename=image8+'.noise.smallks.rebin.masked',
   dropdeg=True,fitsimage=image8+'_noise_ksrebin_masked.fits')

exportfits(imagename=image9+'.noise.smallks.rebin.masked',
   dropdeg=True,fitsimage=image9+'_noise_ksrebin_masked.fits')


# ######################################################
# # removing AGN from center
# # these commands done by hand (not as part of script)
# ######################################################

# cp -rp NGC5257_12CO21_3sig.mom0.pbcor.smallks.rebin.masked NGC5257_12CO21_3sig.mom0.pbcor.agnmask


# ia.open(image2+'.agnmask')
# ia.calcmask(mask=image2+'.agnmask <= 7.0', name='agnmask')
# ia.done()
# ia.close()

# exportfits(imagename=image2+'.agnmask',
#    dropdeg=True,fitsimage=image2+'_ksrebin_masked_noagn.fits')

# # repeat for mom2 (for plotting purposes)

# cp -rp NGC5257_12CO21_4sig.mom2.smallks.rebin.masked NGC5257_12CO21_4sig.mom2.agnmask


# ia.open(image5+'.agnmask')
# ia.calcmask(mask=image5+'.agnmask <= 60.0', name='agnmask')
# ia.done()
# ia.close()

# exportfits(imagename=image5+'.agnmask',
#    dropdeg=True,fitsimage=image5+'_ksrebin_masked_noagn.fits')

# # and for t_dep (for plotting purposes)

# cp -rp NGC5257_12CO21_SigmaSFR_over_SigmaGas.smallks.rebin.masked NGC5257_12CO21_SigmaSFR_over_SigmaGas.agnmask

# ia.open(image7+'.agnmask')
# ia.calcmask(mask=image7+'.agnmask <= 0.0002', name='agnmask')
# ia.done()
# ia.close()

# exportfits(imagename=image7+'.agnmask',
#    dropdeg=True,fitsimage=image7+'_ksrebin_masked_noagn.fits')



# # added: fits of S/N images (helpful for adding errorbars to plots
# # these run by hand

# # need to apply mask to S2N image

# immath(imagename=[image2+'.smallks.rebin.S2N','maskS2Neq4'],
#       mode='evalexpr',outfile=image2+'.smallks.rebin.S2N.masked',
#       expr=' IM0*IM1')

# immath(imagename=[image4+'.smallks.rebin.S2N','maskS2Neq4'],
#       mode='evalexpr',outfile=image4+'.smallks.rebin.S2N.masked',
#       expr=' IM0*IM1')

# immath(imagename=[image5+'.smallks.rebin.S2N','maskS2Neq4'],
#       mode='evalexpr',outfile=image5+'.smallks.rebin.S2N.masked',
#       expr=' IM0*IM1')

# immath(imagename=[image6+'.smallks.rebin.S2N','maskS2Neq4'],
#       mode='evalexpr',outfile=image6+'.smallks.rebin.S2N.masked',
#       expr=' IM0*IM1')

# immath(imagename=[image7+'.smallks.rebin.S2N','maskS2Neq4'],
#       mode='evalexpr',outfile=image7+'.smallks.rebin.S2N.masked',
#       expr=' IM0*IM1')

# immath(imagename=[image8+'.smallks.rebin.S2N','maskS2Neq4'],
#       mode='evalexpr',outfile=image8+'.smallks.rebin.S2N.masked',
#       expr=' IM0*IM1')

# immath(imagename=[image9+'.smallks.rebin.S2N','maskS2Neq4'],
#       mode='evalexpr',outfile=image9+'.smallks.rebin.S2N.masked',
#       expr=' IM0*IM1')

# #write to fits

# exportfits(imagename=image2+'.smallks.rebin.S2N.masked',
#    dropdeg=True,fitsimage=image2+'_ksrebin_S2N_masked.fits')

# exportfits(imagename=image4+'.smallks.rebin.S2N.masked',
#    dropdeg=True,fitsimage=image4+'_ksrebin_S2N_masked.fits')

# exportfits(imagename=image5+'.smallks.rebin.S2N.masked',
#    dropdeg=True,fitsimage=image5+'_ksrebin_S2N_masked.fits')

# exportfits(imagename=image6+'.smallks.rebin.S2N.masked',
#    dropdeg=True,fitsimage=image6+'_ksrebin_S2N_masked.fits')

# exportfits(imagename=image7+'.smallks.rebin.S2N.masked',
#    dropdeg=True,fitsimage=image7+'_ksrebin_S2N_masked.fits')

# exportfits(imagename=image8+'.smallks.rebin.S2N.masked',
#    dropdeg=True,fitsimage=image8+'_ksrebin_S2N_masked.fits')

# exportfits(imagename=image9+'.smallks.rebin.S2N.masked',
#    dropdeg=True,fitsimage=image9+'_ksrebin_S2N_masked.fits')

