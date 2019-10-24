#!/usr/bin/env python
# coding: utf-8

# # To calculate the scale height with CASA

# ## Directory

# In[12]:


Dir='/1/home/heh15/workingspace/Arp240/scatter/'
picDir=Dir+'picture/'
imageDir=Dir+'image/'
logDir=Dir+'log/'
scriptDir=Dir+'script/'


# ## NGC 5257

# ### primary beam correction for moment 0 map 

# In[14]:


imcollapse(imagename=imageDir+'NGC5257_12CO21_combine_noise40.pb/',
           function="mean",
           axes=3,
           outfile=imageDir+'NGC5257_12CO21_combine_mom0.pb')
impbcor(imagename=imageDir+'NGC5257_12CO21_combine_noise40.image.mom0/',
        pbimage=imageDir+'NGC5257_12CO21_combine_mom0.pb/',
        outfile=imageDir+'NGC5257_12CO21_combine_noise40.pbcor.mom0/')


# Rebin the moment 0 map and export the fits image. 

# In[ ]:


imrebin(imagename=imageDir+'NGC5257_12CO21_combine_noise40.pbcor.mom0/',
        factor=[5,5],
        outfile=imageDir+'NGC5257_12CO21_pbcor_rebin.mom0')
exportfits(imagename=imageDir+'NGC5257_12CO21_pbcor_rebin.mom0/',fitsimage='NGC5257_12CO21_pbcor_rebin_mom0.fits')


# ### calculate the scale height

# In[ ]:


im1='NGC5257_12CO21_combine_noise40.pbcor.mom0/'
im2='NGC5257_12CO21_combine_noise40.image.mom2/'
immath(imagename=[imageDir+im1,imageDir+im2],
      expr='IM1*IM1/IM0',
      outfile=imageDir+'NGC5257_height.image')


# rebin the calculated image

# In[ ]:


imrebin(imagename=imageDir+'NGC5257_height.image/',factor=[5,5],outfile=imageDir+'NGC5257_height_rebin.image')
exportfits(imagename=imageDir+'NGC5257_height_rebin.image',fitsimage=imageDir+'NGC5257_height_rebin.fits')

### rebin the moment 2 map. 
imrebin(imagename=imageDir+'NGC5257_12CO21_combine_noise40.image.mom2/',factor=[5,5],outfile=imageDir+'NGC5257_12CO21_combine_mom2_rebin.mom2')  
exportfits(imagename=imageDir+'NGC5257_12CO21_combine_mom2_rebin.mom2',fitsimage=imageDir+'NGC5257_12CO21_combine_mom2_rebin.fits')


### Check the mask



## Calculate the scale height of NGC 5258

imcollapse(imagename=imageDir+'NGC5258_12CO21_combine_noise45.pb',
           function="mean",
           axes=3,
           outfile=imageDir+'NGC5258_12CO21_combine_mom0.pb')
impbcor(imagename=imageDir+'NGC5258_12CO21_combine_noise45.image.mom0/',
        pbimage=imageDir+'NGC5258_12CO21_combine_mom0.pb/',
        outfile=imageDir+'NGC5258_12CO21_combine.pbcor.mom0/')


imrebin(imagename=imageDir+'NGC5258_12CO21_combine.pbcor.mom0/',
        factor=[5,5],
        outfile=imageDir+'NGC5258_12CO21_pbcor_rebin.mom0')
exportfits(imagename=imageDir+'NGC5258_12CO21_pbcor_rebin.mom0/',fitsimage=imageDir+'NGC5258_12CO21_pbcor_rebin_mom0.fits')

im1='NGC5258_12CO21_combine.pbcor.mom0/'
im2='NGC5258_12CO21_combine_noise45.mom2/'
immath(imagename=[imageDir+im1,imageDir+im2],
      expr='IM1*IM1/IM0',
      outfile=imageDir+'NGC5258_height.image')
imrebin(imagename=imageDir+'NGC5258_height.image/',factor=[5,5],outfile=imageDir+'NGC5258_height_rebin.image')
exportfits(imagename=imageDir+'NGC5258_height_rebin.image',fitsimage=imageDir+'NGC5258_height_rebin.fits')
