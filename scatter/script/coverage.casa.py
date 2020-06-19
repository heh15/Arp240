import os
from subprocess import call

Dir='/1/home/heh15/workingspace/Arp240/scatter/'
picDir=Dir+'picture/'
imageDir=Dir+'image/'
logDir=Dir+'log/'
scriptDir=Dir+'script/'
workDir=imageDir+'temp/'


#### for NGC 5257 ####
## moment 2 map

# goes into the directory
os.chdir(imageDir)
call(['mkdir','temp'])
os.chdir(workDir)

rmtables('*')
# start to check the coverage
im1=imageDir+'NGC5257_12CO21_combine_noise40.pbcor.mom0/'
im2=imageDir+'NGC5257_12CO21_combine_noise40.image.mom2/'

# making a mask to remove pixels with limited coverage
immath(imagename=im2,
       mode='evalexpr',outfile='premask1',
       expr='IM0/IM0')
ia.open('premask1')
ia.replacemaskedpixels(0.0, update=True)
ia.done()
ia.close()

imrebin(imagename='premask1',
     outfile='premask1.rebin',
        factor=[5,5])

# require at least 50% coverage

ia.open('premask1.rebin')
ia.calcmask(mask='premask1.rebin >= 0.5', name='mymask1.rebin')
ia.done()
ia.close()

immath(imagename='premask1.rebin',
       mode='evalexpr',outfile='NGC5257_mask50percent',
       expr=' IM0/IM0')

call(['cp','-r','NGC5257_mask50percent','..'])
rmtables('*')

#### for NGC 5258 ####
rmtables('*')

# start to check the coverage
im1=imageDir+'NGC5258_12CO21_combine.pbcor.mom0/'
im2=imageDir+'NGC5258_12CO21_combine_noise45.mom2/'

# making a mask to remove pixels with limited coverage
immath(imagename=im2,
       mode='evalexpr',outfile='premask1',
       expr='IM0/IM0')
ia.open('premask1')
ia.replacemaskedpixels(0.0, update=True)
ia.done()
ia.close()

imrebin(imagename='premask1',
     outfile='premask1.rebin',
        factor=[5,5])

# require at least 50% coverage

ia.open('premask1.rebin')
ia.calcmask(mask='premask1.rebin >= 0.5', name='mymask1.rebin')
ia.done()
ia.close()

immath(imagename='premask1.rebin',
       mode='evalexpr',outfile='NGC5258_mask50percent',
       expr=' IM0/IM0')

call(['cp','-r','NGC5258_mask50percent','..'])

# move back to the image directory
os.chdir(imageDir)
call(['rm','-rf','temp'])

