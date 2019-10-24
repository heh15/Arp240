############################################################
# concatinating data

split(vis='arp240_7m_CalibratedData.ms',
      outputvis='arp240_7m_CO.ms',spw='0',field='1',
      datacolumn='data',keepflags=False)
   
split(vis='arp240_12m_CalibratedData.ms',
      outputvis='arp240_12m_CO.ms',spw='0',field='1',
      datacolumn='data',keepflags=False)

plotms(vis='arp240_7m_CalibratedData.ms',spw='0',xaxis='channel',yaxis='amp',avgtime='1e8',avgscan=True,showgui=True)
plotms(vis='arp240_12m_CalibratedData.ms',spw='0',xaxis='channel',yaxis='amp',avgtime='1e8',avgscan=True,showgui=True)
      
# uvcontsub(vis='arp240_7m_CO.ms',fitorder=1,fitspw='0:0~500,0:700~900')
# uvcontsub(vis='arp240_12m_CO.ms',fitorder=1,fitspw='0:0~500,0:700~900')   
concat(vis=['arp240_7m_CO.ms','arp240_12m_CO.ms'],
       concatvis='arp240_combine_CO.ms')

# #############################################
# Image Parameters


# source parameters
# ------------------

field='1' 
# imagermode='csclean' # uncomment if single field 
imagermode='mosaic'
# phasecenter='J2000 13h39m56s +0d50m00s' # uncomment and set to field number for phase

cell='0.3arcsec' # cell size for imaging.
imsize = [320,320] # size of image in pixels.

# velocity parameters
# -------------------

outframe='bary' # velocity reference frame. See science goals.
veltype='radio' # velocity type. 


# imaging control


weighting = 'briggs'
robust=0.5
niter=1000
threshold = '0.0mJy'

finalvis ='arp240_12m_CO.ms'
linevis = finalvis 

sourcename ='NGC_5257' # name of source
linename = 'CO10' # name of transition (see science goals in OT for name)
typename = '12m'
lineimagename = sourcename+'_'+linename+'_'+typename # name of line image

restfreq='112.73GHz' 

spw='0' 

start='-600km/s' # start velocity. See science goals for appropriate value.
width='10km/s' # velocity width. See science goals.
nchan =120


##############################################
# Image line emission [REPEAT AS NECESSARY]


clearcal(vis=linevis)
delmod(vis=linevis)

for ext in ['.flux','.image','.mask','.model','.pbcor','.psf','.residual','.flux.pbcoverage','.pb','.wtsum']:
    rmtables(lineimagename + ext)

clean(vis=linevis,
      field=field,
      imagename=lineimagename, 
      spw=spw,
#     phasecenter=phasecenter, # uncomment if mosaic.      
      mode='velocity',
      start=start,
      width=width,
      nchan=nchan, 
      outframe=outframe, 
      veltype=veltype, 
      restfreq=restfreq, 
      niter=niter,  
      interactive=False,
      cell=cell,
      imsize=imsize, 
      weighting=weighting, 
      robust=robust,
      imagermode=imagermode)


#1 iteration
# 2.31 e-02 Jy/beam 

# If you'd like to redo your clean, but don't want to make a new mask
# use the following commands to save your original mask. This is an
# optional step.
# linemaskname = 'line.mask'
## rmtables(linemaskname) # uncomment if you want to overwrite the mask.
# os.system('cp -ir ' + lineimagename + '.mask ' + linemaskname)

##############################################
# Apply a primary beam correction

import glob

myimages = glob.glob("*.image")

rmtables('*.pbcor')
for image in myimages:
    pbimage = image.rsplit('.',1)[0]+'.flux'
    outfile = image.rsplit('.',1)[0]+'.pbcor'
    impbcor(imagename=image, pbimage=pbimage, outfile = outfile)

##############################################
# Export the images

import glob

myimages = glob.glob("*.pbcor")
for image in myimages:
    exportfits(imagename=image, fitsimage=image+'.fits',overwrite=True)

myimages = glob.glob("*.flux")
for image in myimages:
    exportfits(imagename=image, fitsimage=image+'.fits',overwrite=True) 

##############################################
# Create Diagnostic PNGs


os.system("rm -rf *.png")
mycontimages = glob.glob("calibrated*.image")
for cimage in mycontimages:
    max=imstat(cimage)['max'][0]
    min=-0.1*max
    outimage = cimage+'.png'
    os.system('rm -rf '+outimage)
    imview(raster={'file':cimage},out=outimage)


# this will have to be run for each sourcename
sourcename='NGC_5257' # insert source here, if it isn't already set
mylineimages = glob.glob(sourcename+"*.image")
for limage in mylineimages:
    rms=imstat(limage,chans='1')['rms'][0]
    mom0=limage+'.mom0'
    os.system("rm -rf "+mom0)
    immoments(limage,moments=[0],outfile=mom0)
    os.system("rm -rf "+mom0+".png")
    imview(raster={'file':mom0},out=mom0+'.png')
