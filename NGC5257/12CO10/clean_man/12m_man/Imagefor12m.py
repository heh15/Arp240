############################################################
# concatinating data

#split(vis='arp240_12m_CalibratedData.ms',
#      outputvis='arp240_12m_CO.ms',spw='0',field='0',
#      datacolumn='data',keepflags=False)


#plotms(vis='arp240_12m_CalibratedData.ms',spw='0',xaxis='channel',yaxis='amp',a vgtime='1e8',avgscan=True,showgui=True)

      
# contsub(vis='arp240_12m_CO.ms',fitorder=1,fitspw='0:0~500,0:700~900')


# #############################################
# Image Parameters


# source parameters
# ------------------

field='0' # science field(s). For a mosaic, select all mosaic fields. DO NOT LEAVE BLANK ('') OR YOU WILL TRIGGER A BUG IN CLEAN THAT WILL PUT THE WRONG COORDINATE SYSTEM ON YOUR FINAL IMAGE.
# imagermode='csclean' # uncomment if single field 
imagermode='mosaic' # uncomment if mosaic or if combining one 7m and one 12m pointing.
# phasecenter='J2000 13h39m56s +0d50m00s' # uncomment and set to field number for phase
                # center. Note lack of ''.  Use the weblog to
                # determine which pointing to use. Remember that the
                # field ids for each pointing will be re-numbered
                # after your initial split. You can also specify the
                # phase center using coordinates, e.g.,
                # phasecenter='J2000 19h30m00 -40d00m00'

# image parameters.
# ----------------




cell='0.5arcsec' # cell size for imaging.
imsize = [320,320] # size of image in pixels.

# velocity parameters
# -------------------

outframe='bary' # velocity reference frame. See science goals.
veltype='radio' # velocity type. 


# imaging control
# ----------------

# The cleaning below is done interactively, so niter and threshold can
# be controlled within clean. 

weighting = 'briggs'
robust=0.5
niter=1000
threshold = '0.0mJy'

##############################################
# Image line emission [REPEAT AS NECESSARY]


finalvis ='arp240_12m_CO.ms'
linevis = finalvis  # uncomment if you neither continuum subtracted nor self-calibrated your data.
# linevis = finalvis + '.contsub' # uncomment if continuum subtracted
# linevis = finalvis + '.contsub.selfcal' # uncommment if both continuum subtracted and self-calibrated
# linevis = finalvis + '.selfcal' # uncomment if just self-calibrated (no continuum subtraction)

sourcename ='NGC_5257' # name of source
linename = 'CO10' # name of transition (see science goals in OT for name)
typename = '12m'
lineimagename = sourcename+'_'+linename+'_'+typename # name of line image

restfreq='112.73GHz' # Typically the rest frequency of the line of
                        # interest. If the source has a significant
                        # redshift (z>0.2), use the observed sky
                        # frequency (nu_rest/(1+z)) instead of the
                        # rest frequency of the
                        # line.

spw='0' # uncomment and replace with appropriate spw if necessary.

start='-600km/s' # start velocity. See science goals for appropriate value.
width='10km/s' # velocity width. See science goals.
nchan =120
  # number of channels. See science goals for appropriate value.


# If necessary, run the following commands to get rid of older clean
# data.

clearcal(vis=linevis)
delmod(vis=linevis)

for ext in ['.flux','.image','.mask','.model','.pbcor','.psf','.residual','.flux.pbcoverage','.pb','.wtsum']:
    rmtables(lineimagename + ext)

clean(vis=linevis,
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
      threshold=threshold, 
      interactive=True,
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
    outimage = cimage+'.png'
    os.system('rm -rf '+outimage)
    imview(raster={'file':cimage},out=outimage)


# this will have to be run for each sourcename
sourcename='NGC_5257' # insert source here, if it isn't already set
mylineimages = glob.glob(sourcename+"*.image")
for limage in mylineimages:
    mom0=limage+'.mom0'
    os.system("rm -rf "+mom0)
    immoments(limage,moments=[0],outfile=mom0)
    os.system("rm "+mom0+".png")
    imview(raster={'file':mom0},out=mom0+'.png')
