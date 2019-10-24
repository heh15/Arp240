thresh=1.5*rms
prename=preName+'_man'

tclean(vis=vis,
       imagename=prename,
       field=field,
       antenna=antenna,
       phasecenter=phasecenter,
       specmode=specmode,
       outframe=outframe,
       restfreq=restfreq,
       width=width,
       nchan=nchan,
       start=start,
       deconvolver=deconvolver,
       cell=cell,
       imsize=imsize,
       weighting=weighting,
       robust=robust,
       gridder=gridder,
       mask=myMask,
       niter=10000,
       threshold=str(thresh)+'Jy/beam',
       restoringbeam=restoringbeam,
       cyclefactor=cyclefactor,
       pblimit=pblimit,
       interactive=True)

myimage=prename+'.image'

immoments(imagename=myimage,moments=0,includepix=[2*rms,100],chans='5~30',outfile=myimage+'.mom0')
