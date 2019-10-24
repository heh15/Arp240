
# NOTE: you can make plotms go faster if you pick a single channel near the
#  center of the spw with your line

# checking uvrange of CO,CN

plotms(vis="i13120band3_115.ms",
   xaxis="uvwave",field="3",spw="0:450")

# continuous from 5 to 500 klambda, then some coverage up to 105 klambda

# and for HCN

plotms(vis="iras13120band3.ms",xaxis="uvwave",field="0",spw="0:200")

# continuous from 15 to 300 klambda

# Easiest to image the line with the larger beam first so you know what to aim for

For HCN cube I used

     uvrange='>15klambda',
     uvtaper=['0.0arcsec','0.8arcsec','71deg'],

# beam 1.03x0.49", PA=71 before uvtaper
# after uvtaper 1.03x0.83", PA=67 with rms=0.95 mJy/beam
# Getting the parameters set right for uvtaper is a bit of trial and error

For CN cube in tclean I used

     uvrange='>15klambda',
     uvtaper=['0.8arcsec','0.9arcsec','-23deg'],

# gave me a beam is 0.95x0.87, PA=-30 rms 1.4 mJy/beam
