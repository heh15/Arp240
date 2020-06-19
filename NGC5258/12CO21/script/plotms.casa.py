import numpy as np
import os
import scipy.ndimage as sni
import sys

vis='/home/heh15/practice/Arp240/' \
      'workingspace/CO12_21/calibrated/calibrated.ms/'
plotpath = '/home/heh15/practice/Arp240/'\
             'workingspace/CO12_21/picture/'
calibrator=['bandpass','Mars','phase','callisto']
fields=['0','1','2','10']

yaxes=['amp','phase']
xaxes=['time','frequency']
avgchannels=['0','1024']
avgtimes=['0','1e8']

observation='0'

xaxis='time'
yaxis='amp'
l=1

for j in range(len(xaxes)):
    xaxis=xaxes[j]
    avgtime=avgtimes[j]
    avgchannel=avgchannels[1-j]
    for k in range(len(yaxes)):
        yaxis=yaxes[k]
        plotfile=plotpath+calibrator[l]+'_'\
                  +yaxis+'_'+xaxis+'.png'
        field=fields[l]
        plotms(vis=vis,
               plotfile=plotfile,
               exprange='all',
               xaxis=xaxis,
               yaxis=yaxis,
               iteraxis='antenna',
               coloraxis='corr',
               observation=observation,
               field=field,
               showgui=False,
               avgchannel=avgchannel,
               avgtime=avgtime)
        
