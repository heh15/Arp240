# check if there is CN emmission. 

Dir='/1/home/heh15/workingspace/Arp240/NGC5257/CN10/'
calDir=Dir+'calibrated/'

CODir='/1/home/heh15/workingspace/Arp240/NGC5257/12CO10/calibrated/'

plotms(vis=CODir+'NGC5257_combine.ms',xaxis='frequency',yaxis='amp',avgscan=True,avgtime='1e6',spw='1,3',avgspw=True)

# uvcontsub(vis='NGC5258_combine.ms',fitorder=1,fitspw='0:0~500,0:700~900,1,2:0~500;700~900,3:',spw='0,2')
