############################################################
# concatinating data
# NGC5257
# 12m
split(vis='/home/heh15/Data/Arp240/Arp240/arp240_all.ms/',
      outputvis='NGC5257_CO12_21_12m.ms',field='1~3',
      datacolumn='data',keepflags=False)
# 7m
split(vis='/home/heh15/Data/Arp240/Arp240/arp240co21_aca.ms/',
      outputvis='NGC5257_CO12_21_7m.ms',spw='0,4,8',field='1',
      datacolumn='data',keepflags=False)

#NGC5258
# 12m
split(vis='/home/heh15/Data/Arp240/Arp240/arp240_all.ms/',
      outputvis='NGC5258_CO12_21_12m.ms',spw='0',field='4~10',
      datacolumn='data',keepflags=False)
# 7m
split(vis='/home/heh15/Data/Arp240/Arp240/arp240co21_aca.ms/',
      outputvis='NGC5258_CO12_21_7m.ms',spw='0,4,8',field='3~5',
      datacolumn='data',keepflags=False)

plotms(vis='NGC5257_CO12_21_combine.ms',spw='0',xaxis='channel',yaxis='amp',avgtime='1e8',avgscan=True,showgui=True)

plotms(vis='arp240_12m_CalibratedData.ms',spw='0',xaxis='channel',yaxis='amp',avgtime='1e8',avgscan=True,showgui=True)
      
uvcontsub(vis='arp240_CO12_21_7m.ms',fitorder=1,fitspw='0:0~200,0:450~800,1:0~200,1:450~800,2:0~200,2:450~800')
uvcontsub(vis='arp240_CO12_21_12m.ms',fitorder=1,fitspw='0:0~200,0:450~800')

concat(vis=['NGC5257_7m.ms','NGC5257_12m.ms'],
       concatvis='NGC5257_combine.ms')

concat(vis=['NGC5258_7m.ms','NGC5258_12m.ms'],
       concatvis='NGC5258_combine.ms')

fitspw='0:450~800,1,2:400~800,3:200:600'\
       '4:450~800,5,6:400~800,7:200:600'\
       '8:450~800,9,10:400~800,11:200:600'\
       '12:450~800,13,14:400~800,15:200:600'

uvcontsub(vis='NGC5257_combine.ms',fitorder=1,fitspw=fitspw,spw='0,4,8,12',combine='spw')


uvcontsub(vis='NGC5257_CO12_21_combine.ms',fitorder=1,fitspw='0:0~200,0:450~800,1:0~200,1:450~800,2:0~200,2:450~800,3:0~200,3:450~800',spw='0,1,2,3',combine='spw')

# check the combined data.

plotms(vis='arp240_CO12_21_combine.ms',spw='0',xaxis='uvdist',yaxis='amp',avgtime='1e8',avgscan=True,showgui=True)
plotms(vis='NGC5257_CO12_21_combine.ms',spw='0',xaxis='uvdist',yaxis='wt',avgtime='1e8',avgscan=True,showgui=True)

