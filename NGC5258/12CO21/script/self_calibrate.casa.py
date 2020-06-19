gaincal(vis="calibrated/NGC5258_CO12_21_7m.ms",
caltable="phase_spw0.cal",
field="0",
solint="30s",
calmode="p",
refant="CM01",
gaintype="G",
spw='1')

plotcal(caltable="phase_spw0.cal",
xaxis="time",
yaxis="phase",
subplot=331,
iteration="antenna",
plotrange=[0,0,-30,30],
markersize=5,
fontsize=10.0,
figfile="sis14_selfcal_phase_scan.png",
showgui = True)
