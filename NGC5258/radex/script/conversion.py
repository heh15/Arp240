D=98;D_cm=98*1e6*3.08e18
alpha=4.3;arcsec2pc=485

intensity=23.46
bmaj=3.8;bmin=2.99;nu=112.73
intensity_K=intensity/(0.0109*bmaj*bmin)/(nu/115.27)**2
mass=intensity_K*alpha

mass_R=10**17.7942*5000*3.32e-24*(3.1e18)**2/2e33*1.4
alpha_true=alpha*mass_R/mass
