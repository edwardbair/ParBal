function LinZ=lapseLongWave(Lin_coarse,Zdiff)
lapse=-.029; %from Hinkelman et al 2015, W m^-2 m^-1
LinZ=Lin_coarse+Zdiff*lapse;