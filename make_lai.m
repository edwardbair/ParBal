function lai=make_lai(type, dateval)
%input: forest type raster where
% 1 coniferous
% 2 deciduous 
% dateval, matlab datetime

%output lai: leaf area index for two veg types

%type
vlai_summer=[2.5 2.5];
vlai_winter=[2.5 0.5];

% Note: A maximum forest LAI of 5.0 will give almost zero (like
% 10 W m^2) incoming solar under the canopy.  Values for Fraser
% Experimental Forest in Colorado are 2-3 (LSOS site = 1.8,
% Kelly's/Gus' site = 2.3).

% Calculate a seasonally varying temperature, assuming a max and
% min temperature and a cos distribution peaking in mid July
% (J_day = 200).  Then use this to define the seasonal lai
% variation.
%compute LAI

doy_1=datenum([year(dateval) 1 1]);
doy_end_yr=datenum([year(dateval) 12 31]);
daysinyr=doy_end_yr-doy_1+1;
doy=floor(dateval)-doy_1+1;
tmax = 298;
tmin = 273;
peak_doy = 200;

dtseason = tmax-tmin;
vtrans = tmin+dtseason./2;

tseason= vtrans + dtseason./2.*cos(2.*pi./daysinyr.*(doy-peak_doy));
fseason=0.0016.*(tmax - tseason).^2;
lai=zeros(size(type));

ntypes=length(vlai_summer);
for n=1:ntypes
    lai(type==n)=(1-fseason).*vlai_summer(n) + fseason.*vlai_winter(n);
end