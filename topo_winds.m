function [u_fine,v_fine]=topo_winds(topo,FOREST,u_coarse,v_coarse,dateval)

% input: topo structure w/ dem and curvature
%forest structure w/ veg info
% u_coarse,v_coarse, u and v components of wind speed (m sec^1)
%output u_fine,v_fine, terrain corrected and downscaled wind speeds to size
%of dem

% from Glen Liston's Snowmodel - Micromet Fortran code
% Liston et al. (2007) Simulating complex snow distrubutions
% in windy environments.
% Journal of Glaciology

slopewt=0.58;
curvewt=0.42;
beta=0.9;

%convert aspect to deg from north
azimuth=180-topo.aspect;
terrain_slope=topo.slope;
curvature=topo.curvature;

u_interp=imresize(u_coarse,size(topo.dem));
v_interp=imresize(v_coarse,size(topo.dem));

% convert back to speed W and direction theta in degrees
windspd=sqrt(u_interp.^2+v_interp.^2);
%fix 0 wind speed problems
windspd(windspd==0)=0.01;

winddir=270-atan2d(v_interp,u_interp);
winddir(winddir >= 360)= winddir(winddir >= 360)-360;

%calc wind slope direction
wind_slope=deg2rad(terrain_slope).*cos(deg2rad(winddir-azimuth));

%scale wind_slope -0.5:0.5 w/ min 1 mm to prevent divide by 0 errors

wslope_max=0.001;
wslope_max=max(wslope_max,max(abs(wind_slope(:))));

wind_slope=wind_slope./(2.*wslope_max);

% Calculate the wind speed and direction adjustments.  The
% curvature and wind_slope values range between -0.5 and +0.5.
% Valid slopewt and curvewt values are between 0 and 1, with
% values of 0.5 giving approximately equal weight to slope and
% curvature.  I suggest that slopewt and curvewt be set such
% that slopewt + curvewt = 1.0.  This will limit the total
% wind weight to between 0.5 and 1.5 (but this is not required).

windwt=1+slopewt.*wind_slope+curvewt.*curvature;

%fine wind speed grid
windspd=windwt.*windspd;

%account for canopy
% Define the canopy wind-weighting factor.  Assume z=0.6*canopy_ht,
% and the canopy_ht equals the vegetation snow-holding depth.
lai=make_lai(FOREST.type.num_val,dateval);

a=beta.*lai;
canopy_windwt=exp((-a).*(1.0-(0.6.*FOREST.shd)./FOREST.shd));
windspd=canopy_windwt.*windspd;

% fix wind dir according to Ryan (1977)
dirdiff=zeros(size(azimuth));

ind1=azimuth > 270 & winddir < 90;
if any(ind1(:))
dirdiff(ind1) = azimuth(ind1) - winddir(ind1) - 360;
end
ind2=azimuth < 90 & winddir > 270;
if any(ind2(:))
dirdiff(ind2)= azimuth(ind2) - winddir(ind2) + 360;
end
if any(~ind1(:) & ~ind2(:))
dirdiff(~ind1 & ~ind2) = azimuth(~ind1 & ~ind2) - winddir(~ind1 & ~ind2);
end

ind=dirdiff < 90;
if any(ind(:))
winddir(ind)=winddir(ind) - 0.5.*min(rad2deg(wind_slope(ind)),45).*...
    sin(deg2rad(2.*dirdiff(ind)));
end

ind2=winddir > 360 & ind;
if any(ind2(:));
winddir(ind2)=winddir(ind2)-360;
end

ind3=winddir < 0 & ind;
if any(ind3(:));
winddir(ind3)=winddir(ind3)+360;
end
%extract u and v
u_fine=-windspd.*sind(winddir);
v_fine=-windspd.*cosd(winddir);



