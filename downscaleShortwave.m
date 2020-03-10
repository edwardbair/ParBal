function [F,B,FinZ,BinZ,albedo,presZ,T_Z] = downscaleShortwave(...
    datevalUTC,pres,T,Gin,gldas_topo,topo,mode,opt_input)
% Calculates incoming shortwave radiation for debris cover
%INPUT
% datevalUTC - UTC mat datetime
% pres - coarse pressure from main source (ie GLDAS) in kPa, 
% reprojected to fine scale
% T - coarse temperature in K, "
% Gin - coarse incoming global shortwave, "
% gldas topo - topo structure for g/n ldas
% topo - fine scale topo structure
% mode - either 'normal' or 'debris', with normal for snow/icemelt and
% debris for icemelt under debris
% then follow with req'd inputs:
% albedo - albedo raster for debris
% opt_input
% for normal
%     FOREST
%     grain_size - grain size, um
%     deltavis - deltavis from MODSCAG, 0-1
% for debris or debris depth
%     albedo - raster of albedo values of debris, 0-1
%OUTPUT
% F, B - downscaled diffuse and beam components
% FinZ, BinZ - " " uncorrected for veg and topo, W/m^2
% albedo - albedo, 0-1
% presZ - downscaled air pressure using lapse rate, kPa
% T_Z - downscaled air temp, K

% sun angles over fine scale topographic grid w/ atmos. correction
Zdiff=gldas_topo.Zdiff;
lapse=-0.0065;% deg K/m
% Adjust air temperature from reference to fine using standard lapse
T_Z = T + (Zdiff.*lapse);
% Pressure at elevation adjusted with lapse and temp
presZ = elevationPressure(pres,lapse,Zdiff,T);
fineTSA  = TopoSunAngle(datevalUTC,topo,presZ,T_Z);

switch mode
    case 'normal'
        FOREST=opt_input{1};
        grain_size=opt_input{2};
        deltavis=opt_input{3};
        %albedo calculation
%         albedo=broadbandSnowAlbedo(grain_size,acosd(fineTSA.mu));
        albedo=spires_albedo(grain_size,fineTSA.mu);

%         albedo=albedo-deltavis/2;
        albedo=albedo-deltavis*0.63;
        %from SMARTS295Main using defaultZSMARTSinput('mlw', 
        %0.6346 (cosZ=0.6) to 0.6370 (cosZ=1)
    otherwise
        albedo=opt_input;
end

%fix negative values
Gin(Gin<0)=0;

% If sun below horizon, all incoming solar is diffuse
if ~any(reshape(fineTSA.mu0,[numel(fineTSA.mu0) 1]))
    FinZ=Gin;
    BinZ=zeros(size(FinZ));
    % if Global incoming all zeros
elseif ~any(reshape(Gin,[numel(Gin) 1]))
    FinZ=zeros(size(Gin));
    BinZ=FinZ;
else
    % Exoatmospheric from solar constant and radius vector
    solarconstant=1367;
    %scale to radius vector
    S0 = solarconstant./(fineTSA.RadiusVector.^2);
    %exoatmospheric flux,compute w/o atmos
    exoflux=S0.*fineTSA.mu0_unrefracted;
    % Total Transmittance at surface
    %fix zeros
    Gin(exoflux==0)=0;
    %assume a max T of 0.95 when Gin>exoflux
    Gin(Gin>exoflux)=0.95*exoflux(Gin>exoflux);
    T=Gin./exoflux;
   % Diffuse component
    Fin = partitionGlobalErbs1982(Gin,T);
    % Optical depth
    tau=-log(T)./fineTSA.airmass;
    % happens when T=0
    tau(~isfinite(tau))=NaN;
    tau(tau < 0)=NaN;
    % or Gin0s>exoflux, set NaN to nanmedian tau0f
    tau(isnan(tau))=nanmedian(tau(:));
    % Adjust optical depth to elevation using pressure   
    tauZ=(presZ./pres).*tau;
    % Beam at elevation
    BinZ = fineTSA.mu0.*S0.*exp(-tauZ.*fineTSA.airmass);
    % diffuse is obtained from an empirical formulation
    FinZ = diffuseAtElevation_modis(Fin,tau,pres,fineTSA.mu0,tauZ,presZ);
end

% Adjust for topography
BinT=fineTSA.mu./fineTSA.mu0.*BinZ;
FinT=FinZ.*topo.view;

% add reflected radiation into diffuse before veg. correction
% terrain config. factor, Dozier and Frew [1990]; 
%Dubayah and Loechel [1997] Eq. 3 
Ct=(1+cosd(topo.slope))./2-topo.view;
%surrounding terrain
RinT=Ct.*albedo.*(FinT.*(1-topo.view)+BinT);
RinT(RinT<0 | isnan(RinT))=0;
%include in diffuse fixed 5/9/16
FinT_R=RinT+FinT;
F=FinT_R;
B=BinT;

if strcmp(mode,'normal')
% Adjust for canopy if normal mode
F=FOREST.cc.*F.*FOREST.tau+(1-FOREST.cc).*F;
B=FOREST.cc.*B.*exp((-FOREST.u.*FOREST.h)./fineTSA.mu0) +...
    (1-FOREST.cc).*B;
end

F(F<0 | isnan(F))=0;
B(B<0 | isnan(B))=0;
