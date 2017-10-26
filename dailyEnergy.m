function dailyEnergy(topo,gldasInterp,gldas_topo,ceresInterp,...,
    ceres_topo,fast_flag,mode,outfile,varargin)
% Calls daily radiation downscaling functions for each hour, then returns
% energy fluxes

% INPUTS
% topo - structure of topo data (DEM, Horizons, Slope, Aspect), all n * m *
% z size, with z = 1 for everything except Horizons
% gldasInterp - Hourly interpolated GLDAS data
% gldas_topo - GLDAS topographic structure
% ceresInterp - Hourly interpolated CERES data
% ceres_topo - ceres topo struct
% fast flag - true - only solve for M; false - solve for all outputs; only set for 'normal'
% mode - either 'normal', 'debris', or 'debris_depth with normal for snow/icemelt and
% debris for icemelt under debris (uses daily averages), and 'debris_depth' to solve for debris
% depth
% outfile - output filename
% then follow with req'd inputs:
% for 'normal';
%   FOREST - forest structure
%   sFile - h5 sca file
% for 'debris'
%   'albedo' - albedo of debris cover, n x m raster
%   'd' - debris depth raster (depth depth in meters)
% for 'debris_depth'
%   'albedo' - same as above
%   'Tsfc' - debris surface temperature, K
% OUTPUT
% for 'normal'
%   M - daily summed snow/ice melt, W/m^2
% for 'debris'
%   G - daily averaged debris cover melt, W/m^2 (should be averaged
% for 'debris_depth'
%   G,d - debris cover depth, mm 

switch mode
    case 'normal'
        FOREST=varargin{1};
        sFile=varargin{2};
        todays_dateval=floor(gldasInterp.datevalsLocal(1));
        %get grain size, in micrometers, for today
        grain_size=GetEndmember(sFile,'grain_size',todays_dateval);
        deltavis=GetEndmember(sFile,'deltavis',todays_dateval);
        % fix min/max grain sizes
        grain_size (grain_size==0)=NaN;
        grain_size(grain_size < 10) = 10;
        grain_size(grain_size > 1100) = 1100;
        %convert to mm
        grain_size=grain_size.*1e-3;
        sw_opt_input(1:3)={FOREST,grain_size,deltavis};
        ebalance_opt_input=FOREST;
        %save diagnostic variables for Sierra but not for other regions
        if contains(topo.topofile,'Sierra','IgnoreCase',true)
            savevars={'M','directZ','diffuseZ','LinZ','presZ','albedo','sensible','latent','windspd','Ta'};
        else
            savevars={'M'};
        end
    case 'debris'
        %albedo
        sw_opt_input=varargin{1};
        d=varargin{2};
        ebalance_opt_input=d;
        savevars={'G'};
    case 'debris depth'
        %albedo
        sw_opt_input=varargin{1};
        %Tsfc
        ebalance_opt_input=varargin{2};
        savevars={'G','d'};
end

num_times=length(gldasInterp.datevalsUTC);
sizes=[topo.hdr.RasterReference.RasterSize num_times];
%determine whether NLDAS or GLDAS based on windfields
windflag=false;
if isfield(gldasInterp,'WIND')
    windflag=true;
end
M = zeros(sizes,'single');
Lin = zeros(sizes,'single');
LinZ = zeros(sizes,'single');
Lout = zeros(sizes,'single');
Ta = zeros(sizes,'single');
diffuse = zeros(sizes,'single');
direct = zeros(sizes,'single');
diffuseZ = zeros(sizes,'single');
directZ = zeros(sizes,'single');
presZ = zeros(sizes,'single');
albedo= zeros(sizes,'single');
sensible = zeros(sizes,'single');
latent = zeros(sizes,'single');
windspd = zeros(sizes,'single');
G = zeros(sizes,'single');
Tsfc = zeros(sizes,'single');
opt_output = zeros(sizes,'single');

%hourly
for h=1:num_times
    UTCstring=datestr(gldasInterp.datevalsUTC(h),'yyyy-mm-dd HH:MM');
    Localstring=datestr(gldasInterp.datevalsLocal(h),'yyyy-mm-dd HH:MM');
    [diffuse(:,:,h),direct(:,:,h),diffuseZ(:,:,h),directZ(:,:,h),...
        albedo(:,:,h),presZ(:,:,h),Ta(:,:,h)] = ...
        downscaleShortwave(gldasInterp.datevalsUTC(h),...
        gldasInterp.PRES(:,:,h)./1000,...
        ceresInterp.aux_surfpress_3h(:,:,h)./10,...
        gldasInterp.TMP(:,:,h),...
        ceresInterp.sfc_comp_sw_down_all_3h(:,:,h),...
        gldas_topo,topo,mode,sw_opt_input);
        fprintf('shortwave & albedo done for %s (%s UTC)\n',Localstring,UTCstring);
        %create wind structure depending on whether its GLDAS (speed only)
        %or NLDAS (U&V)
        %needed for parfor loop
        windS = struct();
    if ~windflag
        windS.U=gldasInterp.UGRD(:,:,h);
        windS.V=gldasInterp.VGRD(:,:,h);
        windS.windflag=false;
    else
        windS.windspd=gldasInterp.WIND(:,:,h);
        windS.windflag=true;
    end
    %solve the energy balance
        [M(:,:,h),Tsfc(:,:,h),Lin(:,:,h),LinZ(:,:,h),...
            Lout(:,:,h),sensible(:,:,h),latent(:,:,h),...
            G(:,:,h),windspd(:,:,h),opt_output(:,:,h)]=...
            run_ebalance(floor(gldasInterp.datevalsLocal(1)),...
            ceresInterp.sfc_comp_lw_down_all_3h(:,:,h),...
            gldasInterp.SPFH(:,:,h),ceres_topo.Zdiff,gldasInterp.TMP(:,:,h),...
            Ta(:,:,h),gldasInterp.PRES(:,:,h)./1000,presZ(:,:,h),...
            albedo(:,:,h),direct(:,:,h)+diffuse(:,:,h),topo,windS,fast_flag,mode,...
            ebalance_opt_input);
        fprintf('E-balance solved for %s (%s UTC)\n',Localstring,UTCstring);    
end

%use daily average surface temperatures to compute G for debris melt
% MATLABdates=gldasInterp.datevalsLocal;

switch mode
case 'normal'
    %sum for the day to reduce storage
    %Jepsen et al 2012 cold content treatment, Eq. 4 p3
    cold_content=cumsum(M,3);
    x=min(cold_content(:),M(:));
    x=reshape(x,size(M));
    x(x<0)=0;
    out.M=sum(x,3);
case 'debris'
    out.Tsfc=Tsfc; 
    %Assume Kd = 1 W/(m*K) and 1/2 of debris depth to 0 deg debris temp
    out.G=(Tsfc-273.15)./(0.5.*d);
case 'debris depth'
    out.Tsfc=Tsfc;
    out.d=opt_output; %depth in m
    out.G=(Tsfc-273.15)./(0.5.*out.d);
    out.d=out.d*1000; %convert to mm
    %for both of these, G should be summed and multipled by 0.0108 mm / W
    %m^-2 hr or averaged and multipled by 0.2592 mm/ W m^-2 to damp out
    %diurnal Tsfc swings
end
    
%save
if exist(outfile,'file')
   delete(outfile);        
end
m = matfile(outfile,'Writable',true);
for i=1:length(savevars)
%     if ~strcmp(savevars{i},'M')
%         out.(savevars{i})=eval(savevars{i});
%     end
    m.(savevars{i})=int16(out.(savevars{i}));
end