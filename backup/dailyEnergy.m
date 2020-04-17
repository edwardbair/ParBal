<<<<<<< HEAD
function dailyEnergy(topo,gldasInterp,gldas_topo,ceresInterp,...,
    ceres_topo,fast_flag,metvars_flag,mode,outfile,varargin)
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
% metvars_flag - true - output hourly met vars, false - don't output hourly
% metvars
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
% then met vars if metvars_flag=true
% 'directZ' - flat surface direct sw, w/m^2
% 'diffuseZ' - flat surface diffuse sw, w/m^2
% 'LinZ' - flat surface incoming longwave, w/m^2
% 'presZ' - air pressure, mb
% 'albedo' - albedo, percent
% 'windspd' - wind speed, m/s *10
% 'Ta' - air temp, K * 10
% 'Td' - dew point temp, K * 10
% 'ea' - vapor pressure of air, mb

switch mode
    case 'normal'
        FOREST=varargin{1};
        sFile=varargin{2};
        todays_dateval=floor(gldasInterp.datevalsLocal(1));
        %get grain size, in micrometers, for today
        [grain_size,~,hdr]=GetEndmember(sFile,'grain_size',todays_dateval);
        deltavis=GetEndmember(sFile,'deltavis',todays_dateval);
        % fix min/max grain sizes
        grain_size (grain_size==0)=NaN;
        grain_size(grain_size < 10) = 10;
        grain_size(grain_size > 1100) = 1100;
        %convert to mm
        grain_size=grain_size.*1e-3;
        %if the topo hdr doesn't match the fsca hdr, reproject
        if ~isequal(topo.hdr,hdr)
            grain_size=reprojectRaster(grain_size,hdr.RefMatrix,...
                hdr.ProjectionStructure,topo.hdr.ProjectionStructure,...
                'rasterref',topo.hdr.RasterReference);
            deltavis=reprojectRaster(deltavis,hdr.RefMatrix,...
                hdr.ProjectionStructure,topo.hdr.ProjectionStructure,...
                'rasterref',topo.hdr.RasterReference);
        end      
        sw_opt_input(1:3)={FOREST,grain_size,deltavis};
        ebalance_opt_input=FOREST;
        if metvars_flag
            savevars={'M','directZ','diffuseZ','LinZ','presZ','albedo',...
                'windspd','Ta','Td','ea'};
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

%get LDASOnlyFlag
LDASOnlyFlag=false;
if isempty(ceresInterp)
    LDASOnlyFlag=true;
end

num_times=length(gldasInterp.datevalsUTC);
sizes=[topo.hdr.RasterReference.RasterSize num_times];
%determine whether NLDAS or GLDAS based on windfields
gldasflag=false;
if isfield(gldasInterp,'Wind_f_inst')
    gldasflag=true;
end
out.M = zeros(sizes,'single');
out.Lin = zeros(sizes,'single');
out.LinZ = zeros(sizes,'single');
out.Lout = zeros(sizes,'single');
out.Ta = zeros(sizes,'single');
out.diffuse = zeros(sizes,'single');
out.direct = zeros(sizes,'single');
out.diffuseZ = zeros(sizes,'single');
out.directZ = zeros(sizes,'single');
out.presZ = zeros(sizes,'single');
out.albedo= zeros(sizes,'single');
out.sensible = zeros(sizes,'single');
out.latent = zeros(sizes,'single');
out.windspd = zeros(sizes,'single');
out.G = zeros(sizes,'single');
out.Td = zeros(sizes,'single');
out.ea = zeros(sizes,'single');
out.Tsfc = zeros(sizes,'single');
out.opt_output = zeros(sizes,'single');
windS = struct();

% clean up input names
if gldasflag %gldas name
    in.pres=gldasInterp.Psurf_f_inst./1000; %Pa to KPa
    in.Ta=gldasInterp.Tair_f_inst; %Air temp, K
    in.Q=gldasInterp.Qair_f_inst; %specific humidity, Kg/Kg
else %nldas names 
   in.pres=gldasInterp.PRES./1000; %Pa to KPa
   in.Ta=gldasInterp.TMP; %Air temp, K
   in.Q=gldasInterp.SPFH; %specific humidity, Kg/Kg 
end
    
if LDASOnlyFlag
    in.sw=gldasInterp.SWdown_f_tavg;% Downward shortwave radiation flux W m-2
    in.lw=gldasInterp.LWdown_f_tavg;% Downward longwave radiation flux W m-2
    ceres_topo.Zdiff=gldas_topo.Zdiff;
else
    %ceres eliminated the surface pressure variable in ed 4a
in.sw=ceresInterp.adj_atmos_sw_down_all_surface_1h;
in.lw=ceresInterp.adj_atmos_lw_down_all_surface_1h;
% ceres.var={'sfc_comp_sw-down_all_3h','sfc_comp_lw-down_all_3h',...
% 'aux_surfpress_3h'};
%     in.sw=ceresInterp.sfc_comp_sw_down_all_3h;
%     in.lw=ceresInterp.sfc_comp_lw_down_all_3h;
    %ceres stopped supplying pressure w/ ed 4a
%     in.aux_pres=ceresInterp.aux_surfpress_3h./10; %hPa to kPa
    
end

%hourly
for h=1:num_times
    UTCstring=datestr(gldasInterp.datevalsUTC(h),'yyyy-mm-dd HH:MM');
    Localstring=datestr(gldasInterp.datevalsLocal(h),'yyyy-mm-dd HH:MM');
    [out.diffuse(:,:,h),out.direct(:,:,h),out.diffuseZ(:,:,h),out.directZ(:,:,h),...
        out.albedo(:,:,h),out.presZ(:,:,h),out.Ta(:,:,h)] = ...
        downscaleShortwave(gldasInterp.datevalsUTC(h),...
        in.pres(:,:,h),...
        in.Ta(:,:,h),...
        in.sw(:,:,h),...
        gldas_topo,topo,mode,sw_opt_input);
        fprintf('shortwave & albedo done for %s (%s UTC)\n',Localstring,UTCstring);
        %create wind structure depending on whether its GLDAS (speed only)
        %or NLDAS (U&V)
        %needed for parfor loop
    if ~gldasflag
        windS.U=gldasInterp.UGRD(:,:,h);
        windS.V=gldasInterp.VGRD(:,:,h);
        windS.windflag=false;
    else
        windS.windspd=gldasInterp.Wind_f_inst(:,:,h);
        windS.windflag=true;
    end
    %solve the energy balance
        [out.M(:,:,h),out.Tsfc(:,:,h),out.Lin(:,:,h),out.LinZ(:,:,h),...
            out.Lout(:,:,h),out.sensible(:,:,h),out.latent(:,:,h),...
            out.G(:,:,h),out.windspd(:,:,h),out.Td(:,:,h),out.ea(:,:,h),...
            out.opt_output(:,:,h)]=...
            run_ebalance(floor(gldasInterp.datevalsLocal(1)),...
            in.lw(:,:,h),...
            in.Q(:,:,h),...
            ceres_topo.Zdiff,...
            in.Ta(:,:,h),...
            out.Ta(:,:,h),...
            in.pres(:,:,h),...
            out.presZ(:,:,h),...
            out.albedo(:,:,h),...
            out.direct(:,:,h)+out.diffuse(:,:,h),...
            topo,windS,fast_flag,mode,ebalance_opt_input);
        fprintf('E-balance solved for %s (%s UTC)\n',Localstring,UTCstring);    
end

%use daily average surface temperatures to compute G for debris melt
% MATLABdates=gldasInterp.datevalsLocal;

switch mode
case 'normal'
    %sum for the day to reduce storage
    %Jepsen et al 2012 cold content treatment, Eq. 4 p3
    cold_content=cumsum(out.M,3);
    x=min(cold_content(:),out.M(:));
    x=reshape(x,size(out.M));
    x(x<0)=0;
    out.M=sum(x,3);
case 'debris'
%     out.Tsfc=Tsfc; 
    %Assume Kd = 1 W/(m*K) and 1/2 of debris depth to 0 deg debris temp
    out.G=(out.Tsfc-273.15)./(0.5.*d);
case 'debris depth'
%     out.Tsfc=Tsfc;
    out.d=out.opt_output; %depth in m
    out.G=(out.Tsfc-273.15)./(0.5.*out.d);
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
svar=savevars{i};
switch svar
    case 'albedo' 
        m.(svar)=int8(out.(svar).*100); %0-100 pct
    case 'presZ'
        m.(svar)=int16(out.(svar).*10); %(kPa to mb)
    case 'Ta'
        m.(svar)=int16(out.(svar).*10); %deg K * 10
    case 'windspd'
        m.(svar)=int16(out.(svar).*10); %m/sec * 10
    case 'ea'
        m.(svar)=int16(out.(svar).*10); %(Pa*10)
    case 'Td'
        m.(svar)=int16(out.(svar).*10); %deg K * 10
    otherwise
        m.(svar)=int16(out.(svar));
end
=======
function dailyEnergy(topo,gldasInterp,gldas_topo,ceresInterp,...,
    ceres_topo,fast_flag,metvars_flag,mode,outfile,varargin)
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
% metvars_flag - true - output hourly met vars, false - don't output hourly
% metvars
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
% then met vars if metvars_flag=true
% 'directZ' - flat surface direct sw, w/m^2
% 'diffuseZ' - flat surface diffuse sw, w/m^2
% 'LinZ' - flat surface incoming longwave, w/m^2
% 'presZ' - air pressure, mb
% 'albedo' - albedo, percent
% 'windspd' - wind speed, m/s *10
% 'Ta' - air temp, K * 10
% 'Td' - dew point temp, K * 10
% 'ea' - vapor pressure of air, mb

switch mode
    case 'normal'
        FOREST=varargin{1};
        sFile=varargin{2};
        todays_dateval=floor(gldasInterp.datevalsLocal(1));
        %get grain size, in micrometers, for today
        [grain_size,~,hdr]=GetEndmember(sFile,'grain_size',todays_dateval);
        deltavis=GetEndmember(sFile,'deltavis',todays_dateval);
        % fix min/max grain sizes
        grain_size (grain_size==0)=NaN;
        grain_size(grain_size < 10) = 10;
        grain_size(grain_size > 1100) = 1100;
        %convert to mm
        grain_size=grain_size.*1e-3;
        %if the topo hdr doesn't match the fsca hdr, reproject
        if ~isequal(topo.hdr,hdr)
            grain_size=reprojectRaster(grain_size,hdr.RefMatrix,...
                hdr.ProjectionStructure,topo.hdr.ProjectionStructure,...
                'rasterref',topo.hdr.RasterReference);
            deltavis=reprojectRaster(deltavis,hdr.RefMatrix,...
                hdr.ProjectionStructure,topo.hdr.ProjectionStructure,...
                'rasterref',topo.hdr.RasterReference);
        end      
        sw_opt_input(1:3)={FOREST,grain_size,deltavis};
        ebalance_opt_input=FOREST;
        if metvars_flag
            savevars={'M','directZ','diffuseZ','LinZ','presZ','albedo',...
                'windspd','Ta','Td','ea'};
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

%get LDASOnlyFlag
LDASOnlyFlag=false;
if isempty(ceresInterp)
    LDASOnlyFlag=true;
end

num_times=length(gldasInterp.datevalsUTC);
sizes=[topo.hdr.RasterReference.RasterSize num_times];
%determine whether NLDAS or GLDAS based on windfields
gldasflag=false;
if isfield(gldasInterp,'Wind_f_inst')
    gldasflag=true;
end
out.M = zeros(sizes,'single');
out.Lin = zeros(sizes,'single');
out.LinZ = zeros(sizes,'single');
out.Lout = zeros(sizes,'single');
out.Ta = zeros(sizes,'single');
out.diffuse = zeros(sizes,'single');
out.direct = zeros(sizes,'single');
out.diffuseZ = zeros(sizes,'single');
out.directZ = zeros(sizes,'single');
out.presZ = zeros(sizes,'single');
out.albedo= zeros(sizes,'single');
out.sensible = zeros(sizes,'single');
out.latent = zeros(sizes,'single');
out.windspd = zeros(sizes,'single');
out.G = zeros(sizes,'single');
out.Td = zeros(sizes,'single');
out.ea = zeros(sizes,'single');
out.Tsfc = zeros(sizes,'single');
out.opt_output = zeros(sizes,'single');
windS = struct();

% clean up input names
if gldasflag %gldas name
    in.pres=gldasInterp.Psurf_f_inst./1000; %Pa to KPa
    in.Ta=gldasInterp.Tair_f_inst; %Air temp, K
    in.Q=gldasInterp.Qair_f_inst; %specific humidity, Kg/Kg
else %nldas names 
   in.pres=gldasInterp.PRES./1000; %Pa to KPa
   in.Ta=gldasInterp.TMP; %Air temp, K
   in.Q=gldasInterp.SPFH; %specific humidity, Kg/Kg 
end
    
if LDASOnlyFlag
    in.sw=gldasInterp.SWdown_f_tavg;% Downward shortwave radiation flux W m-2
    in.lw=gldasInterp.LWdown_f_tavg;% Downward longwave radiation flux W m-2
    ceres_topo.Zdiff=gldas_topo.Zdiff;
else
    %ceres eliminated the surface pressure variable in ed 4a
in.sw=ceresInterp.adj_atmos_sw_down_all_surface_1h;
in.lw=ceresInterp.adj_atmos_lw_down_all_surface_1h;
% ceres.var={'sfc_comp_sw-down_all_3h','sfc_comp_lw-down_all_3h',...
% 'aux_surfpress_3h'};
%     in.sw=ceresInterp.sfc_comp_sw_down_all_3h;
%     in.lw=ceresInterp.sfc_comp_lw_down_all_3h;
    %ceres stopped supplying pressure w/ ed 4a
%     in.aux_pres=ceresInterp.aux_surfpress_3h./10; %hPa to kPa
    
end

%hourly
for h=1:num_times
    UTCstring=datestr(gldasInterp.datevalsUTC(h),'yyyy-mm-dd HH:MM');
    Localstring=datestr(gldasInterp.datevalsLocal(h),'yyyy-mm-dd HH:MM');
    [out.diffuse(:,:,h),out.direct(:,:,h),out.diffuseZ(:,:,h),out.directZ(:,:,h),...
        out.albedo(:,:,h),out.presZ(:,:,h),out.Ta(:,:,h)] = ...
        downscaleShortwave(gldasInterp.datevalsUTC(h),...
        in.pres(:,:,h),...
        in.Ta(:,:,h),...
        in.sw(:,:,h),...
        gldas_topo,topo,mode,sw_opt_input);
        fprintf('shortwave & albedo done for %s (%s UTC)\n',Localstring,UTCstring);
        %create wind structure depending on whether its GLDAS (speed only)
        %or NLDAS (U&V)
        %needed for parfor loop
    if ~gldasflag
        windS.U=gldasInterp.UGRD(:,:,h);
        windS.V=gldasInterp.VGRD(:,:,h);
        windS.windflag=false;
    else
        windS.windspd=gldasInterp.Wind_f_inst(:,:,h);
        windS.windflag=true;
    end
    %solve the energy balance
        [out.M(:,:,h),out.Tsfc(:,:,h),out.Lin(:,:,h),out.LinZ(:,:,h),...
            out.Lout(:,:,h),out.sensible(:,:,h),out.latent(:,:,h),...
            out.G(:,:,h),out.windspd(:,:,h),out.Td(:,:,h),out.ea(:,:,h),...
            out.opt_output(:,:,h)]=...
            run_ebalance(floor(gldasInterp.datevalsLocal(1)),...
            in.lw(:,:,h),...
            in.Q(:,:,h),...
            ceres_topo.Zdiff,...
            in.Ta(:,:,h),...
            out.Ta(:,:,h),...
            in.pres(:,:,h),...
            out.presZ(:,:,h),...
            out.albedo(:,:,h),...
            out.direct(:,:,h)+out.diffuse(:,:,h),...
            topo,windS,fast_flag,mode,ebalance_opt_input);
        fprintf('E-balance solved for %s (%s UTC)\n',Localstring,UTCstring);    
end

%use daily average surface temperatures to compute G for debris melt
% MATLABdates=gldasInterp.datevalsLocal;

switch mode
case 'normal'
    %sum for the day to reduce storage
    %Jepsen et al 2012 cold content treatment, Eq. 4 p3
    cold_content=cumsum(out.M,3);
    x=min(cold_content(:),out.M(:));
    x=reshape(x,size(out.M));
    x(x<0)=0;
    out.M=sum(x,3);
case 'debris'
%     out.Tsfc=Tsfc; 
    %Assume Kd = 1 W/(m*K) and 1/2 of debris depth to 0 deg debris temp
    out.G=(out.Tsfc-273.15)./(0.5.*d);
case 'debris depth'
%     out.Tsfc=Tsfc;
    out.d=out.opt_output; %depth in m
    out.G=(out.Tsfc-273.15)./(0.5.*out.d);
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
svar=savevars{i};
switch svar
    case 'albedo' 
        m.(svar)=int8(out.(svar).*100); %0-100 pct
    case 'presZ'
        m.(svar)=int16(out.(svar).*10); %(kPa to mb)
    case 'Ta'
        m.(svar)=int16(out.(svar).*10); %deg K * 10
    case 'windspd'
        m.(svar)=int16(out.(svar).*10); %m/sec * 10
    case 'ea'
        m.(svar)=int16(out.(svar).*10); %(Pa*10)
    case 'Td'
        m.(svar)=int16(out.(svar).*10); %deg K * 10
    otherwise
        m.(svar)=int16(out.(svar));
end
>>>>>>> a551fb190757d0789821ebee09340b8365d27973
end