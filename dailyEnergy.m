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

%adjust variables to save here---------------
savevars={'M','direct','diffuse','Lin','presZ','albedo',...
    'windspd','Ta','ea'};
%default number of times during the day to process (24)
num_times=length(gldasInterp.datevalsUTC);
switch mode
    case 'normal'
        FOREST=varargin{1};
        sFile=varargin{2};
        todays_dateval=floor(gldasInterp.datevalsLocal(1));
       %check for deltavis 
       group = findMODISh5group(sFile,'500m');
       i=h5info(sFile,group);
       D=i.Datasets;
       dvisflag=false;ii=1;
       while ii < length(D) && ~dvisflag
           if strcmp(D(ii).Name,'deltavis')
           dvisflag=true;
           end
           ii=ii+1;
       end

        %get grain size, in micrometers, for today
        [grain_size,~,hdr]=GetEndmember(sFile,'grain_size',todays_dateval);
        % fix min/max grain sizes
        grain_size (grain_size==0)=NaN;
        grain_size(grain_size < 10) = 10;
        grain_size(grain_size > 1100) = 1100;
        %either use deltavis (reflectance,MODDRFS) or dust (ppmw, SPIRES)
        if dvisflag
            deltavis=GetEndmember(sFile,'deltavis',todays_dateval);
        else
            dust=GetEndmember(sFile,'dust',todays_dateval);
        end

        %if the topo hdr doesn't match the fsca hdr, reproject
        if ~(isequal(hdr.RefMatrix,topo.hdr.RefMatrix) &&...
            isequal(hdr.RasterReference.RasterSize,...
            topo.hdr.RasterReference.RasterSize))
        
            grain_size=reprojectRaster(grain_size,hdr.RefMatrix,...
                hdr.ProjectionStructure,topo.hdr.ProjectionStructure,...
                'rasterref',topo.hdr.RasterReference);
            if dvisflag
            deltavis=reprojectRaster(deltavis,hdr.RefMatrix,...
                hdr.ProjectionStructure,topo.hdr.ProjectionStructure,...
                'rasterref',topo.hdr.RasterReference);
            else
            dust=reprojectRaster(dust,hdr.RefMatrix,...
               hdr.ProjectionStructure,topo.hdr.ProjectionStructure,...
               'rasterref',topo.hdr.RasterReference);
            end
        end
        if dvisflag
            sw_opt_input(1:4)={FOREST,grain_size,dvisflag,deltavis};
        else
            sw_opt_input(1:4)={FOREST,grain_size,dvisflag,dust};
        end
        
        ebalance_opt_input=FOREST;
        if ~metvars_flag 
            savevars={'M'};
        end
    case 'debris'
        %albedo
        sw_opt_input=varargin{1};
        d=varargin{2};
        ebalance_opt_input=d;
        savevars={'G'};
        %daily averages
        num_times=1;
        fn=fieldnames(gldasInterp);  
        for i=1:length(fn)
            fni=fn{i};
            switch fni
                case {'datevalsUTC' , 'datevalsLocal'}
                    %1d mean
                    gldasInterp.(fni)=mean(gldasInterp.(fni),'omitnan');
                otherwise %3d mean
                    gldasInterp.(fni)=mean(gldasInterp.(fni),3,'omitnan');
            end
        end
        fn=fieldnames(ceresInterp);
        for i=1:length(fn)
            fni=fn{i};
            switch fni
                case {'datevalsUTC' , 'datevalsLocal'}
                    %1d mean
                    ceresInterp.(fni)=mean(ceresInterp.(fni),'omitnan');
                otherwise %3d mean
                    ceresInterp.(fni)=mean(ceresInterp.(fni),3,'omitnan');
            end
        end
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

in.sw=ceresInterp.adj_atmos_sw_down_all_surface_1h;
in.lw=ceresInterp.adj_atmos_lw_down_all_surface_1h;
    
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
    Kd=1; %W/(m*K)
    out.G=Kd*(out.Tsfc-273.15)./(0.5*d); %positive values for heat going out of debris
    %and into ice below
case 'debris depth'
%     out.Tsfc=Tsfc;
    out.d=out.opt_output; %depth in m
    %switched to daily avg
%     out.G=(out.Tsfc-273.15)./out.d;
    out.d=out.d*1000; %convert to mm
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
        t=isnan(out.albedo); 
        x=uint8(out.albedo.*100); %0-100 pct
        x(t)=intmax('uint8');
    case 'presZ'
        t=isnan(out.presZ); 
        x=uint8(out.presZ); %kPa
        x(t)=intmax('uint8');
    case 'Ta'
        t=isnan(out.Ta); 
        x=out.(svar)-273.15 ; %K to deg C
        x=int8(x);
        x(t)=intmin('int8');
    case 'windspd'
        t=isnan(out.windspd); 
        x=uint8(out.windspd); %m/s
        x(t)=intmax('int8');
%     case 'rh' %different since rh has to be computed
%         t=isnan(out.ea); %null values from ea (Pa)
%         es=SaturationVaporPressure(out.Ta,'ice'); 
%         %input in K, output in kPa
%         rh=(out.ea./1000)./es*100;
%         x=uint8(rh); %0-100 pct
%         x(t)=intmax(x);      
    otherwise % 'M','direct','diffuse','Lin'
        x=out.(svar);
        t=isnan(x);
        x=uint16(x);
        x(t)=intmax('uint16');
end
m.(svar)=x;
end