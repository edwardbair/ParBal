function LDAS = read_LDAS(v,filename)
%LDAS = read_LDAS(v,filename)
% Read variable from NLDAS or GLDAS,requires read_grib package for NLDAS
% List variables in NLDAS or GLDAS: grib_struct=read_grib(filename,'inv')
%INPUT
% v - character array describing variable (see list of variables below)
% filename - .grb filename for NLDAS or GLDAS, character array
%OUTPUT
% LDASconus.matrix - NLDAS 2d, y by x
% LDASconus.Rmap - referencing matrix for NLDAS
% LDASconus.units - units of variable
% LDASconus.datevals - MATLAB date values in local time

%NLDAS
% Accumulation - 1 hr
% Rec  Name   Units  Quantity       Level               Description
% -------------------------------------------------------------------------
% 1  - TMP    K       Valid at P1   2 m above gnd       Temp.
% 2  - SPFH   kg/kg   Valid at P1   2 m above gnd       Specific humidity
% 3  - PRES   Pa      Valid at P1   surface             Pressure
% 4  - UGRD   m/s     Valid at P1   10 m above gnd      u wind
% 5  - VGRD   m/s     Valid at P1   10 m above gnd      v wind
% 6  - DLWRF  W/m^2   Valid at P1   surface             Downward long wave flux
% 7  - CLWMR  kg/kg   Accumulation  surface             Cloud water
% 8  - CAPE   J/kg    Valid at P1   180-0 mb above gnd  Convective Avail. Pot. Energy
% 9  - PEVAP  kg/m^2  Accumulation  surface              Pot. evaporation
% 10 - APCP   kg/m^2  Accumulation  surface              Total precipitation
% 11 - DSWRF  W/m^2   Valid at P1   surface              Downward short wave flux

% GLAS, description for netCDF4 format
% Short Name Description Unit
% Swnet_tavg Net short wave radiation flux W m-2
% Lwnet_tavg Net long-wave radiation flux W m-2
% Qle_tavg Latent heat net flux W m-2
% Qh_tavg Sensible heat net flux W m-2
% Qg_tavg Heat flux W m-2
% Snowf_tavg Snow precipitation rate kg m-2 s-1
% Rainf_tavg Rain precipitation rate kg m-2 s-1
% Evap_tavg Evapotranspiration kg m-2 s-1
% Qs_acc Storm surface runoff kg m-2
% Qsb_acc Baseflow-groundwater runoff kg m-2
% Qsm_acc Snow melt kg m-2
% AvgSurfT_inst Average Surface Skin temperature K
% Albedo_inst Albedo %
% SWE_inst Snow depth water equivalent kg m-2
% SnowDepth_inst Snow depth m
% SoilMoi0_10cm_inst Soil moisture kg m-2
% SoilMoi10_40cm_inst Soil moisture kg m-2
% SoilMoi40_100cm_inst Soil moisture kg m-2
% SoilMoi100_200cm_inst Soil moisture kg m-2
% SoilTMP0_10cm_inst Soil temperature K
% SoilTMP10_40cm_inst Soil temperature K
% SoilTMP40_100cm_inst Soil temperature K
% SoilTMP100_200cm_inst Soil temperature K
% PotEvap_tavg Potential evaporation rate W m-2
% ECanop_tavg Canopy water evaporation W m-2
% Tveg_tavg Transpiration W m-2
% ESoil_tavg Direct Evaporation from Bare Soil W m-2
% RootMoist_inst Root zone soil moisture kg m-2
% CanopInt_inst Plant canopy surface water kg m-2
% Wind_f_inst Wind speed m/s
% Rainf_f_tavg Total precipitation rate kg m-2 s-1
% Tair_f_inst Temperature K
% Qair_f_inst Specific humidity kg/kg
% Psurf_f_inst Pressure Pa
% SWdown_f_tavg Downward short-wave radiation flux W m-2
% LWdown_f_tavg Downward long-wave radiation flux W m-2
%
% The short names with extension “_tavg” are past 3-hr averaged variables.
% The short names with extension “_acc” are past 3-hr accumulated variables.
% The short names with extension “_inst” are instantaneous variables.
% The short names with “_f” are forcing variables.


nldas_flag=false;
gldas_flag=false;

if contains(filename,'NLDAS')
    nldas_flag=true;
elseif contains(filename,'GLDAS')
    gldas_flag=true;
end

%Read file
if nldas_flag
    grib_struct=read_grib(filename,{v},'ParamTable','NCEPREAN','ScreenDiag',0);
    if isempty(grib_struct)
        error('empty grib for file:%s variable: %s',filename,v);
    end
    %Get Matrix Dimensions, create referencing matrix and bounding box
    gds=grib_struct.gds;
    Ni=gds.Ni;
    Nj=gds.Nj;
    Di=gds.Di;
    Dj=gds.Dj;
    % Make referencing matrix for NLDAS, it appears rounded so adjust
    %http://www.emc.ncep.noaa.gov/mmb/nldas/LDAS8th/LDASspecs/LDASspecs.shtml
    %LDASconus.Rmap = makerefmat(-124.9375,52.9375, Dj, -Di);
    Rmap = makerefmat(gds.Lo1+0.0005,gds.La2-0.0005, Dj, -Di);
    %Reshape,rotate, and scale matrix
    dt=rot90(reshape(grib_struct.fltarray,[Ni Nj]));
    %9.9990e+20 is commonly the no data value
    max_dt=max(max(dt));
    if max_dt>0
        dt(dt==max_dt)=NaN;
    end
    matrix=dt;
    
    %Universal Time
    pds=grib_struct.pds;
    
    %Use 2 digit year to get 4 digit year, (hack?)
    if pds.year < 100 && pds.year > 50
        y1=19;
    elseif pds.year==100
        y1 = 20;
        pds.year=00;
    else
        y1=20;
    end
    yr=str2double([sprintf('%02s',num2str(y1)) ...
        sprintf('%02s',num2str(pds.year))]);
    
    % Add datevals to output structure
    dateval=datenum([yr pds.month pds.day pds.hour pds.min 0]);
    units=grib_struct.units;
    
elseif  gldas_flag
    ncid = netcdf.open(filename);
    varid = netcdf.inqVarID(ncid,v);
    x=netcdf.getVar(ncid,varid);
    x=rot90(x);
    nullval=-9999;
    x(x==nullval)=NaN;
    matrix=x;
    units=netcdf.getAtt(ncid,varid,'units');
    varid = netcdf.inqVarID(ncid,'lat');
    lat=netcdf.getVar(ncid,varid);
    varid = netcdf.inqVarID(ncid,'lon');
    lon=netcdf.getVar(ncid,varid);
    Rmap = makerefmat(double(lon(1)),double(lat(end)),double(lon(2)-lon(1)),...
        double(lat(1)-lat(2)));
    %http://ldas.gsfc.nasa.gov/gldas/data/0.25deg/lis_elev.ctl
    % -179.875 0.25; 59.875 -0.25 (starting from upper left corner)
    %
    %begin_date time attribute eliminated at some point in 2019
    % now use filename for this, NB 12/18/20
    idx=regexp(filename,'[0-9]{8}\.[0-9]{4}');
    dateval=datenum(filename(idx:(idx+12)),'yyyymmdd.HHMM');
%     varid=netcdf.inqVarID(ncid,'time');
%     bd=netcdf.getAtt(ncid,varid,'begin_date');
%     bt=netcdf.getAtt(ncid,varid,'begin_time');
%     dateval=datenum([bd,bt],'yyyymmddHHMMSS');
    netcdf.close(ncid);
end

LDAS.matrix=matrix;
LDAS.dateval=dateval;
LDAS.units=units;
LDAS.Rmap=Rmap;

LDAS=orderfields(LDAS);