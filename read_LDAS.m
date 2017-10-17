function LDAS = read_LDAS(v,filename)
%LDAS = read_LDAS(v,filename)
% Read variable from GRIB NLDAS or GLDAS,requires read_grib package
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

%GLDAS
% Average - 3hr average 
% Rec  Name   Units  Quantity     Level        Description
% -------------------------------------------------------------------------
% 1  - NSWRS  W/m^2  Average      surface      Net short wave (surface) 
% 2  - NLWRS  W/m^2  Average      surface      Net long wave (surface) 
% 3  - LHTFL  W/m^2  Average      surface      Latent heat flux 
% 4  - SHTFL  W/m^2  Average      surface      Sensible heat flux 
% 5  - GFLUX  W/m^2  Average      surface      Ground heat flux 
% 6  - LFTX   K      Valid at P1  surface      Surface lifted index 
% 7  - 4LFTX  K      Valid at P1  surface      Best (4-layer) lifted index 
% 8  - EVP    kg/m^2 Valid at P1  surface      Evaporation 
% 9  - SSRUN  kg/m^2 Valid at P1  surface      Storm surface runoff 
% 10 - BGRUN  kg/m^2 Valid at P1  surface      Baseflow-groundwater runoff 
% 11 - SNOM   kg/m^2 Average      surface      Snow melt 
% 12 - BVF2   1/s^2  Valid at P1  surface      Brunt-Vaisala frequency^2 
% 13 - WEASD  kg/m^2 Valid at P1  surface      Accum. snow 
% 14 - TSOIL  K      Valid at P1  0-4 cm down  Soil temp. 
% 15 - TSOIL  K      Valid at P1  0-3 cm down  Soil temp. 
% 16 - TSOIL  K      Valid at P1  0-2 cm down  Soil temp. 
% 17 - TSOIL  K      Valid at P1  0-1 cm down  Soil temp. 
% 18 - SOILM  kg/m^2 Valid at P1  0-4 cm down  Soil moisture content 
% 19 - SOILM  kg/m^2 Valid at P1  0-3 cm down  Soil moisture content 
% 20 - SOILM  kg/m^2 Valid at P1  0-2 cm down  Soil moisture content 
% 21 - SOILM  kg/m^2 Valid at P1  0-1 cm down  Soil moisture content 
% 22 - TCDC   %      Valid at P1  surface      Total cloud cover 
% 23 - WIND   m/s    Valid at P1  surface      Wind speed 
% 24 - TMP    K      Valid at P1  surface      Temp. 
% 25 - SPFH   kg/kg  Valid at P1  surface      Specific humidity 
% 26 - PRES   Pa     Valid at P1  surface      Pressure 
% 27 - DSWRF  W/m^2  Average      surface      Downward short wave flux 
% 28 - DLWRF  W/m^2  Average      surface      Downward long wave flux

%Read file
grib_struct=read_grib(filename,{v});
if isempty(grib_struct)
    error('empty grib for file:%s variable: %s',filename,v);
end
%Get Matrix Dimensions, create referencing matrix and bounding box
gds=grib_struct.gds;
Ni=gds.Ni;
Nj=gds.Nj;
Di=gds.Di;
Dj=gds.Dj;

% Check filename to determine if NLDAS or GLDAS, adjust NLDAS Rmap
if ~isempty(strfind(filename,'NLDAS'));
    % Make referencing matrix for NLDAS, it appears rounded so adjust
    %http://www.emc.ncep.noaa.gov/mmb/nldas/LDAS8th/LDASspecs/LDASspecs.shtml
    %LDASconus.Rmap = makerefmat(-124.9375,52.9375, Dj, -Di);
    LDAS.Rmap = makerefmat(gds.Lo1+0.0005,gds.La2-0.0005, Dj, -Di);
elseif  ~isempty(strfind(filename,'GLDAS'));
    LDAS.Rmap = makerefmat(gds.Lo1,gds.La2, Dj, -Di);
    %http://ldas.gsfc.nasa.gov/gldas/data/0.25deg/lis_elev.ctl
    % -179.875 0.25; -59.875 0.25 (left,lower)
end

%Reshape,rotate, and scale matrix
dt=rot90(reshape(grib_struct.fltarray,[Ni Nj]));

%9.9990e+20 is commonly the no data value but does not appear anywhere
max_dt=max(max(dt));
if max_dt>0
    dt(dt==max_dt)=NaN;
end

%Fill date cube with all x,y, for each hour
LDAS.matrix=dt;

%Universal Time
pds=grib_struct.pds;

%Use 2 digit year to get 4 digit year, (hack?)
if pds.year < 100 && pds.year > 50
    y1=19;
elseif pds.year==100
    y1 = 20;
    pds.year=00;
else y1=20;
end
yr=str2double([sprintf('%02s',num2str(y1)) ...
    sprintf('%02s',num2str(pds.year))]);

% Add datevals to output structure
LDAS.dateval=datenum([yr pds.month pds.day pds.hour pds.min 0]);

% Add units to output structure
LDAS.units=grib_struct.units;
LDAS=orderfields(LDAS);