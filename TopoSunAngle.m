function T  = TopoSunAngle(matdate,topo,press,temp)
%computes sun angles over a topographic grid
%INPUT -
%   matdate - MATLAB date in UTC over which to do the computation
%   topo -  topo structure including filename for slope,aspect,& horizon in hd5
%   Press - grid of pressures or scalar (Pa)
%   Temp - grid of temperatures or scalar (K)
%
%OUTPUT -
%   T - topographic structure with elements over the grid
%       scalars
%           matdate
%           Declination
%           RadiusVector
%           SolarLongitude
%       grids
%           mu0 - cosine of sun angle on horizontal surface
%           mu0_unrefracted -top of atmosphere " ", not corr. for refract.
%           phi0 - sun azimuth, degrees, +ccw from south
%           airmass - 1.0 means sun overhead, others larger
%           mu - cosine of sun angle on slope (set to zero if shaded)
S=topo.hdr;
T.gridtype = S.gridtype;
if isfield(S,'RefMatrix')
    T.RefMatrix = S.RefMatrix;
end
if isfield(S,'Projection')
    T.Projection = S.Projection;
elseif isfield(S,'Geoid')
    T.Geoid = S.Geoid;
end

% sun location
T.product = 'sun angles';
T.date = matdate;

% [T.Declination,T.RadiusVector,T.SolarLongitude] = Ephemeris(matdate);
[T.Declination,T.RadiusVector,T.SolarLongitude] = EarthEphemeris(matdate);

% lat-lon if geolocated
switch S.gridtype
    case 'geolocated'
        T.Lat = S.Lat;
        T.Lon = S.Lon;
        lat = S.Lat;
        lon = S.Lon;
    case 'projected'
        [x,y] = pixcenters(S.RefMatrix,size(topo.slope),'makegrid');
        [lat,lon] = minvtran(S.ProjectionStructure,x,y);
    case 'geographic'
        [lon,lat] = pixcenters(S.RefMatrix,size(topo.slope),'makegrid');
    otherwise
        error('S.gridtype unknown')
end

% sun angles w/ amospheric refraction
[mu0,mu0_unrefracted, phi0, airmass] = sunang2(lat,lon,T.Declination,T.SolarLongitude,...
        true,press/1000,temp);
% sun on slopes
mu = sunslope(mu0,phi0,topo.slope,topo.aspect);
%mask for slopes above the horizon being illuminated
hmask = GetHorizon(topo.topofile,phi0,mu0);
%set slopes below horizon to 0
mu(~hmask) = 0;

% put into structure
T.mu0 = mu0;
T.mu0_unrefracted = mu0_unrefracted;
T.phi0 = phi0;
T.airmass = airmass;
T.mu = mu;

end