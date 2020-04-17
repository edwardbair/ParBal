<<<<<<< HEAD
function a=scagd_albedo(r,cosZ)
%statsitical fit to Jeff Dozier's SCAGD albedo model, based off G-S approach
%input r - grain radius, um
%cosZ - cosine of the solar zenith angle

%albedo as function of grain radius at cosZ=1, r^2=1.00
% a=-0.1745;
% b=0.1456;
% c=1.131;
a=-0.1963;
b=0.1312;
c=1.144;

%use as as input into polynomial
as=a.*r.^b+c;
x=as;
%and cosZ
y=cosZ;
%5th degree Horner polynomial

p00 =       1.277;
p01 =       6.839;
p02 =      -1.798;
p03 =       3.181;
p04 =      -1.376;
p05 =      0.2043;
p10 =         -17;
p11 =      -34.81;
p12 =      -1.163;
p13 =      -3.613;
p14 =      0.9432;
p20 =       69.45;
p21 =       71.91;
p22 =       6.238;
p23 =       0.716;
p30 =      -125.7;
p31 =      -68.02;
p32 =      -3.348;
p40 =       106.4;
p41 =       24.14;
p50 =      -34.49;

dtheta=p00 + p01.*y + p02.*y.^2 + p03.*y.^3 + p04.*y.^4 + p05.*y.^5 + x.*(p10 + ...
    x.*(p20 + p21.*y + p22.*y.^2 + p23.*y.^3 + x.*(p30 + x.*(p40 ...
    + p50.*x + p41.*y) + p31.*y + p32.*y.^2)) + p11.*y + p12.*y.^2 ...
    + p13.*y.^3 + p14.*y.^4);
a=as+dtheta;

=======
function a=scagd_albedo(r,mu0,varargin)
%statsitical fit to Jeff Dozier's SCAGD, which is  a MATLAB version of
%the Warren and Wiscombe (1980) albedo model
%input r - grain radius (scalar, vector, or matrix), um
%mu0 - cosine of the solar zenith angle (scalar, vector, or matrix), 0-1
%optional - 'mlw' (mid-lat. winter @ 3km elevation, default) or
%'sas' (sub-Arctic summer @ sea level) for atmospheric profile
%impurities to come later
%see Bair et al (in review) WRR, submitted 1/2019
% Ned Bair 2/2019

p = inputParser;
addRequired(p,'r',@(x) all((x >= 1 & x < 2000) | isnan(x),'all'))
addRequired(p,'mu0',@(x) all(x >= 0 & x <= 1 | isnan(x),'all'))
defaultatm='mlw';
addOptional(p,'atm',defaultatm,@(s) ischar(s));
parse(p,r,mu0,varargin{:});
R=p.Results.r;
M=p.Results.mu0;
atm=p.Results.atm;

assert(strcmp(atm,'mlw') || strcmp(atm,'sas'),...
    'atmosphere has to be ''mlw'' or ''sas''');
% cutoff albedo at 85 deg/cosd(85)=0.09
M(M<0.09)=0.09;

switch atm
    case 'mlw'
        P=[-9.025001 -6.853901 -6.360441;...
            0.05785986 0.273218 0.1890732;...
            0.07632736 1.017243 0.4149719];
        Q=[1 92.35081 27.87415;...
           1 1.28665 1.53981;...
           0 1 0.3373872];
    case 'sas'
        P=[-0.6458545 -0.1641362 -0.4793498;...
            0.1143997 0.05545726 0.1315713;...
            0.06805254 0.991294 0.5284415];
        Q=[1 5.093014 2.773746;...
           1 0.0001412588 0.9561747;...
           0 1 0.446777];
end

A=(P(1,1).*M.^2+P(1,2).*M+P(1,3))./(M.^2+Q(1,2).*M+Q(1,3));
B=(P(2,1).*M.^2+P(2,2).*M+P(2,3))./(M.^2+Q(2,2).*M+Q(2,3));
D=(P(3,1).*M.^2+P(3,2).*M+P(3,3))./(M+Q(3,3));
a=A.*R.^B+D;
>>>>>>> a551fb190757d0789821ebee09340b8365d27973
