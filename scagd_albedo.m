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

