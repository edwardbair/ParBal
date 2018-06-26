function a=scagd_albedo(r,cosZ)
%statsitical fit to Jeff Dozier's SCAGD albedo model, based of G-S approach
%input r - grain radius, um
%cosZ - cosine of the solar zenith angle

%albedo as function of grain radius at cosZ=1, r^2=1.00
a=-0.1745;
b=0.1456;
c=1.131;
as=a.*r.^b+c;


x=as;y=cosZ;

%4th degree Horner polynomial

p00=-0.9757;
p01=-2.6268;
p02=2.7072;
p03=-1.9832;
p04=0.5199;
p10=8.3210;
p11=4.9436;
p12=-1.9170;
p13=0.8584;
p20=-19.1899;
p21=-4.4519;
p22=0.0366;
p30=18.3292;
p31=1.8901;
p40=-6.4690;

dtheta=p00 + p10.*x + p01.*y + y.^2.*(p02 + y.*(p03 + p04.*y)) ...
    + x.*(y.*(p11 + y.*(p12 + p13.*y)) + x.*(p20 + y.*(p21 + p22.*y) ...
    + x.*(p30 + p40.*x + p31.*y)));

a=as+dtheta;
                
