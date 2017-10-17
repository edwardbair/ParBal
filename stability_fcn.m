function [stability,De_h]=stability_fcn(windspd,Tair,Tsfc)
stability=ones(size(Tair));
%exchange coeffs for sensible and latent heat
% Liston (1995) Local Advection of Momentum, 
% Heat and Moisture during the Melt of Patchy Snow Covers. JAM
% and Price and Dunne (1976)

%constants
ht_wind_obs=10; %m
xkappa=0.4;% von Karmans constant
gravity=9.8; %m/s
z_0=0.0005; %snow roughness length, m

C1=5.3.*9.4.*(xkappa./(log(ht_wind_obs./z_0))).^2.*...
    sqrt(ht_wind_obs./z_0);
C2=gravity.*ht_wind_obs./(Tair.*windspd.^2);
B1=9.4.*C2;
B2=C1.*sqrt(C2);

%unstable
ind=Tsfc > Tair;
if any(ind(:));
    B3= 1 + B2(ind).*sqrt(Tsfc(ind)-Tair(ind));
    stability(ind) = 1 + B1(ind).*(Tsfc(ind)-Tair(ind))./B3;
end
%stable
ind2=Tsfc < Tair;
if any(ind2(:));
    B8=B1(ind2)./2;
    stability(ind2) = 1./((1+B8.*(Tair(ind2)-Tsfc(ind2))).^2);
end
% else neutrally stable = 1
%compute turbulent coef
De_h=(xkappa.^2.*windspd)./((log(ht_wind_obs./z_0)).^2);