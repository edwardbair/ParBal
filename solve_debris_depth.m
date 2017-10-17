function [d,LinZ,Lout,sensible,G]=solve_debris_depth(Lin_coarse,Zdiff,T_fine,Vf,...
    albedo,Sin,Tsfc,sig,ems,rho_air,Cp,De_h,K)
%solve for dry debris depth d at one time using inputs:
%Lin_coarse, coarse  incoming longwave, W/m^2
%Zdiff, difference between fine and coarse elevation, corresponding to
%Lin_coarse
%T_fine, fine air temp, K
%Vf, view factor, 0-1
%albedo, debris albedo, 0-1
%pres_fine, fine_scale atmos pressure, kPa
%ea, vapor pressure of air
%Tsfc, debris surface temperature
%sig, steffan boltzman constant
%ems, emissivity of debris (i.e. 0.94 for metamorphic rock &
%granite)
%Cp, specific heat of air, J/(kg deg)
%De_h, sensible heat exchange coefficient, m/sec
%K, conducitivity of debris, W/(m*K), e.g. 1 W/(m*K) is mean from
%Schauwecker et al 2015 Table 3
%output
%d - debris depth, m

LinZ=lapseLongWave(Lin_coarse,Zdiff);
% adjust for surrounding terrain
Lout = -(ems.*sig.*Tsfc.^4);
LinT = Vf.* LinZ + (1-Vf).*Lout;
Lnet=LinT+Lout;
Snet=Sin*(1-albedo);
% assume a neutral atmosphere (stability = 1). See Rounce and McKiney 2014
stability = single(1);
sensible=rho_air.*Cp.*De_h.*stability.*(T_fine-Tsfc);
% debris is dry during Tsfc aquisitions so latent = 0.
% use iterative approach of Schauwecker et al 2015
num_tries = 20;
d0 = 0.45; %m, initial guess
thresh = 1e-3;
for i=1:num_tries
    F = 6.17*d0+1;
    d = ((1+F).*K.*(Tsfc-273.15))./(Snet+Lnet+sensible).*1/0.5;
    if abs(d-d0) < thresh
        break;
    end
    d0 = d;
end
%if d is negative or didn't converge
if d < 0 || i==num_tries 
    d=NaN;
end
%Schauwecker et al 2016
G=(K.*(Tsfc-273.15))/(0.5*d);