function [d,LinZ,Lout,sensible,G]=solve_debris_depth(Lin_coarse,Zdiff,T_fine,Vf,...
    albedo,Sin,Tsfc,sig,emd,rho_air,Cp,De_h,K)
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
%emd, emissivity of debris (i.e. 0.94 for metamorphic rock &
%granite)
%Cp, specific heat of air, J/(kg deg)
%De_h, sensible heat exchange coefficient, m/sec
%K, conducitivity of debris, W/(m*K), e.g. 1 W/(m*K) is mean from
%Schauwecker et al 2015 Table 3
%output
%d - debris depth, m

LinZ=lapseLongWave(Lin_coarse,Zdiff);
% adjust for surrounding terrain
Lout = -(emd.*sig.*Tsfc.^4);
LinT = Vf.* LinZ + (1-Vf).*-Lout;
Lnet=LinT+Lout;
Snet=Sin*(1-albedo);
% assume a neutral atmosphere (stability = 1). See Rounce and McKiney 2014
stability = single(1);
sensible=rho_air.*Cp.*De_h.*stability.*(T_fine-Tsfc);
% debris is dry during Tsfc aquisitions so latent = 0.

Tsfc_c=Tsfc-273.15;

%Schauwecker approach
% a=6.17;b=1;
% F=@(d0) a*d0*b;
% i_d=0.5;
% ddepth=@(d0)(((1+F(d0))*K*Tsfc_c)/(Snet+Lnet+sensible))*(1/i_d);
% minfct=@(d0) abs(d0-ddepth(d0));
% options = optimset('Display','off');
% lo=-2;
% hi=2;
% x=fminbnd(minfct,lo,hi,options);
% d=x;

%Rounce and McKinney approach
Gratio=2.7;
d=(Gratio*K*Tsfc_c)/(Snet+Lnet+sensible);
G=(K.*(Tsfc-273.15))/(d);