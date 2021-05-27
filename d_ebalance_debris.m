function dM=d_ebalance_debris(Tsfc,T_fine,Vf,sig,ems,...
    Cp,rho_air,De_h,B1,B2,K,d)
%analytical deriative of L
dL=- 4*Tsfc^3*ems*sig - 4*Tsfc^3*ems*sig*(Vf - 1) - ...
    4*Tsfc^3*ems*sig*(Vf - 1)*(ems - 1);
%analytical derivative of sensible heat given unstable
%conditions (Tsfc > Tair)
if Tsfc > T_fine
dS = Cp*De_h*rho_air*((B1*(T_fine - Tsfc))/...
    (B2*(Tsfc - T_fine)^(1/2) + 1) - 1) + Cp*De_h*rho_air*(T_fine - Tsfc)*...
    (B1/(B2*(Tsfc - T_fine)^(1/2) + 1) - (B1*B2*(Tsfc - T_fine)^(1/2))/...
    (2*(B2*(Tsfc - T_fine)^(1/2) + 1)^2));
    %analytical derivative of sensible heat given stable conditions
elseif Tsfc < T_fine
    dS = (B1*Cp*De_h*rho_air*(T_fine - Tsfc))/((B1*(T_fine - Tsfc))/2 + 1)^3 - ...
        (Cp*De_h*rho_air)/((B1*(T_fine - Tsfc))/2 + 1)^2;
%     analytical derivative of sensible heat given neutral conditions
else
    dS=-Cp*De_h*rho_air;
end
%analytical derivative of G
dG = -K/d ;
%total melt derivative
dM=dL+dS+dG;
end