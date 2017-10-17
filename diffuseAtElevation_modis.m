function FinZ = diffuseAtElevation_modis(Fin,tau,pres,mu0,tauZ,presZ)
% Dubayah and Loechel (1997). J. Applied Met. Appendix C
% INPUT
% Fin - diffuse incoming at reference elevation
% tau -optical depth at reference elevation
% tauZ - optical depth, adjusted to fine scale elevation
% pres - air pressure at reference elevation
% presZ - air pressure, adjusted to fine scale elevation
% mu0 - fine scale cos solar zenith
%OUTPUT
% FinZ - diffuse incoming corrected for elevation

% Fraction of unabsorbed exoatmospheric flux
M0 = fracUnAbsExoAtmo(pres,mu0);
MZ = fracUnAbsExoAtmo(presZ,mu0);

% ratio of elevation diffuse to reference diffuse
r0 = M0 - exp(-tau./mu0);
rz = MZ - exp(-tauZ./mu0);
ratio = rz./r0;
% Correct diffuse for elevation
ratio(ratio<0.7)=0.7;
ratio(ratio>1.2)=1.2;
FinZ=Fin.*ratio;

% Set diffuse at elevation equal to spatially integrated diffuse where
% where r close to zero (small changes result in huge corrections)
% FinZ(rz<0.05)=Fin(rz<0.05);
% FinZ(r0<0.05)=Fin(r0<0.05);