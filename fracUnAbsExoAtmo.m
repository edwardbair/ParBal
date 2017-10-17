function Mj = fracUnAbsExoAtmo(Pj,mu0)
%fraction of unabsorbed exoatmospheric flux
%Dubayah and Loechel (1997). J. Applied Met. Appendix C
%INPUT
% Pj - pressure at height j (kPa)
% mu0 - cosine of solar zenith for flat surface
%OUTPUT
% Mj - fraction of unabsorbed exoatmospheric flux at height j

% P0 - pressure at sea level
P0=101.325; %kPa

% fraction of unabsorbed exoatmospheric flux at height j
Mj=(1-0.027.*exp(2.*Pj./P0)).*(1.075-0.105.*log(1./mu0));
