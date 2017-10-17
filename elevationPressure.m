function presZ  = elevationPressure(pres0f,l_rate,Zdiff,tmp0f)
% Estimate pressure at elevation and ratio to reference pressure using
% hydrostatic equation
%INPUT
% pres0f - pressure at reference elevation
% l_rate - lapse rate
% Zdiff - elevation difference (fine - coarse)
% tmp0f - temperature at reference
%OUTPUT
% presZ - pressure at elevation, Pa

%Constants 
R=8.31432;%Universal gas constant
g0=9.80665;%Gravitational acceleration
M=0.0289644;%Molecular weight of air

%Exponent
pow = (-g0*M)/(l_rate*R);

%Base
base = 1 + (l_rate.*Zdiff)./tmp0f;

%Pressure ratio Pz/Pz0
pRatio = base.^pow;

%Fine resolution pressure
presZ = pRatio.*pres0f;