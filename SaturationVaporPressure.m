function [ vp ] = SaturationVaporPressure( T, substance )
%SatVP saturation vapor pressure over water or ice
% [ vp ] = SaturationVaporPressure( T, substance )
%
% input
%   T - Kelvin temperature, can be a vector or matrix
%   substance - either 'ice' or 'water'
%
% output
%   vp - saturation vapor pressure(s) in kPa corresponding to temperature(s)
%   (set to NaN where substance is 'ice' and T>273.16
%
% Equations are from Bohren, C. F., and B. A. Albrecht (1998),
% Atmospheric Thermodynamics, Oxford University Press.


assert(strcmpi(substance,'ice') || strcmpi(substance,'water'),...
    'substance (2nd argument) must be either ''ice'' or ''water''')
assert(single(nnz(T(T<0))) == 0, 'Temperatures must be in Kelvin')

T0 = 273.16;
es0 = 0.611; % sat vp in kPa at 273.16K
if strcmpi(substance,'ice')
    logs = 6293.*(1/T0 - 1./T) - 0.555*log(T/T0);
    ok = T <= T0;
    if single(nnz(~ok))>0
        logs(~ok) = NaN;
    end
else % water
    logs = 6808.*(1/T0 - 1./T) - 5.09*log(T/T0);
end
vp = es0 * exp(logs);

end

