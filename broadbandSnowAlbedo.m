function albedo = broadbandSnowAlbedo(radius,Zenith,varargin)
% albedo = broadbandSnowAlbedo(radius,Zenith,{carbon,tau})
% broadband albedo as a function of grain radius and solar zenith angle Z
% (optionally including carbon concentration and cloud optical thickness
%
% Input
%   radius - grain radius, mm
%   Zenith - illumination angle, degrees
% Optional input
%   carbon - ppmw 0-2
%   tau - cloud optical thickness, 0-30
%
% Output
%   albedo - broadband
%
% Adapted from
% Gardner, A. S., and M. J. Sharp (2010), A review of snow and ice albedo
% and the development of a new physically based broadband albedo
% parameterization, Journal of Geophysical Research, 115, F01009,
% doi: 10.1029/2009jf001444.
% Comments refer to the equations (numbered) in that paper

% check range of inputs
assert(~any(radius(:)<0.009) && ~any(radius(:)>=10),...
    'radius in mm, not um (your min=%g, max=%g)',...
    nanmin(radius(:)),nanmax(radius(:)))
assert(~any(Zenith(:)<0) && ~any(Zenith(:)>90),...
    'arg Zenith out of range [%g %g]', [nanmin(Zenith(:)) nanmax(Zenith(:))])

% initial values for optional arguments
carbon = 0;
tau = 0;
dac = 0;
dcloud = 0;

% Grid of solar zenith, or scalar
if ~xor(isscalar(Zenith),isscalar(radius))
    assert(isequal(size(radius),size(Zenith)),...
        'if not scalars, radius and Zenith arrays must be same size')
end

% radius where there is no snow is NaN, not zero, so reset
radius(radius==0) = NaN;

for k=1:size(varargin,2)
    switch k
        case 1
            carbon = varargin{k};
            if ~isscalar(carbon) && ~(isscalar(radius) && isscalar(Zenith))
                assert(isequal(size(carbon),size(radius)) ||...
                    isequal(size(carbon),size(zenith)),...
                    'carbon dataset must be size of radius or zenith')
            end
        case 2
            tau = varargin{k};
            if ~isscalar(tau) && ~(isscalar(radius) && isscalar(Zenith))
                assert(isequal(size(tau),size(radius)) ||...
                    isequal(size(tau),size(zenith)),...
                    'tau dataset must be size of radius or zenith')
            end
    end
end

% ssa in Gardner & Sharp is in cm^2/g
radius = radius/10; % convert grain size to cm
rhoIce = .917; % ice density
ssa = 3./(rhoIce*radius);

% equation 7 - clean snow at 0 degree Z
as = 1.48-ssa.^(-0.07);

% equation 8 - modification because of carbon (can mimic dust?)
if carbon>0
    dac = max(0.04-as, -carbon.^0.55./(0.16+0.6*sqrt(ssa)+1.8*carbon.^0.6.*ssa.^(-.25)));
end
ac = as+dac;

% equation 9 - modification for solar zenith
mu0 = cosd(Zenith);
dtheta = 0.53*as.*(1-ac).*(1-mu0).^1.2;

% equations 10/11 - modification for cloud thickness
if tau>0
    x = min(sqrt(tau/(3*mu0)),1);
    uprime = 0.64*x+(1-x)*mu0;
    asd = as+0.53*as*(1-ac)*(1-uprime)^1.2;
    dcloud = 0.1*tau*ac^1.3/((1.5*tau)^asd);
end

% equation 12 - final result
albedo = as+dac+dtheta+dcloud;
% for really really small grains, result too large, so correct
tooBright = 0.96;
if nnz(albedo>tooBright)
    albedo(albedo>tooBright) = tooBright-rand(nnz(albedo>tooBright),1)*0.02;
end

end