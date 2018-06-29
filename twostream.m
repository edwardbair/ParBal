function [Refl, Trans, Beam, method] = twostream( mu0, omega, g, tau, varargin )
% [Refl, Trans, Beam, method] = twostream( mu0, omega, g, tau, ...)
%twostream: Two-stream solution of radiative transfer for single layer
%
% Input (omega, g & tau must be same size, mu0 can be scalar or same size)
%   mu0 - cosine of illumination angle
%   omega - single-scattering albedo
%   g - asymmetry factor
%   tau - optical depth
% Optional input - name value pairs
%   'R0', then number - reflectance of substrate, set to zero if omitted
%       (scalar or same size as omega, g, tau)
%   'airmass', then number - compensated for Earth curvature, set to 1/mu0
%       if omitted, otherwise must be same size as mu0
%   'method' - either 'delta-Eddington' or 'hybrid' (default)
%
% Output
%   Refl - Reflectance of layer, including substrate
%   Trans - Transmittance of layer
%   Beam - Direct transmittance of layer (i.e. to beam illumination)
%   method - calculation method used
%
% source: Meador, W. E., and W. R. Weaver (1980), Two-stream approximations
% to radiative transfer in planetary atmospheres – A unified description of
% existing methods and a new improvement, Journal of the Atmospheric
% Sciences, 37, 630-643, doi:
% 10.1175/1520-0469(1980)037<0630:TSATRT>2.0.CO;2.

% inputs
p = inputParser;
addRequired(p,'mu0',@isnumeric)
addRequired(p,'omega',@isnumeric)
addRequired(p,'g',@isnumeric)
addRequired(p,'tau',@isnumeric)
defaultR0 = 0;
defaultAirmass = [];
defaultMethod = 'hybrid';
addParameter(p,'R0',defaultR0,@isnumeric)
addParameter(p,'airmass',defaultAirmass,@isnumeric)
addParameter(p,'method',defaultMethod,@ischar)
parse(p,mu0, omega, g, tau, varargin{:})

% check inputs
assert(isequal(size(omega),size(g),size(tau)),...
    'omega, g, and tau must be same size')
if ~isscalar(mu0)
    assert(isequal(size(omega),size(mu0)),...
        'if not scalar, mu0 must be same size as omega, g, & tau')
else
    mu0 = repmat(mu0,size(omega));
end

% to run with matrix inputs, save size and convert to vectors, then resize
% at end
origSize = size(omega);
omega = omega(:);
g = g(:);
tau = tau(:);
mu0 = mu0(:);

assert(all(mu0>=0 & mu0<=1),'some mu0[%g %g]  out of range, must be 0-1',...
    min(mu0(:)),max(mu0(:)))
assert(all(g>=0 & g<=1),'some g [%g %g] out of range, must be 0-1',...
    min(g(:)),max(g(:)))
assert(all(omega>=0 & omega<=1),'some omega [%g %g] out of range, must be 0-1',...
    min(omega(:)),max(omega(:)))
assert(all(tau>0),'tau=%g out of range, must be > 0', min(tau(:)))
% optional inputs - set defaults
if isempty(p.Results.airmass)
    airmass = 1./mu0;
else
    [airmass,~] = checkSizes(p.Results.airmass,omega);
    assert(all(airmass>0),'airmass=% out of range, must be >0',min(airmass))
end
[R0,~] = checkSizes(p.Results.R0(:),omega);
assert(all(R0>=0 & R0<=1),'R0 [%g %g] out of range, must be 0-1',...
    min(R0(:)),max(R0(:)))
method = p.Results.method;
assert(strcmp(method,'hybrid') || strcmp(method,'delta-Eddington'),...
    '''method'' must be ''hybrid'' or ''delta-Eddington''')

% Meador-Weaver gamma values
[gam, method] = mwgamma( mu0, omega, g, method);

% hack: gam(1) = gam(2) with conservative scattering
t = gam(:,1)<=gam(:,2);
if any(t)
    gam(t,2) = gam(t,1)*(1-eps);
end
% intermediate variables
alph1 = gam(:,1) .* gam(:,4) + gam(:,2) .* gam(:,3);
alph2 = gam(:,2) .* gam(:,4) + gam(:,1) .* gam(:,3);
xi = sqrt((gam(:,1) - gam(:,2)) .* (gam(:,2) + gam(:,1)));
em = exp(-tau .* xi);
et = exp(-tau .* airmass);
ep = exp(tau .* xi);
gpx = xi + gam(:,1);
opx = mu0 .* xi + 1;

Refl = zeros(size(omega));
Trans = zeros(size(omega));
Beam = zeros(size(omega));
% t = et==0 | em==0 | isinf(ep); % semi-infinite
t = isinf(ep);
if any(t)
    Refl(omega==1) = 1;
    tx = omega~=1;
    if any(tx)
        Refl(tx) = omega(tx) .* (gam(tx,3) .* xi(tx) + alph2(tx)) ./...
            (gpx(tx) .* opx(tx));
    end
end
if all(~t)
    % more intermediate variables, needed only for finite case
    omx = 1 - mu0 .* xi;
    gmx = gam(:,1) - xi;
    rm = gam(:,2) - gmx .* R0;
    rp = gam(:,2) - gpx .* R0;
    
    % denominator for reflectance and transmittance
    denrt = ep .* gpx .* rm - em .* gmx .* rp;
    
    % reflectance
    Refl = (omega.*(ep.*rm.*(gam(:,3).*xi+alph2)./opx -...
        em.*rp.*(alph2-gam(:,3).*xi)./omx) + 2*et.*gam(:,2).*...
        (R0-((alph1.*R0-alph2).*mu0+gam(:,4).*R0+gam(:,3)).*omega./...
        (omx.*opx)).*xi)./denrt;
    
    % transmittance
    Trans = (et.*(ep.*gpx.*(gam(:,2)-omega.*(alph2-gam(:,3).*xi)./omx)-...
        em.*gmx.*(gam(:,2)-omega.*(gam(:,3).*xi+alph2)./opx))+...
        2*gam(:,2).*(alph1.*mu0+gam(:,4)).*omega.*xi./(omx.*opx))./denrt;
    Trans(Trans>1) = 1-eps;
    
    % direct transmittance
    Beam = et;
elseif any(~t)
    mu0 = mu0(~t);
    xi = xi(~t);
    omx = 1-mu0.*xi;
    gmx = gam(~t,1)-xi;
    R0 = R0(~t);
    gpx = gpx(~t);
    rm = gam(~t,2)-gmx.*R0;
    rp = gam(~t,2)-gpx.*R0;
    
    % denominator for reflectance and transmittance
    ep = ep(~t);
    em = em(~t);
    denrt = ep .* gpx .* rm - em .* gmx .* rp;
    
    % reflectance
    omega = omega(~t);
    alph1 = alph1(~t);
    alph2 = alph2(~t);
    opx = opx(~t);
    et = et(~t);
    Refl(~t) = (omega.*(ep.*rm.*(gam(~t,3).*xi+alph2)./opx -...
        em.*rp.*(alph2-gam(~t,3).*xi)./omx) + 2*et.*gam(~t,2).*...
        (R0-((alph1.*R0-alph2).*mu0+gam(~t,4).*R0+gam(~t,3)).*omega./...
        (omx.*opx)).*xi)./denrt;
    
    % transmittance
    Trans(~t) = (et.*(ep.*gpx.*(gam(~t,2)-omega.*(alph2-gam(~t,3).*xi)./omx)-...
        em.*gmx.*(gam(~t,2)-omega.*(gam(~t,3).*xi+alph2)./opx))+...
        2*gam(~t,2).*(alph1.*mu0+gam(~t,4)).*omega.*xi./(omx.*opx))./denrt;
    Trans(Trans>1) = 1-eps;
end

if ~isequal(size(Refl),origSize)
    Refl = reshape(Refl,origSize);
    Trans = reshape(Trans,origSize);
    Beam = reshape(Beam,origSize);
end

end

function [ gam, method ] = mwgamma( mu0, omega, g, method )
% [ gam, method ] = mwgamma( mu0, omega, g, method )
%mwgamma - Meador-Weaver gamma values
%
% call as
%   [gam, method] = mwgamma( mu0, omega, g, method) where method is
%   optional
%
% Input
%   mu0 - cosine of illumination angle
%   omega - single-scattering albedo
%   g - asymmetry factor
%   method  - either 'delta-Eddington' or 'hybrid' (default)
%
% Output
%   gamma - values as a column vector
%   method - calculation method used

deltaEdd = strcmpi(method,'delta-Eddington'); %set flag to zero if case

% check input ranges
assert(isequal(size(mu0),size(omega),size(g)),...
    'input mu0, omega, g must be same size')
assert(all(mu0>=0 & mu0<=1), 'mu0 must be 0-1')
assert(all(omega>=0 & omega<=1), 'omega must be 0-1')
assert(all(g>=0 & g<=1), 'g=%g, must be 0-1')

gam = zeros(length(g),4);
b0 = betanaught(mu0, g);

if deltaEdd
    t = omega==1;
    if any(t)
        gam(t,1) = 3*(1-g(t))/4;
        gam(t,2) = gam(1);
        gam(~t,1) = (7-(3*g(~t)+4).*omega(~t))/4;
        gam(~t,2) = ((4-3*g(~t)).*omega(~t)-1)/4;
    else
        gam(:,1) = (7-(3*g+4).*omega)/4;
        gam(:,2) = ((4-3*g).*omega-1)/4;
    end
else
    method = 'Meador-Weaver hybrid';
    hd = 4*(1-g.*g.*(1-mu0)); % denominator for both gam(1:2)
    
    % Horner expressions for Meador-Weaver Table1
    t = omega==1;
    if any(t)
        gam(t,1) = (g(t).*((3*(g(t)-1)+4*b0(t)).*g(t)-3)+3)./hd(t);
        gam(t,2) = gam(t,1);
        gam(~t,1) = (g(~t).*(g(~t).*((3*g(~t)+4*b0(~t)).*omega(~t)-3)...
            -3*omega(~t))-4*omega(~t)+7)./hd(~t);
        gam(~t,2) = (g(~t).*(g(~t).*((3*g(~t)+4*(b0(~t)-1)).*omega(~t)+1)...
            -3*omega(~t))+4*omega(~t)-1)./hd(~t);
    else
        gam(:,1) = (g.*(g.*((3*g+4*b0).*omega-3)-3*omega)-4*omega+7)./hd;
        gam(:,2) = (g.*(g.*((3*g+4*(b0-1)).*omega+1)-3*omega)+4*omega-1)./hd;
    end
end

gam(:,3) = b0;
gam(:,4) = 1-gam(:,3);

end

function beta0 = betanaught( mu0, g )
% beta0 = betanaught( mu0, g )
%compute Meador-Weaver beta_0 value
%   Equations 3 and 4
%
% Input
%   mu0 - cosine of illumination angle
%   g - scattering asymmetry parameter

% parameters
MAXNO = 2048;
TOL = eps('single')*10;

% check inputs
assert(all(mu0 >= 0 & mu0 <= 1))
assert(all(g >= 0 & g <= 1))

% sum until convergence; we use the even terms only for the recursive
% calculation of the Legendre polynomials
bsum = zeros(size(g));
for k=1:length(g)
    % Legendre polynomials of degree 0 and 1
    pnm2 = 1;
    pnm1 = mu0(k);
    % first coefficients and initial sum
    fm = -1/8;
    gn = 7 * g(k)^3;
    bsum(k) = 3 * g(k) * mu0(k) / 2;
    if g(k) ~= 0 && mu0(k) ~= 0
        for n=2:MAXNO
            % order n Legendre polynomial
            pn = ((2 * n - 1) * mu0(k) * pnm1 + (1 - n) * pnm2) / n;
            if mod(n,2) == 1 % odd n
                last = bsum(k);
                bsum(k) = last + gn * fm * pn;
                if abs((bsum(k) - last) / bsum(k)) < TOL && bsum(k) <= 1
                    break;
                end
                % recursively find next f(m) and gn coefficients n = 2 * m + 1
                m = (n - 1) / 2;
                fm = fm * (-(2 * m + 1) / (2 * (m + 2)));
                gn = gn * (g(k)^2 * (4 * m + 7) / (4 * m + 3));
            end
            % ready to compute next Legendre polynomial
            pnm2 = pnm1;
            pnm1 = pn;
        end
        if n >= MAXNO
            conv = (bsum(k)-last)/bsum(k);
            if conv>TOL*10
                warning('%s: %s - mu0=%g g=%g bsum=%g last=%g conv=%g fm=%g',...
                    'betanaught', 'no convergence',...
                    mu0(k), g(k), bsum(k), last, conv, fm);
            end
            if (bsum(k) > 1)
                bsum(k) = 1;
            end
        else
            assert(bsum(k) >= 0 && bsum(k) <= 1,'numerical error');
        end
    else
        bsum(k) = 0;
    end
end

beta0 = (1-bsum)/2;

end