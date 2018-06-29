function M = lookupMie(substance,radius,radiusUnits,lambda,lambdaUnits,varargin)
% M = lookupMie(substance,radius,radiusUnits,lambda,lambdaUnits [,waterFraction])
%lookupMie look up Mie variables for ice, water, dust, or soot
%
%Input
%   substance - either 'ice' or 'snow' (same), 'water', 'wetSnow',
%       'dust', or 'soot' (various different dusts to be added)
%   radius - optical radius of particle
%   radiusUnits - typically 'um', 'mm', etc
%   lambda - wavelengths
%   lambdaUnits = typically 'nm', 'um'
%   (if both radius and lambda are vectors or matrices, they must be same
%   size)
%
%Optional input needed if 'wetSnow' is specified
%   waterFraction - mass fraction of liquid water
%
%Output
%   M - structure with variables
%       Qext - extinction efficiency
%       Qabs - absorption efficiency
%       Qsca - scattering efficiency
%       Qpr - radiation pressure efficiency
%       omega - single scattering albedo
%       g - Mie asymmetry parameter

persistent Fice Fwater Fdust Fsoot mievariables

p = inputParser;
rangeValidation = @(x) isnumeric(x) && all(x(:)>=0);
positiveValidation = @(x) isnumeric(x) && all(x(:)>0);
addRequired(p,'substance',@ischar)
addRequired(p,'radius',positiveValidation)
addRequired(p,'radiusUnits',@ischar)
addRequired(p,'lambda',positiveValidation)
addRequired(p,'lambdaUnits',@ischar)
addOptional(p,'waterFraction',0,rangeValidation)
parse(p,substance,radius,radiusUnits,lambda,lambdaUnits,varargin{:})

% load lookup table
if strcmpi(p.Results.substance,'snow')
    filename = ['LUT_Mie_' 'ice' '.mat'];
else
    filename = ['LUT_Mie_' substance '.mat'];
end
switch p.Results.substance
    case {'ice','snow'}
        if isempty(Fice)
            S = load(filename);
            Fice = S.F;
            mievariables = S.variables;
        end
        F = Fice;
    case 'water'
        if isempty(Fwater)
            S = load(filename);
            Fwater = S.F;
            mievariables = S.variables;
        end
        F = Fwater;
    case 'dust'
        if isempty(Fdust)
            S = load(filename);
            Fdust = S.F;
            mievariables = S.variables;
        end
        F = Fdust;
    case 'soot'
        if isempty(Fsoot)
            S = load(filename);
            Fsoot = S.F;
            mievariables = S.variables;
        end
        F = Fsoot;
    case 'wetsnow'
        error('''wetSnow'' not yet implemented')
    otherwise
        error('''substance'' %s unrecognized',p.Results.substance)
end

[radius,lambda] = checkSizes(p.Results.radius,p.Results.lambda);

%lookup tables use 'mum' for wavelength and radius
radius = convertUnits(radius,p.Results.radiusUnits,'mum');
lambda = convertUnits(lambda,p.Results.lambdaUnits,'mum');

% lookup tables use sqrt of radius and log of wavelength
radius = sqrt(radius);
lambda = log(lambda);

% populate output structure with interpolation function
for k=1:length(mievariables)
    thisF = F{k};
    % for omega near 1.0, lookup uses 1-omega
    if strcmp(mievariables{k},'omega') &&...
            (strcmp(p.Results.substance,'ice') ||...
            strcmp(p.Results.substance,'snow') ||...
            strcmp(p.Results.substance,'water'))
        M.(mievariables{k}) = 1-thisF(radius,lambda);
    else
        M.(mievariables{k}) = thisF(radius,lambda);
    end
end

% warn if maybe out of range
for k=1:length(mievariables)
    if any(isnan(M.(mievariables{k})))
        warning('variable %s, %d NaN values, inputs maybe out of range, check SnowCloudLimits for permissible values, and check your inputs',...
            mievariables{k},nnz(isnan(M.(mievariables{k}))))
    end
end

end