function [ R ] = SnowCloudSpectralReflectance(cosZ,substance,radius, radiusUnits, lambda, lambdaUnits, varargin)
% R = SnowCloudSpectralReflectance(cosZ,substance,radius,radiusUnits,lambda,lambdaUnits [, ...)
%spectral reflectance of snow
%
%Inputs (if not scalars, must be same size)
% cosZ - cosine of illumination angle
% substance - 'ice' or 'snow' (same), 'water', or 'wetSnow'
%   (to model a mixed phase cloud, specifiy 'ice' or 'water' as the
%   substance, and the opposite as the contaminant, but wet snow cannot be
%   mimicked in this way because water and ice form clusters rather than
%   independent scatterers)
% radius - either scalar or vector of size(lambda)
% radiusUnits - any metric length ('angstrom', 'nm', 'um' or 'mum', 'mm', 'cm', 'm')
% lambda - wavelength of light
% lambdaUnits - any length like radiusUnits, or 'MHz', 'GHz', 'THz'
%
%Optional input, name-value pairs in any order
% 'WE' - snow or cloud water equivalent (Inf if not specified)
% 'waterConc' - mass fraction of water, needed only if 'wetSnow' is the
%   substance
% 'weUnits' - any metric length (typically 'mm', 'cm', or 'm'), default 'mm'
% 'R0' - reflectance of underlying ground (ignored unless WE specified),
%   can be scalar or same size as lambda, default 0.0,
% 'contam', either a character string or cell vector, options are 'soot' or
%   'dust' for dirty snow, or 'water' for mixed-phase cloude, default is
%   'neither', in which case snow or cloud is assumed clean and not of
%   mixed phase
% 'contamRadius', followed by effective radii of contaminants, same
%   radiusUnits as for grain size and vector of same size as 'contam'
% 'contamConc', concentration of contaminants, as a mass fraction, vector
%	of same size as 'contam'
% 'lookup', if true, lookup tables are called for Mie calculations if they
%   exist, with a warning if they don't or if input values are out of range,
%   if false, Mie calculations are done (default true)
%
%Output
% R - snow or cloud reflectance, dimensionless, same size dimensions as
%   radius or lambda

narginchk(6,22)
nargoutchk(0,1)

% parse inputs (converts radii to mum, wavelengths to nm, WE to mm,
sizeUnit = 'mum';
waveUnit = 'nm';
optargin=size(varargin,2);
% assert (mod(optargin,2)==0,'must be even number of optional arguments')
iStruct = parseInput(cosZ,substance,radius,radiusUnits,lambda,lambdaUnits,varargin{:});

% Mie calculations
if iStruct.lookup % Mie parameters by lookup table
    M = lookupMie(iStruct.substance,iStruct.radius,sizeUnit,...
        iStruct.lambda,waveUnit,iStruct.waterConc);
    if ~iStruct.cleanSnowCloud
        for k=1:length(iStruct.contaminant)
            D(k) = lookupMie(iStruct.contaminant{k},iStruct.contamRadius(k),...
                sizeUnit,iStruct.lambda,waveUnit); %#ok<AGROW>
        end
    end
else % Mie scattering parameters (will do the scalar to vector conversion if
    % needed)
    if iStruct.wetSnow
        CRefIn = RefractiveIndex(iStruct.lambda,'ice',waveUnit)*(1-iStruct.waterConc) +...
            RefractiveIndex(iStruct.lambda,'water',waveUnit)*iStruct.waterConc;
        M = MieSphere(iStruct.radius,sizeUnit,iStruct.lambda,waveUnit,...
            'refindex',CRefIn);
    else
        M = MieSphere(iStruct.radius,sizeUnit,iStruct.lambda,...
            waveUnit,'substance',iStruct.substance);
    end
    if ~cleanSnowCloud
        for k=1:length(iStruct.contaminant)
            D(k) = MieSphere(iStruct.contamRadius(k),sizeUnit,iStruct.lambda,...
                waveUnit,'substance',iStruct.contaminant{k}); %#ok<AGROW>
        end
    end
end

% averaged Mie parameters if dirty snow, or mixed phase
% cloud
if ~iStruct.cleanSnowCloud || iStruct.mixedPhase
    density = zeros(1,1+length(iStruct.contaminant));
    allM = cell(size(density));
    conc = zeros(size(density));
    allRad = zeros(size(density));
    allM{1} = M;
    allRad(1) = unique(iStruct.radius);
    assert(length(unique(iStruct.radius))==1,...
        'for mixtures, radius must be scalar')
    if contains(iStruct.substance,'water','IgnoreCase',true)
        density(1) = 1000;
    else
        density(1) = 917;
    end
    for k=1:length(iStruct.contaminant)
        n = 1+k;
        switch(iStruct.contaminant{k})
            case 'ice'
                density(n) = 917;
            case 'water'
                density(n) = 1000;
            case 'dust'
                density(n) = 2800;
            case 'soot'
                density(n) = 2800;
            otherwise
                error('''contam'' %s not recognized',iStruct.contaminant{k})
        end
        allM{n} = D(k);
        conc(n) = iStruct.contamConc(k);
        allRad(n) = iStruct.contamRadius(k);
    end
    conc(1) = 1-sum(conc(2:end));
    mix = MieMixture(allM,allRad,density,conc);
    M = mix;
end
% vector/matrix sizes of mu0 and omega must be equal
[mu0,omega] = checkSizes(iStruct.mu0,M.omega);

if all(isinf(iStruct.WE))
    tau = Inf(size(omega));
else
    switch lower(iStruct.substance)
        case 'ice'
            if iStruct.mixedPhase
                k = find(contains(iStruct.contaminant,'water'));
                mfrac = [1-iStruct.waterConc iStruct.waterConc];
                vfrac = mfrac./[density(1) density(k)];
                vfrac = vfrac/sum(vfrac);
                [ri,rw,WE,Qext,wQext] =...
                    checkSizes(convertUnits(iStruct.radius,sizeUnit,'m'),...
                    convertUnits(iStruct.contamRadius(k),sizeUnit,'m'),...
                    iStruct.WE,M.Qext,D(k).Qext);
                tau = tauSnow(ri,WE,Qext)*vfrac(1)+tauCloud(rw,WE,wQext)*vfrac(2);
            else
                [r,WE,Qext] = checkSizes(convertUnits(iStruct.radius,sizeUnit,'m'),iStruct.WE,M.Qext);
                tau = tauSnow(r,WE,Qext);
            end
            
        case 'water'
            if iStruct.mixedPhase
                k = find(contains(iStruct.contaminant,'ice'));
                mfrac = [iStruct.waterConc 1-iStruct.waterConc];
                vfrac = mfrac./[density(1) density(k)];
                vfrac = vfrac/sum(vfrac);
                [rw,ri,WE,Qext,iQext] =...
                    checkSizes(convertUnits(iStruct.radius,sizeUnit,'m'),...
                    convertUnits(iStruct.contamRadius(k),sizeUnit,'m'),...
                    iStruct.WE,M.Qext,D(k).Qext);
                tau = tauSnow(ri,WE,iQext)*vfrac(2)+tauCloud(rw,WE,Qext)*vfrac(1);
            else
                [r,WE,Qext] =...
                    checkSizes(convertUnits(iStruct.radius,sizeUnit,'m'),...
                    iStruct.WE,M.Qext);
                tau = tauCloud(r,WE,Qext);
            end
        otherwise
            error('''substance'' ''%s'' not recognized',p.Results.substance)
    end
end

[g,R0,~] = checkSizes(M.g,iStruct.R0,omega);
R = twostream(mu0,omega,g,tau,'R0',R0);
if iStruct.changeSize
    R = reshape(R,iStruct.origSize);
end

end

function iStruct = parseInput(cosZ,substance,radius,radiusUnits,lambda,lambdaUnits,varargin)
%parse input values to produce input/output raster references

% some parameters
sizeUnit = 'mum';
waveUnit = 'nm';
S = SnowCloudLimits;
defaultWE = Inf;
defaultWEunits = 'mm';
defaultR0 = 0;
defaultContam = 'neither';

p = inputParser;
rangeValidation = @(x) isnumeric(x) && all(x(:)>=0 & x(:)<=1);
positiveValidation = @(x) isnumeric(x) && all(x(:)>0);
contamValidation = @(x) all(isnan(x)) || (isnumeric(x) && all(x(:)>=0 & x(:)<1));
contamID = @(x) ischar(x) || iscell(x);
addRequired(p,'cosZ',rangeValidation)
addRequired(p,'substance',@ischar);
addRequired(p,'radius',@isnumeric)
addRequired(p,'radiusUnits',@ischar)
addRequired(p,'lambda',positiveValidation)
addRequired(p,'lambdaUnits',@ischar)
addParameter(p,'we',defaultWE,positiveValidation)
addParameter(p,'weunits',defaultWEunits,@ischar)
addParameter(p,'r0',defaultR0,rangeValidation)
addParameter(p,'contam',defaultContam,contamID)
addParameter(p,'contamradius',0,positiveValidation)
addParameter(p,'contamconc',0,contamValidation)
addParameter(p,'lookup',true,@islogical)
addParameter(p,'waterconc',0,rangeValidation)
parse(p,cosZ,substance,radius,radiusUnits,lambda,lambdaUnits,varargin{:})

% some logical variables
iStruct.wetSnow = contains(p.Results.substance,'wetsnow','IgnoreCase',true);
iStruct.mixedPhase = (contains(p.Results.substance,'ice','IgnoreCase',true) &&...
    any(contains(p.Results.contam,'water','IgnoreCase',true))) ||...
    (contains(p.Results.substance,'water','IgnoreCase',true) &&...
    any(contains(p.Results.contam,'ice','IgnoreCase',true)));

% begin transfer to data structure
substance = lower(p.Results.substance);
assert(strcmp(substance,'ice') || strcmp(substance,'snow') ||...
    strcmp(substance,'water') || strcmp(substance,'wetsnow'),...
    'substance ''%s'' not recognized',p.Results.substance)
iStruct.substance = substance;

% check for some incompatibilities
if ischar(p.Results.contam)
    thisContam = {p.Results.contam};
else
    thisContam = p.Results.contam;
end
[contaminant,contamRadius,contamConc] =...
    checkSizes(thisContam,p.Results.contamradius,p.Results.contamconc);
assert(length(contaminant)==length(unique(contaminant)),...
    '''contam'' cell vector must not contain duplicates')
if iStruct.wetSnow
    assert(~isempty(p.Results.waterconc) && p.Results.waterconc<=0.15,...
        'if ''substance'' is ''wetSnow'', ''waterConc'' must be <= 0.15')
    iStruct.waterConc = p.Results.waterconc;
elseif iStruct.mixedPhase
    if contains(substance,'ice')
        k = contains(contaminant,'water');
        iStruct.waterConc = contamConc(k);
    elseif contains(substance,'water')
        iStruct.waterConc = 1-sum(contamConc);
    end
else
    iStruct.waterConc = 0;
end

% check and adjust sizes
[mu0,radius,R0,lambda] = checkSizes(p.Results.cosZ,p.Results.radius,...
    p.Results.r0,p.Results.lambda);
if iscolumn(lambda)
    iStruct.changeSize = false;
else
    iStruct.changeSize = true;
    iStruct.origSize = size(lambda);
end
iStruct.mu0 = mu0(:);
iStruct.radius = convertUnits(radius(:),p.Results.radiusUnits,sizeUnit);
iStruct.R0 = R0(:);
iStruct.lambda = convertUnits(lambda(:),p.Results.lambdaUnits,waveUnit);

if isinf(p.Results.we)
    iStruct.WE = p.Results.we;
else
    iStruct.WE = convertUnits(p.Results.we,p.Results.weunits,'mm');
end

if any(contains(contaminant,defaultContam))
    cleanSnowCloud = true;
else
    cleanSnowCloud = false;
    for k=1:length(contaminant)
        switch contaminant{k}
            case 'dust'
                if contamRadius(k)==0
                    contamRadius(k) = convertUnits(S.defaultDustRadius,S.unitsSize,...
                        p.Results.radiusUnits);
                end
            case 'soot'
                if contamRadius(k)==0
                    contamRadius(k) = convertUnits(S.defaultSootRadius,S.unitsSize,...
                        p.Results.radiusUnits);
                end
            case {'water','ice'}
                if contamRadius(k)==0
                    contamRadius(k) = p.Results.radius(1);
                end
            otherwise
                error('''contam'' ''%s'' not recognized',contaminant{k})
        end
    end
    iStruct.contamRadius = convertUnits(contamRadius,p.Results.radiusUnits,sizeUnit);
end

iStruct.cleanSnowCloud = cleanSnowCloud;
if ~cleanSnowCloud
    iStruct.contaminant = contaminant;
    iStruct.contamConc = contamConc;
    iStruct.contamRadius = contamRadius;
end

% use Mie Lookup tables
iStruct.lookup = p.Results.lookup;
end