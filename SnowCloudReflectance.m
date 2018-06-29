function [ R ] = SnowCloudReflectance(cosZ,substance,radius, radiusUnit, varargin)
% [ R ] = SnowCloudReflectance(cosZ,substance,radius, radiusUnit, varargin)
%reflectance of snow or cloud across band passes or sensor bands
%
%Inputs (if not scalars, must be of same size)
% cosZ - cosine of illumination angle
% substance - must be scalar, 'ice' or 'snow' (same), 'water', or 'wetSnow'
%   (to model a mixed phase cloud, specifiy 'ice' or 'water' as the
%   substance, and the opposite as the contaminant, but wet snow cannot be
%   mimicked in this way because water and ice form clusters rather than
%   independent scatterers)
% radius - either scalar or vector
% radiusUnit - any metric length ('angstrom', 'nm', 'um' or 'mum', 'mm', 'cm', 'm')
%
%Optional input, name-value pairs in any order
%Arguments about the bands (either 'bandPass' or 'sensor'/'bands' must
%be specified, but not both
% 'bandPass' - matrix of size Nx2, specifying the wavelength ranges of N band
%   passes
% 'waveUnit' - wavelength units of bandPass (not needed when 'sensor' is
%   specified)
% 'sensor' - followed by character string, any sensor from SensorTable.m
% 'band' - either numeric vector, cell vector, or categorical vector of bands,
%   of, if omitted, all bands for that sensor
%   (numeric vector can be used for those sensors whose bands are designated just by number)
% 'ignoreSolar', if true ignores solar radiation and just provides
%   band-average reflectivity, default false, in which case solar radiation
%   accounted for unless outside range of SolarScale.m
%   (this is needed to calculate emissivities around 4 um)
%Arguments about the snow or cloud
% 'WE' - snow or cloud water equivalent (Inf if not specified)
% 'weUnits' - any metric length (typically 'mm', 'cm', or 'm'), default 'mm'
% 'waterConc' - mass fraction of water, needed only if 'wetSnow' is the
%   substance
% 'R0' - reflectance of underlying ground (ignored unless WE specified),
%   can be scalar or same size as number of bands, default 0.0,
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
%   radius or cosZ

% hang onto wavelength values corresponding to refractive indices
persistent wvice

narginchk(5,25)
nargoutchk(0,1)

% parse inputs (convert radii to mum, wavelengths to nm, WE to mm,
sizeUnit = 'mum';
waveUnit = 'nm';
optargin=size(varargin,2);
assert (mod(optargin,2)==0,'must be even number of optional arguments')
iStruct = parseInput(cosZ,substance,radius,radiusUnit,varargin{:});

% wavelengths to consider is the full range of all band passes at
% resolution of refractive index data
if isempty(wvice)
    [~,wvice] = RefractiveIndex([],'ice',waveUnit);
end

% option single band, many radii
if size(iStruct.bandPass,1)==1
    k1 = find(wvice<min(iStruct.bandPass(1,:)),1,'last');
    k2 = find(wvice>max(iStruct.bandPass(1,:)),1,'first');
    thisWave = wvice(k1:k2);
    R = zeros(size(iStruct.radius));
    parfor k=1:numel(iStruct.radius)
        if iStruct.cleanSnowCloud %#ok<PFBNS>
            thisR = SnowCloudSpectralReflectance(iStruct.mu0(k),iStruct.substance,...
                iStruct.radius(k),sizeUnit,...
                thisWave,waveUnit,...
                'WE',iStruct.WE,...
                'R0',iStruct.R0(k),...
                'waterConc',iStruct.waterConc);
            if any(thisR(:)<0)
                thisR(thisR<0) = 0;
            end
        else
            thisR = SnowCloudSpectralReflectance(iStruct.mu0(k),iStruct.substance,...
                iStruct.radius(k),sizeUnit,...
                thisWave,waveUnit,...
                'WE',iStruct.we,...
                'R0',iStruct.R0(k),...
                'contam',iStruct.contaminant,...
                'contamConc',iStruct.contamConc,...
                'contamRadius',iStruct.contamRadius,...
                'waterConc',iStruct.waterConc);
        end
        R(k) = bandPassReflectance(thisWave,waveUnit,thisR,...
            'bandPass',iStruct.bandPass,'IgnoreSolar',iStruct.ignoreSolar);
    end
else
    % snow reflectance to cover all bands, number of radii <= number of bands
    snowR = cell(size(iStruct.bandPass,1),1);
    wave = cell(size(snowR));
    for k=1:size(iStruct.bandPass,1)
        k1 = find(wvice<min(iStruct.bandPass(k,:)),1,'last');
        k2 = find(wvice>max(iStruct.bandPass(k,:)),1,'first');
        wave{k} = wvice(k1:k2);
        if iStruct.cleanSnowCloud
            snowR{k} = SnowCloudSpectralReflectance(iStruct.mu0(k),iStruct.substance,...
                iStruct.radius(k),sizeUnit,...
                wave{k},waveUnit,...
                'WE',iStruct.WE,...
                'R0',iStruct.R0(k),...
                'waterConc',iStruct.waterConc);
        else
            snowR{k} = SnowCloudSpectralReflectance(iStruct.mu0(k),iStruct.substance,...
                iStruct.radius(k),sizeUnit,...
                wave{k},waveUnit,...
                'WE',iStruct.WE,...
                'R0',iStruct.R0(k),...
                'contam',iStruct.contaminant,...
                'contamConc',iStruct.contamConc,...
                'contamRadius',iStruct.contamRadius,...
                'waterConc',iStruct.waterConc);
        end
    end
    
    R = zeros(size(iStruct.bandPass,1),1);
    
    for k=1:size(R,1)
        R(k) = bandPassReflectance(wave{k},waveUnit,snowR{k},...
            'bandPass',iStruct.bandPass(k,:),'IgnoreSolar',iStruct.ignoreSolar);
    end
end
end
function iStruct = parseInput(cosZ,substance,radius,radiusUnit,varargin)
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
bandValidation = @(x) (isrow(x) || iscolumn(x)) &&...
    ((isnumeric(x) && all(x(:)>0)) || iscell(x) || iscategorical(x));
bpValidation = @(x) isnumeric(x) && all(x(:)>=0) && size(x,2)==2;
contamID = @(x) ischar(x) || iscell(x);
addRequired(p,'cosZ',rangeValidation)
addRequired(p,'substance',@ischar);
addRequired(p,'radius',@isnumeric)
addRequired(p,'radiusUnit',@ischar)
addParameter(p,'we',defaultWE,positiveValidation)
addParameter(p,'weunits',defaultWEunits,@ischar)
addParameter(p,'waveunit',waveUnit,@ischar)
addParameter(p,'r0',defaultR0,rangeValidation)
addParameter(p,'contam',defaultContam,contamID)
addParameter(p,'contamradius',0,positiveValidation)
addParameter(p,'contamconc',0,contamValidation)
addParameter(p,'lookup',true,@islogical)
addParameter(p,'waterconc',0,rangeValidation)
addParameter(p,'bandpass',[],bpValidation)
addParameter(p,'sensor','',@ischar)
addParameter(p,'band',[],bandValidation)
addParameter(p,'ignoresolar',false,@islogical)
parse(p,cosZ,substance,radius,radiusUnit,varargin{:})

% some logical variables
iStruct.wetSnow = contains(p.Results.substance,'wetsnow','IgnoreCase',true);
iStruct.mixedPhase = (contains(p.Results.substance,'ice','IgnoreCase',true) &&...
    any(contains(p.Results.contam,'water','IgnoreCase',true))) ||...
    (contains(p.Results.substance,'water','IgnoreCase',true) &&...
    any(contains(p.Results.contam,'ice','IgnoreCase',true)));

% begin transfer to data structure
% sensor/band characteristics
iStruct.bandPass = getBands(p);
iStruct.ignoreSolar = p.Results.ignoresolar;

% snow/cloud characteristics
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
if ~isscalar(p.Results.radius) || ~isscalar(p.Results.cosZ)
    assert(size(iStruct.bandPass,1)==1,...
        'if radius or cosZ not scalar, bandPass can have just one band')
    [mu0,radius,R0,WE]= checkSizes(p.Results.cosZ,p.Results.radius,...
        p.Results.r0,p.Results.we);
else
    [mu0,radius,R0,WE,~] = checkSizes(p.Results.cosZ,p.Results.radius,...
        p.Results.r0,p.Results.we,iStruct.bandPass(:,1));
end
iStruct.mu0 = mu0;
iStruct.radius = convertUnits(radius,p.Results.radiusUnit,sizeUnit);
iStruct.R0 = R0;
iStruct.WE = WE;

if ~isinf(p.Results.we)
    iStruct.WE = convertUnits(iStruct.WE,p.Results.weunit,'mm');
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
                        p.Results.radiusUnit);
                end
            case 'soot'
                if contamRadius(k)==0
                    contamRadius(k) = convertUnits(S.defaultSootRadius,S.unitsSize,...
                        p.Results.radiusUnit);
                end
            case {'water','ice'}
                if contamRadius(k)==0
                    contamRadius(k) = p.Results.radius(1);
                end
            otherwise
                error('''contam'' ''%s'' not recognized',contaminant{k})
        end
    end
    iStruct.contamRadius = convertUnits(contamRadius,p.Results.radiusUnit,sizeUnit);
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

function bandPass = getBands(p)
% check consistence of band designations
assert(xor(isempty(p.Results.bandpass),isempty(p.Results.sensor)),...
    'either ''bandPass'' or ''sensor''/''bands'' must be specified')
waveUnit = 'nm';
if ~isempty(p.Results.bandpass)
    bandPass = convertUnits(p.Results.bandpass,p.Results.waveunit,waveUnit);
else
    T = SensorTable(p.Results.sensor,waveUnit);
    if isempty(p.Results.band)
        % all bands
        bandPass = [T.LowerWavelength T.UpperWavelength];
    else
        x = p.Results.band;
        if isnumeric(x)
            band = categorical(x);
        elseif iscategorical(x)
            band = x;
        else % cell
            band = categorical(x);
        end
        bandPass = zeros(length(x),2);
        for k=1:length(band)
            b = find(T.Band==band(k));
            if isempty(b)
                warning(['band ' band(k) ' not found'])
            end
            bandPass(k,:) = [T.LowerWavelength(b) T.UpperWavelength(b)];
        end
    end
end
% make sure bigger number is in column 2
bandPass = sort(bandPass,2);
end