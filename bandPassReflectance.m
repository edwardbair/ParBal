function [ R,varargout ] = bandPassReflectance(lambda, lambdaUnit, spectralReflectance, varargin)
% [ R ] = bandPassReflectance(lambda, lambdaUnit, spectralReflectance, 'bandPass')
% [ R ] = bandPassReflectance(lambda, lambdaUnit, spectralReflectance, 'sensor', 'bands')
% [ R ] = bandPassReflectance(_____________, 'ignoreSolar', true or false
% [R, bandPass] = bandPassReflectance(___________)
%spectral reflectance converted to integration over wavelength ranges
%
%Inputs (if not scalars, lambda and spectralReflectance must be same size)
% lambda - wavelengths of spectralReflectance values
% lambdaUnit - any length unit, typically 'um', 'mum', or 'nm'
% spectralReflectance - vector of same size as lambda
%Optional input, name-value pair (either 'bandPass' or 'sensor'/'bands'
%   must be specified, but not both
% 'bandPass' - matrix of size Nx2, specifying the wavelength ranges of N band
%   passes, must be in the same units as lambdaUnit
% 'sensor' - followed by character string, any sensor from SensorTable.m
% 'band' - either numeric vector, cell vector, or categorical vector of bands,
%   of, if omitted, all bands for that sensor
%   (numeric vector can be used for those sensors whose bands are designated just by number)
% 'ignoreSolar', if true ignores solar radiation and just provides
%   band-average reflectivity, default false, in which case solar radiation
%   accounted for unless outside range of SolarScale.m
%   (this is needed to calculate emissivities around 4 um)
%
%Output
% R - band-averaged reflectance values, dimensionless matrix of size Nx1
% bandPass (optional) - band passes of the input

% parse inputs
narginchk(4,8)
nargoutchk(0,2)
p = inputParser;
nonNegativeFcn = @(x) isnumeric(x)&&all(x(:)>=0);
bpValidation = @(x) isnumeric(x) && all(x(:)>=0) && size(x,2)==2;
bandValidation = @(x) (isrow(x) || iscolumn(x)) &&...
    ((isnumeric(x) && all(x(:)>0)) || iscell(x) || iscategorical(x));
addRequired(p,'lambda',nonNegativeFcn)
addRequired(p,'lambdaUnit',@ischar)
addRequired(p,'spectralReflectance',nonNegativeFcn)
addParameter(p,'bandpass',[],bpValidation)
addParameter(p,'sensor','',@ischar)
addParameter(p,'band',[],bandValidation)
addParameter(p,'ignoresolar',false,@islogical)
parse(p,lambda,lambdaUnit,spectralReflectance,varargin{:})

% check sizes
lambda = p.Results.lambda;
spectralReflectance = p.Results.spectralReflectance;
assert(iscolumn(lambda) || iscolumn(spectralReflectance) ||...
    isrow(lambda) || isrow(spectralReflectance),...
    'lambda and spectralReflectance must be vectors, not matrices')
assert(isequal(numel(lambda),numel(spectralReflectance)),...
    'lambda and spectralReflectance vectors must be same size')

% check consistence of band designations
bandPass = getBands(p);

% use column vectors for both wavelength and reflectance
lambda = lambda(:);
spectralReflectance = spectralReflectance(:);

passUnit = p.Results.lambdaUnit;

% if within solar spectrum, weight by solar scaled value
% if outsize solar spectrum (i.e. infrared) or if ignoreSolar is true, use average
issolar = ~p.Results.ignoresolar;

% interpolation function for spectral reflectance
Frefl = fit(lambda,spectralReflectance,'pchip');

% make sure bandPass values are within spectrum
assert(all(bandPass(:)>=min(lambda)),...
    'check input - some bandPass wavelengths < min(lambda) %g',min(lambda))
assert(all(bandPass(:)<=max(lambda)),...
    'check input - some bandPass wavelengths > max(lambda) %g',max(lambda))

R = zeros(size(bandPass,1),1);
for k=1:size(bandPass,1)
    if issolar
        denom = integral(@integrandS,bandPass(k,1),bandPass(k,2));
        if denom==0
            R(k) = integral(@avgR,bandPass(k,1),bandPass(k,2))/...
                (bandPass(k,2)-bandPass(k,1));
        else
            num = integral(@integrandR,bandPass(k,1),bandPass(k,2));
            R(k) = num/denom;
        end
    else
        R(k) = integral(@avgR,bandPass(k,1),bandPass(k,2))/...
            (bandPass(k,2)-bandPass(k,1));
    end
end

if nargout>1
    varargout{1} = bandPass;
end

    function X = integrandR(wavelength)
        refl = Frefl(wavelength)';
        S = SolarScale(wavelength,'units',passUnit,'Location','surface');
        X = (refl.*S);
    end
    function S = integrandS(wavelength)
        S = SolarScale(wavelength,'units',passUnit,'Location','surface');
    end
    function X = avgR(wavelength)
        X = Frefl(wavelength)';
    end
end

function bandPass = getBands(p)
% check consistence of band designations
assert(xor(isempty(p.Results.bandpass),isempty(p.Results.sensor)),...
    'either ''bandPass'' or ''sensor''/''bands'' must be specified')
if ~isempty(p.Results.bandpass)
    bandPass = p.Results.bandpass;
else
    T = SensorTable(p.Results.sensor,p.Results.lambdaUnit);
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