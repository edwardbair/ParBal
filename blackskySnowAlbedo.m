function [ albedo ] = blackskySnowAlbedo(cosZ,grainSize,radiusUnits,varargin )
% [ albedo ] = blackskySnowAlbedo(cosZ,grainSize,radiusUnits,varargin )
%calculates broadband albedo of snow for clear-sky conditions
%(cloudy option to be added)
%Method: calculate spectral reflectance from Mie theory and twostream
%radiative transfer
%
%Input
% cosZ - cosine of illumination angle, scalar or vector or matrix
% grainSize - optically equivalent radius, scalar or vector or matrix
% radiusUnits - typically 'um', 'mum', or 'mm'
% (if cosZ and grainSize are both non-scalar, they must be same size)
%
%Optional input, name-value pair, if vector or matrix, must be same size as
%   cosZ or grainSize
% 'wavelengthRange' - of measurement, vector of length 2, default whole
%       solar spectrum 0.28 to 5.0 um, but can be specified to match other
%       radiometers
% 'waveUnits' - units for wavelength, default um if not specified
% 'contam', followed by 'soot' or 'dust' (program deals with either, but
%       not both in combination, default is 'neither', in which case snow is
%       assumed clean)
% 'contamConc', concentration of particulates, as a mass fraction
% 'particSize', followed by effective radius of particulates, same 'units'
%       as for grain size (default 1 um for dust, 10 nm for soot)
% 'WE' - snow water equivalent (Inf if not specified), scalar or
%   vector
% 'weUnits' - any metric length (typically 'mm', 'cm', or 'm')
% 'R0' - reflectance of underlying ground (ignored unless WE specified),
%   scalar but choose value appropriate for visible wavelengths
%
%Output
% albedo - same size as input cosZ or grainSize

% parse inputs
defaultWE = Inf;
defaultWEunits = 'mm';
defaultR0 = 0;
defaultContam = 'neither';
defaultWavelengthRange = [.28 5];
defaultWavelengthUnits = 'um';
p = inputParser;
rangeValidation = @(x) isnumeric(x) && all(x(:)>=0 & x(:)<=1);
positiveValidation = @(x) isnumeric(x) && all(x(:)>=0);
contamValidation = @(x) all(isnan(x)) || (isnumeric(x) && all(x(:)>=0 & x(:)<1));
addRequired(p,'cosZ',rangeValidation)
addRequired(p,'grainSize',@isnumeric)
addRequired(p,'radiusUnits',@ischar)
addParameter(p,'wavelengthrange',defaultWavelengthRange,positiveValidation)
addParameter(p,'waveunits',defaultWavelengthUnits,@ischar)
addParameter(p,'we',defaultWE,positiveValidation)
addParameter(p,'weunits',defaultWEunits,@ischar)
addParameter(p,'r0',defaultR0,rangeValidation)
addParameter(p,'contam',defaultContam,@ischar)
addParameter(p,'particsize',0,positiveValidation)
addParameter(p,'contamconc',[],contamValidation)
parse(p,cosZ,grainSize,radiusUnits,varargin{:})
cleanSnow = strcmpi(p.Results.contam,defaultContam);
deepSnow = isinf(p.Results.we);
if cleanSnow
    [mu0,radius] = checkSizes(p.Results.cosZ,p.Results.grainSize);
else
    [mu0,radius,contamConc] = checkSizes(p.Results.cosZ,...
        p.Results.grainSize,p.Results.contamconc);
end
if ~deepSnow
    [~,WE,R0] = checkSizes(radius,p.Results.we,p.Results.r0);
end

albedo = zeros(size(radius));

% trapezoidal integration (do in nm)
range = convertUnits(sort(p.Results.wavelengthrange),p.Results.waveunits,'nm');
wv = floor(range(1)):5:ceil(range(2));
wv = linspace(range(1),range(2),length(wv));
wvUnits = 'nm';
S = SolarScale(wv,'units',wvUnits,'location','surface');
if isscalar(radius)
    if cleanSnow && ~deepSnow
        R = SnowCloudSpectralReflectance(mu0,'snow',radius,radiusUnits,...
            wv,wvUnits,'WE',WE,'WEunits',p.Results.weunits,'R0',R0);
    elseif cleanSnow && deepSnow
        R = SnowCloudSpectralReflectance(mu0,'snow',radius,radiusUnits,...
            wv,wvUnits);
    elseif ~cleanSnow && deepSnow
        R = SnowCloudSpectralReflectance(mu0,'snow',radius,radiusUnits,...
            wv,wvUnits,'contam',p.Results.contam,'contamConc',contamConc,...
            'contamSize',p.Results.particsize);
    else %(~cleanSnow && ~deepSnow)
        R = SnowCloudSpectralReflectance(mu0,'snow',radius,radiusUnits,...
            wv,wvUnits,'contam',p.Results.contam,'contamConc',contamConc,...
            'contamSize',p.Results.particsize,'WE',WE,...
            'WEunits',p.Results.weunits);
    end
    albedo = trapz(wv,R.*S)/trapz(wv,S);
else
    if deepSnow
        WE = inf(size(radius));
        R0 = zeros(size(radius));
    end
    if cleanSnow
        contamConc = zeros(size(radius));
    end
    parfor k=1:numel(radius) % works okay if radius is a matrix
        if ~isnan(radius(k))
            if cleanSnow && ~deepSnow
                warning('no lookup table for shallow snow, doing Mie calculations')
                R = SnowCloudSpectralReflectance(mu0(k),'snow',radius(k),radiusUnits,...
                    wv,wvUnits,'WE',WE(k),'WEunits',p.Results.weunits,'R0',R0(k)); %#ok<PFBNS>
            elseif cleanSnow && deepSnow
                R = SnowCloudSpectralReflectance(mu0(k),'snow',radius(k),radiusUnits,...
                    wv,wvUnits,'lookup',true);
            elseif ~cleanSnow && deepSnow
                if p.Results.particsize ~= 0
                    warning('lookup tables use the default particulate size, so doing Mie calculations')
                    lookup = false;
                else
                    lookup = true;
                end
                R = SnowCloudSpectralReflectance(mu0(k),'snow',radius(k),radiusUnits,...
                    wv,wvUnits,'contam',p.Results.contam,'contamConc',contamConc(k),...
                    'contamSize',p.Results.particsize,'lookup',lookup);
            else %(~cleanSnow && ~deepSnow)
                warning('no lookup table for shallow snow, doing Mie calculations')
                R = SnowCloudSpectralReflectance(mu0(k),'snow',radius(k),radiusUnits,...
                    wv,wvUnits,'contam',p.Results.contam,'contamConc',contamConc(k),...
                    'contamSize',p.Results.particsize,'WE',WE(k),...
                    'WEunits',p.Results.weunits);
            end
            albedo(k) = trapz(wv,R.*S)/trapz(wv,S);
        end
    end
end

end