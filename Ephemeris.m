function [ declin, radiusvector, solar_lon ] = Ephemeris( matdates )
%Ephemeris: [ declin, radiusvector, solar_lon ] = Ephemeris( matdates )
% declination, radius vector, and solar longitude
%   fit to calculations from JPL Horizons
%   (http://ssd.jpl.nasa.gov/?horizons)
%
% for dates spanning 1/1/2000 thru 12/31/2018
% The Fourier coefficients and other values
% (A0D,A0R,AL,AR,BL,wD,wR,A0L,AD,ALLYEARS,BD,BR,wL)
% are loaded as globals on first pass from the file Ephem_coefficients.mat

%input
%   matdates - scalar or vector (or matrix) of MATLAB datenums (UTC)
%
%output
%   declination & solar longitude in degrees
%   radius vector in AU

persistent A0D A0R AL AR BL wD wR A0L AD ALLYEARS BD BR wL
persistent already
if isempty(already)
    load('Ephem_coefficients.mat')
    already = true;
end

declin = ifourfit( matdates, ALLYEARS, A0D, wD, AD, BD );
radiusvector = ifourfit( matdates, ALLYEARS, A0R, wR, AR, BR );
lon_offset = ifourfit( matdates, ALLYEARS, A0L, wL, AL, BL );

% solar longitude depending on time of day
[~,~,~,hour,minute,second] = datevec(matdates);
hr = hour+minute/60+second/3600;
solar_lon = lon_offset + 15*(12-hr);

function values = ifourfit( matdates, allyears, A0, w, A, B )
%ifourfit: return values associated with a Fourier fit
%
% matdates - dates at which we need function evaluated
% allyears - years in the coefficients
% A0, w, A, B - coefficients from MakeFourierFit

[yr,~,~,~,~,~] = datevec(matdates);
yu = unique(yr); % if all in same year, then length(yu)=1

values = zeros(size(matdates));
for n=1:length(yu)
    k = find(yu(n) == allyears);
    if isempty(k)
        % use 2008 if leap year, otherwise 2007
        if eomday(yu(n),2)==29
            y = 2008;
        else
            y = 2007;
        end
        warning('coefficients developed for years %d to %d, so using %d instead of %d',...
            min(allyears), max(allyears), y, yu(n))
        m = find(y == allyears);
        t = yr==yu(n);
        p = find(t,1,'first');
        difference = addtodate(matdates(p),y-yu(n),'year')-matdates(p);
        matdates(t) = matdates(t)+difference;
    else
        t = yr==yu(n);
        m = k;
    end
    wx = w(m) * matdates(t);
    swx = sin(wx);
    cwx = cos(wx);
    values(t)=A0(m)+cwx.*(A(1,m)+cwx.*(A(2,m)+cwx.*(A(3,m)+cwx.*(A(4,m)+...
        cwx.*(A(5,m)+cwx.*(A(6,m)+cwx.*(A(7,m)+A(8,m).*cwx)))))))+...
        swx.*(B(1,m)+cwx.*(2.*B(2,m)+cwx.*(3.*B(3,m)+cwx.*(4.*B(4,m)+...
        cwx.*(5.*B(5,m)+cwx.*(6.*B(6,m)+cwx.*(7.*B(7,m)+...
        8.*B(8,m).*cwx))))))+...
        swx.*(-A(2,m)+cwx.*(-3.*A(3,m)+cwx.*(-6.*A(4,m)+...
        cwx.*(-10.*A(5,m)+cwx.*(-15.*A(6,m)+cwx.*(-21.*A(7,m)-...
        28.*A(8,m).*cwx)))))+swx.*(-B(3,m)+cwx.*(-4.*B(4,m)+...
        cwx.*(-10.*B(5,m)+cwx.*(-20.*B(6,m)+cwx.*(-35.*B(7,m)-...
        56.*B(8,m).*cwx))))+swx.*(A(4,m)+cwx.*(5.*A(5,m)+...
        cwx.*(15.*A(6,m)+cwx.*(35.*A(7,m)+70.*A(8,m).*cwx)))+...
        swx.*(B(5,m)+cwx.*(6.*B(6,m)+cwx.*(21.*B(7,m)+56.*B(8,m).*cwx))+...
        swx.*(-A(6,m)+cwx.*(-7.*A(7,m)-28.*A(8,m).*cwx)+swx.*(-B(7,m)-...
        8.*B(8,m).*cwx+A(8,m).*swx)))))));
end

end
end