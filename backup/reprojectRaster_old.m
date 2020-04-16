function [ B, RB, RRB] = reprojectRaster_old(A,RA,InStruct,OutStruct,varargin )
% [ B, RB, RRB] = reprojectRaster(A,RA,InStruct,OutStruct,... )
% reproject raster from one map projection (or lat-lon) to another
% projection (or lat-lon)
%
%INPUT
%   A - input raster (can be 3D), any numeric type, or logical
%   RA - referencing matrix for A, must be empty if input data are
%       geolocated, in which case lat-lon are specified below
%   InStruct - input projection structure, [] if geographic or geolocated
%   OutStruct - output projection structure, [] if geographic
% OPTIONAL INPUT
%   name-value pairs in any order specifying . . .
%       'planet' - planet name as a character string, case insensitive,
%           defaults to 'earth'
%       'wholeHemisphere' - specify if input image is entire hemisphere or
%           whole planet, choices are 'north', 'south', 'both'
%       'method' - 'linear' (default), 'nearest' (fastest), 'cubic', or
%           'spline'
%       'rasterref' - output map or geo raster reference object
%       'lat' and 'lon' - matrices of latitude and longitude, same size as
%           first 2 dimensions of A, if input data are geolocated
%
%       The following arguments are ignored if 'rasterref' is used
%       'pixelsize' - height and width of output pixels (if scalar, both
%           same, default is to match size of input image)
%       'XLimit' and 'YLimit' - each vectors of length 2: min & max of output
%           x- and y-coordinates (default is to cover extent of A)
%       'Origin' - 'ul' (default unless 'rasterref' specified), 'll', 'ur', or 'lr'
%       'adjust' - true or false to adjust x- and y-limits to be a multiple
%           of the pixelsize (default true unless 'rasterref' specified, in
%           which case default is false)
%
%OUTPUT
%   B output reprojected raster, same class as input A
%   RB referencing matrix for B
%   RRB raster reference object for B

assert(ismatrix(A) || ndims(A)==3,...
    'input array must have 2 or 3 dimensions')

% parse inputs
optargin=size(varargin,2);
assert (mod(optargin,2)==0,'must be even number of optional arguments')
[InRasterRef,OutRasterRef,planet,method,inLat,inLon] =...
    parseInput(size(A),RA,InStruct,OutStruct,varargin{:});

% coarsen input image first if output is at significantly coarser
% resolution, so that the output is averaged over multiple input pixels
if ~isempty(InRasterRef) && strcmpi(InRasterRef.RasterInterpretation,'cells')
    [InRasterRef,A] = coarsenInput(A,InRasterRef,OutRasterRef,planet,method);
end

% world coordinates in output image
[XIntrinsic,YIntrinsic] =...
    meshgrid(1:OutRasterRef.RasterSize(2),1:OutRasterRef.RasterSize(1));
if contains(class(OutRasterRef),'MapCellsReference')
    [XWorld, YWorld] = intrinsicToWorld(OutRasterRef,XIntrinsic,YIntrinsic);
    [lat,lon] = minvtran(OutStruct,XWorld,YWorld);
elseif contains(class(OutRasterRef),'GeographicCellsReference')
    [lat,lon] = intrinsicToGeographic(OutRasterRef,XIntrinsic,YIntrinsic);
else
    error('OutRasterRef class %s unrecognized',class(OutRasterRef))
end

% input coordinates that correspond to all output coordinates
if isempty(InStruct) % input grid is lat-lon
    Xq = lon;
    Yq = lat;
else
    [Xq,Yq] = mfwdtran(InStruct,lat,lon);
end

geolocated = ~isempty(inLat); % otherwise geographic or projected
if geolocated
    X = inLon;
    Y = inLat;
else
    % all input coordinates
    [XIntrinsic,YIntrinsic] =...
        meshgrid(1:InRasterRef.RasterSize(2),1:InRasterRef.RasterSize(1));
    if contains(class(InRasterRef),'MapCellsReference')
        [X,Y] = intrinsicToWorld(InRasterRef,XIntrinsic,YIntrinsic);
    elseif contains(class(InRasterRef),'GeographicCellsReference')
        [Y,X] = intrinsicToGeographic(InRasterRef,XIntrinsic,YIntrinsic);
    else
        error('InRasterRef class %s unrecognized',class(InRasterRef))
    end
end

% interpolate to output points specified in terms of input coordinates
% input values must be either double or single, so convert if necessary
origtype = class(A);
A = double(A);

if ismatrix(A)
    if geolocated
        F = scatteredInterpolant(X(:),Y(:),A(:),method,'none');
        B = F(Xq,Yq);
    else
        B = interp2(X,Y,A,Xq,Yq,method);
    end
    %     B = F(Xq,Yq);
elseif ndims(A)==3
    for k=1:size(A,3)
        if geolocated
            V = A(:,:,k);
            F = scatteredInterpolant(X(:),Y(:),V(:),method,'none');
            B1 = F(Xq,Yq);
        else
            B1 = interp2(X,Y,A(:,:,k),Xq,Yq,method);
        end
        % allocate on first pass
        if k==1
            B = zeros(size(B1,1),size(B1,2),size(A,3),'like',A);
        end
        B(:,:,k) = B1;
    end
else
    error('arrays of more than 3 dimensions not supported')
end
% keep NaNs from propagating
t = isnan(B);
if any(t(:)) && ~strcmp(method,'nearest')
    if ismatrix(A)
        if geolocated
            F = scatteredInterpolant(X(:),Y(:),A(:),'nearest','none');
            NN = F(Xq,Yq);
        else
            NN = interp2(X,Y,A,Xq,Yq,'nearest');
        end
        B(t) = NN(t);
    else
        for k=1:size(A,3)
            B1 = B(:,:,k);
            t = isnan(B1);
            if any(t(:))
                if geolocated
                    V = A(:,:,k);
                    F = scatteredInterpolant(X(:),Y(:),V(:),'nearest','none');
                    NN = F(Xq,Yq);
                else
                    NN = interp2(X,Y,A(:,:,k),Xq,Yq,'nearest');
                end
                B1(t) = NN(t);
                B(:,:,k) = B1;
            end
        end
    end
end

% reset to original type if not double or single
if ~(strcmpi(origtype,'double') || strcmpi(origtype,'single'))
    B(isnan(B)) = 0;
    B = round(B);
end
B = cast(B,origtype);

% set the output referencing matrix
RB = RasterRefToRefMat(OutRasterRef);
RRB = OutRasterRef;
end

function [xnew,ynew] = adjustLimits(xlimit,ylimit,pixelsize)
% pixelsize is [height width]
xnew = xlimit;
ynew = ylimit;

if mod(xlimit(1),pixelsize(2))~=0
    xnew(1) = pixelsize(2)*floor(xlimit(1)/pixelsize(2));
end
if mod(xlimit(2),pixelsize(2))~=0
    xnew(2) = pixelsize(2)*ceil(xlimit(2)/pixelsize(2));
end
if mod(ylimit(1),pixelsize(1))~=0
    ynew(1) = pixelsize(1)*floor(ylimit(1)/pixelsize(1));
end
if mod(ylimit(2),pixelsize(1))~=0
    ynew(2) = pixelsize(1)*ceil(ylimit(2)/pixelsize(1));
end
% preserve aspect ratio
aspectRatio = (max(xlimit)-min(xlimit))/(max(ylimit)-min(ylimit));
newAspectRatio = (max(xnew)-min(xnew))/(max(ynew)-min(ynew));
thresh = 1.e-6;
if abs(aspectRatio-newAspectRatio)>thresh
    if newAspectRatio>aspectRatio
        ynew(2) = ynew(2)+pixelsize(1);
    else
        xnew(2) = xnew(2)+pixelsize(2);
    end
end
end

function [ newInR, newRaster ] = coarsenInput(raster,InR,OutR,planet,method)
%coarsen input raster if resolution of output is significantly coarser
%  inputs are the raster itself and the input and output raster reference
%  objects
%  output is the new input raster reference and the new (coarser) raster

thresholdRatio = 1.2;
[x11,y11,dx,dy,inProj] = cornerCoords(InR);
[~,~,odx,ody,outProj] = cornerCoords(OutR);
hdy = dy;
hdx = dx;
hody = ody;
hodx = odx;

if xor(inProj,outProj) % one projected, one geographic
    S = referenceSphere(planet);
    R = S.Radius;
    % convert the geographic one to projected
    if inProj
        hody = degtorad(ody)*R;
        hodx = degtorad(odx)*R*(cosd(mean(OutR.LatitudeLimits)));
    else
        hdy = degtorad(dy)*R;
        hdx = degtorad(dx)*R*(cosd(mean(InR.LatitudeLimits)));
    end
end
minRatio = min(abs([hody hodx]./[hdy hdx]));
iterations = floor(log2(minRatio/thresholdRatio));

% reduce by about half every interation
for k=1:iterations
    if strcmpi(method,'nearest')
        raster = medfilt2(raster,'symmetric');
        newRaster = imresize(raster,1/2,method);
    else
        newRaster = impyramid(raster,'reduce');
        % use nearest neighbor to keep NaNs from propagating
        if any(isnan(newRaster(:)))
            X = imresize(raster,[size(newRaster,1) size(newRaster,2)],'nearest');
            t = isnan(newRaster);
            newRaster(t) = X(t);
        end
    end
    adjustment = size(raster)./size(newRaster);
    newdy = dy*adjustment(1);
    newdx = dx*adjustment(2);
    x11 = x11+(newdx-dx);
    y11 = y11+(newdy-dy);
    raster = newRaster;
    dx = newdx;
    dy = newdy;
end

% adjust raster reference if we did anything
if iterations>=1
    RM = makerefmat(x11,y11,dx,dy);
    if inProj
        newInR = refmatToMapRasterReference(RM,size(newRaster));
    else
        newInR = refmatToGeoRasterReference(RM,size(newRaster));
    end
else
    newInR = InR;
    newRaster = raster;
end
end

function [x11,y11,dx,dy,inProj] = cornerCoords(R)
% corner coordinates from raster reference
if contains(class(R),'MapCellsReference')
    inProj = true;
    [x11,y11] = intrinsicToWorld(R,1,1);
    if strcmp(R.ColumnsStartFrom,'north')
        dy = -R.CellExtentInWorldY;
    else
        dy = R.CellExtentInWorldY;
    end
    if strcmp(R.RowsStartFrom,'west')
        dx = R.CellExtentInWorldX;
    else
        dx = -R.CellExtentInWorldX;
    end
elseif contains(class(R),'GeographicCellsReference')
    inProj = false;
    [y11,x11] = intrinsicToGeographic(R,1,1);
    if strcmp(R.ColumnsStartFrom,'north')
        dy = -R.CellExtentInLatitude;
    else
        dy = R.CellExtentInLatitude;
    end
    if strcmp(R.RowsStartFrom,'west')
        dx = R.CellExtentInLongitude;
    else
        dx = -R.CellExtentInLongitude;
    end
else
    error('raster reference class %s unrecognized',class(R));
end
end

function [ InRR,OutRR,planet,method,inLat,inLon ] =...
    parseInput( sizeA,R,InS,OutS,varargin )
%parse input values to produce input/output raster references

p = inputParser;
defaultPixelSize = NaN;
defaultXLimit = NaN;
defaultYLimit = defaultXLimit;
defaultAdjust = true;
defaultOrigin = 'ul';
defaultPlanet = 'earth';
defaultWholeHemisphere = 'no';
defaultMethod = 'linear';
defaultRR = [];

addRequired(p,'sizeA',@isnumeric);
addRequired(p,'R',@isnumeric)
addRequired(p,'InS',@isstruct)
addRequired(p,'OutS',@isstruct)
addParameter(p,'pixelsize',defaultPixelSize,@isnumeric)
addParameter(p,'XLimit',defaultXLimit,@isnumeric)
addParameter(p,'YLimit',defaultYLimit,@isnumeric)
addParameter(p,'adjust',defaultAdjust,@islogical)
addParameter(p,'origin',defaultOrigin,@ischar)
addParameter(p,'planet',defaultPlanet,@ischar)
addParameter(p,'wholeHemisphere',defaultWholeHemisphere,@ischar)
addParameter(p,'method',defaultMethod,@ischar)
addParameter(p,'rasterref',defaultRR,@isobject)
addParameter(p,'lat',[],@isnumeric)
addParameter(p,'lon',[],@isnumeric)
if isempty(InS)
    InS = struct([]);
end
if isempty(OutS)
    OutS = struct([]);
end

% only need the x- & y-size
sizeA = [sizeA(1) sizeA(2)];
parse(p,sizeA,R,InS,OutS,varargin{:});

% which planet
planet = lower(p.Results.planet);

% interpolation method
method = lower(p.Results.method);

% bounding input coordinates, geographic or projected
if strcmpi(p.Results.wholeHemisphere,'north')
    lonlim = [-180 180];
    latlim = [0 90];
elseif strcmpi(p.Results.wholeHemisphere,'south')
    lonlim = [-180 180];
    latlim = [-90 0];
elseif strcmpi(p.Results.wholeHemisphere,'both')
    lonlim = [-180 180];
    latlim = [-90 90];
else
    latlim = NaN;
    lonlim = NaN;
end

inLat = double(p.Results.lat);
inLon = double(p.Results.lon);
if isempty(p.Results.R)
    assert(~isempty(p.Results.lat) && isequal(size(p.Results.lat),size(p.Results.lon)),...
        'if input RefMatrix is empty, lat/lon of same size must be specified')
    InRR = [];
    latlim = double([min(p.Results.lat(:)) max(p.Results.lat(:))]);
    lonlim = double([min(p.Results.lon(:)) max(p.Results.lon(:))]);
else
    validateattributes(R,{'numeric'},{'size',[3 2]});
    assert(isempty(p.Results.lat) && isempty(p.Results.lon),...
        'if input RefMatrix is specified, then input ''lat'' and ''lon'' must not be specified')
    if isempty(InS)
        InRR = refmatToGeoRasterReference(p.Results.R,p.Results.sizeA);
        if isnan(latlim)
            [XIntrinsic,YIntrinsic] =...
                meshgrid([1 InRR.RasterSize(2)],[1 InRR.RasterSize(1)]);
            [latlim,lonlim] = intrinsicToGeographic(InRR,XIntrinsic,YIntrinsic);
        end
    else
        InRR = refmatToMapRasterReference(p.Results.R,p.Results.sizeA);
        if isnan(latlim)
            [XIntrinsic,YIntrinsic] =...
                meshgrid([1 InRR.RasterSize(2)],[1 InRR.RasterSize(1)]);
            [xWorld,yWorld] = intrinsicToWorld(InRR,XIntrinsic, YIntrinsic);
            [latlim,lonlim] = minvtran(InS,xWorld,yWorld);
        end
    end
end

% if raster reference specified, other values default to it
% otherwise, parse other inputs
if isempty(p.Results.rasterref)
    % start of rows and cols
    if strncmpi(p.Results.origin,'l',1)
        startcol = 'south';
    else
        startcol = 'north';
    end
    if strcmpi(p.Results.origin(2),'l')
        startrow = 'west';
    else
        startrow = 'east';
    end
    
    % x- and y-limits
    if any(isnan(p.Results.XLimit(:))) || any(isnan(p.Results.YLimit(:)))
        if isempty(OutS)
            xlimit = [min(lonlim(:)) max(lonlim(:))];
            ylimit = [min(latlim(:)) max(latlim(:))];
        else
            [x,y] = mfwdtran(OutS,latlim,lonlim);
            xlimit = [min(x(:)) max(x(:))];
            ylimit = [min(y(:)) max(y(:))];
        end
    end
    if isnan(p.Results.XLimit)
        XLimit = xlimit;
    else
        assert(length(p.Results.XLimit)==2,...
            'if specified, XLimit must be vector of length 2')
        XLimit = p.Results.XLimit;
    end
    if isnan(p.Results.YLimit)
        YLimit = ylimit;
    else
        assert(length(p.Results.YLimit)==2,...
            'if specified, YLimit must be vector of length 2')
        YLimit = p.Results.YLimit;
    end
    
    % pixel size
    if isnan(p.Results.pixelsize)
        % pixel size is [height width], default is preserve image size
        pixelsize(1) = (YLimit(2)-YLimit(1))/p.Results.sizeA(1);
        pixelsize(2) = (XLimit(2)-XLimit(1))/p.Results.sizeA(2);
    else
        assert(length(p.Results.pixelsize)<=2,...
            'if specified, pixelsize must be scalar or vector length 2')
        pixelsize = p.Results.pixelsize;
        if isscalar(pixelsize)
            pixelsize(2) = pixelsize(1);
        end
    end
    
    % adjust limits to multiple of pixel size
    if p.Results.adjust
        [XLimit,YLimit] = adjustLimits(XLimit,YLimit,pixelsize);
    end
    nRows = round((YLimit(2)-YLimit(1))/pixelsize(1));
    nCols = round((XLimit(2)-XLimit(1))/pixelsize(2));
    
    % output raster references
    if isempty(OutS)
        OutRR = georasterref('LatitudeLimits',YLimit,...
            'LongitudeLimits',XLimit,'RasterSize',[nRows nCols],...
            'ColumnsStartFrom',startcol,'RowsStartFrom',startrow);
    else
        OutRR = maprasterref('YLimWorld',YLimit,...
            'XLimWorld',XLimit,'RasterSize',[nRows nCols],...
            'ColumnsStartFrom',startcol,'RowsStartFrom',startrow);
    end
else
    OutRR = p.Results.rasterref;
end

end

function RefMat = RasterRefToRefMat(R)
[x11,y11,dx,dy,~] = cornerCoords(R);
RefMat = makerefmat(x11,y11,dx,dy);
end