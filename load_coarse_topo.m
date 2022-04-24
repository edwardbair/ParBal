function coarse_topo=load_coarse_topo(...
    topofile,topo_names,fine_topo)
%loads coarse topography
%inputs:
%topofile - h5 or matfile with coarse topograpy
%topo_names - names of coarse topography vars
%fine_topo - fine topo struct

%output
%coarse_topo - reprojected coarse topo w/ Zdiff
%ReferencingMatrix - coarse Referencing Matrix
%load topofile
geolocated=false;

[~,~,ext]=fileparts(topofile);

switch ext
    case '.mat'
        m=load(topofile);
        if isempty(m)
            error('could not read ldas file %s',topofile);
        end
        if isfield(m,'lat')
            lat=m.lat;
            lon=m.lon;
            geolocated=true;
        else
            ReferencingMatrix=m.ReferencingMatrix;
        end
        for i=1:length(topo_names)
            coarse_topo.(topo_names{i})=m.(topo_names{i});
        end
        
    case '.h5'
        for i=1:length(topo_names)
            loc=['/Grid/',topo_names{i}];
            coarse_topo.(topo_names{i})=h5read(topofile,loc);
        end
        ReferencingMatrix=h5readatt(topofile,'/Grid',...
            'ReferencingMatrix');
    otherwise
        error('topo file %s not .mat or .h5',topofile);
end

% reproject coarse elevation, slope, and aspect to fine scale
for i=1:length(topo_names)
    if geolocated
    coarse_topo.(topo_names{i})=reprojectRaster(coarse_topo.(topo_names{i}),...
        [],[],fine_topo.hdr.ProjectionStructure,'lat',lat,'lon',lon,...
    'rasterref',fine_topo.hdr.RasterReference);
    else
    coarse_topo.(topo_names{i})=reprojectRaster(coarse_topo.(topo_names{i}),...
        ReferencingMatrix,[],fine_topo.hdr.ProjectionStructure,...
    'rasterref',fine_topo.hdr.RasterReference);
    end
end

% Difference between fine and reprojected  DEM
coarse_topo.Zdiff=fine_topo.dem-coarse_topo.Z;
% coarse_topo=rmfield(coarse_topo,'Z');
if ~geolocated
    coarse_topo.CoarseRefMatrix=ReferencingMatrix;
end
end