function ceres=addGLDAStoCERES(ceres,ceres_hdr,gldas,topo,addvars)
%adds in  GLDAS variables to ceres data structure
%requires scaling up gldas data to 1 deg
%input
%ceres - ceres data struct
%ceres_topo - ceres topo struct
%gldas - gldas struct
%gldas_topo - gldas topo struct
%topo - fine scale topographic struct
%addvars - Nx1 cell of names of variables to be added from gldas to ceres, e.g. 
% {'TMP','PRES'}
%output
%appended ceres struct

%go up to 1 deg
for i=1:length(addvars)
    x=reprojectRaster(gldas.(addvars{i}),topo.hdr.RefMatrix,...
        topo.hdr.ProjectionStructure,[],'rasterref',ceres_hdr.RasterReference);
    %back down to fine
    x=reprojectRaster(x,ceres_hdr.RefMatrix,[],topo.hdr.ProjectionStructure,...
        'rasterref',topo.hdr.RasterReference);
    ceres.(addvars{i})=x;
end

