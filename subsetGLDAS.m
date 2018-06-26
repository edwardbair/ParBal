function gldas = subsetGLDAS(GLDAS,coarseRefMat,fine_dem_hdr,vars,mask,...
    inpaint)
%creates reprojected subset of gldas data in fine dem bounding box and projection
%INPUTS
%GLDAS - GLDAS structure from stackGLDAS_modis
%coarseRefMat - coarse GLDAS referencing matrix, e.g. 1/4 deg, geog, world
%fine_dem_hdr - header for fine scale dem with
%RefMatrix,ProjectionStructure, and RasterReference
%vars - list of gldas variables you want to reproject
%mask - mask of areas not to interpolate
%inpaint - use inpaint nans on GLDAS variables (true/false)
%OUTPUT
% gldas - structure of subsetted GLDAS variables  
gldas=GLDAS;
len=length(GLDAS.datevalsUTC);
gldastmp=cell(length(vars),1);
% sca=GetEndmember(sFile,'snow',floor(GLDAS.datevalsLocal(1)));
mask=repmat(mask,[1 1 len]);
% for i=1:length(vars);
parfor i=1:length(vars)
    t=reprojectRaster(GLDAS.(vars{i}),coarseRefMat,...
        [],fine_dem_hdr.ProjectionStructure,'rasterref',...
        fine_dem_hdr.RasterReference,'method','linear');
    gldastmp{i}=single(t);
    %dont bother interpolating for NaNs in areas w/o snow
    ind = mask > 0 & isnan(t);
    if inpaint && any(ind(:))
        for j=1:len
            V=double(gldastmp{i}(:,:,j));
            [X,Y]=meshgrid(1:size(V,1),1:size(V,2));
            ind2=ind(:,:,j);
            f=scatteredInterpolant(X(~ind2),Y(~ind2),V(~ind2));
            V(ind2)=f(X(ind2),Y(ind2));
            gldastmp{i}(:,:,j)=single(V);
%             x=interp2(X(~ind2),Y(~ind2),V(~ind2),X(ind2),Y(ind2));
%             gldastmp{i}(:,:,j)=single(x);
        end
    end
end
for i=1:length(vars)
    gldas.(vars{i})=gldastmp{i};
end
