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

%crop GLDAS first to reduce memory usage during reprojection
%note fliplr to get upper left corner for y
[lat,lon]=minvtran(fine_dem_hdr.ProjectionStructure,...
    fine_dem_hdr.RasterReference.XWorldLimits,...
    fliplr(fine_dem_hdr.RasterReference.YWorldLimits));

% [x,y]=mfwdtran(hdr.ProjectionStructure,lat,lon);

% [r,c]=map2pix(hdr.RefMatrix,x,y);
[r,c]=latlon2pix(coarseRefMat,lat,lon);


border=1.25;
by=border*fine_dem_hdr.RasterReference.RasterSize(1);
bx=border*fine_dem_hdr.RasterReference.RasterSize(2);

r=round([r(1)-by r(2)+by]);
r(r<1)=1;r(r>size(GLDAS.(vars{1}),1))=size(GLDAS.(vars{1}),1);

c=round([c(1)-bx c(2)+bx]);
c(c<1)=1;c(c>size(GLDAS.(vars{1}),2))=size(GLDAS.(vars{1}),2);

[x,y]=pixcenters(coarseRefMat,size(GLDAS.(vars{1})));

x=x(c(1):c(2));
y=y(r(1):r(2));
hdr.RefMatrix=makerefmat(x(1),y(1),x(2)-x(1),y(2)-y(1));
% hdr.RasterReference=refmatToGeoRasterReference(hdr.RefMatrix,...
%     size(GLDAS.(vars{1})));
% v=v(r(1):r(2),c(1):c(2));
% v=single(v);
parfor i=1:length(vars)

    %crop
    GLDAS.(vars{i})=GLDAS.(vars{i})(r(1):r(2),c(1):c(2));    
    GLDAS.(vars{i})=single(GLDAS.(vars{i}));

    gldastmp{i}=reprojectRaster(GLDAS.(vars{i}),hdr.RefMatrix,...
        [],fine_dem_hdr.ProjectionStructure,'rasterref',...
        fine_dem_hdr.RasterReference,'method','linear');

%    gldastmp{i}=single(t);
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
