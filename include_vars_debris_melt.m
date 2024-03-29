function [datevalsDay,fsca,d,topo,ldas,ldas_topo,ceres,ceres_topo,merra,...
    merra_topo,tz]=include_vars_debris_melt(sFileDay,sFile,debris_depth_file,...
topofile,ldas_dir,ldas_topofile,ceres_dir,ceres_topofile,merra_dir,merra_topofile,...
LDASOnlyFlag)

%sca file
datevals=h5readatt(sFile,'/','MATLABdates');
datevalsDay=datevals(sFileDay);
[fsca,~,fsca_hdr]=GetEndmember(sFile,'snow',datevalsDay);

%debris_depth
m=matfile(debris_depth_file);
d=single(m.debris_depth)/1000;%convert from mm to m

%topofile
[topo.slope,topo.hdr]=GetTopography(topofile,'slope');
topo.aspect=GetTopography(topofile,'aspect');
topo.dem=GetTopography(topofile,'elevation');
topo.view=GetTopography(topofile,'viewfactor');
topo.topofile=topofile;

% compute curvature
curve_len_scale=2000;%m
topo.curvature=curvature(topo.dem,topo.hdr.RasterReference,...
    curve_len_scale);

%reproject fsca if it doesn't match topo
if any(topo.hdr.RefMatrix~=fsca_hdr.RefMatrix,'all') 
    fsca=reprojectRaster(fsca,fsca_hdr.RefMatrix,...
        fsca_hdr.ProjectionStructure,...
            topo.hdr.ProjectionStructure,...
            'rasterref',topo.hdr.RasterReference);
end

%use midpoints
x=mean(topo.hdr.RasterReference.XWorldLimits);
y=mean(topo.hdr.RasterReference.YWorldLimits);
[~,lon]=minvtran(topo.hdr.ProjectionStructure,x,y);

%convert to - for west of PM and as fraction of 24 hr
tz=-timezone(lon)/24;
ldas=make_ldas_filelist(datevalsDay,ldas_dir,tz);

%LDAS variable list
ldas.var={'Tair_f_inst','Psurf_f_inst','Qair_f_inst'};
ldas_topo_names={'Z','aspect','slope'};
ldas_topo=load_coarse_topo(ldas_topofile,ldas_topo_names,topo);


if ~LDASOnlyFlag
    ceres.var={'adj_atmos_sw_down_all_surface_1h',...
        'adj_atmos_lw_down_all_surface_1h'};
    ceres.ceres_dir=ceres_dir;
    ceres_topo=load_coarse_topo(ceres_topofile,ldas_topo_names,topo);   
    
    merra.merra_dir=merra_dir;
    merra.var={'ULML','VLML'};
    merra_topo_names={'Z','aspect','slope'};
    merra_topo=load_coarse_topo(merra_topofile,merra_topo_names,topo);
else
    ceres=[];
    ceres_topo=[];
    merra=[];
    merra_topo=[];
end