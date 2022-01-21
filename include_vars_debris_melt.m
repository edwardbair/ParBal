function [datevalsDay,fsca,d,topo,ldas,ldas_topo,ceres,ceres_topo,tz]=...
include_vars_debris_melt(sFileDay,sFile,debris_depth_file,...
topofile,ldas_dir,ceres_dir,ceres_topofile,ldas_topofile,LDASOnlyFlag)


%sca file
datevals=h5readatt(sFile,'/','MATLABdates');
datevalsDay=datevals(sFileDay);
[fsca,~,fsca_hdr]=GetEndmember(sFile,'snow',datevalsDay);

%debris_depth
m=matfile(debris_depth_file);
d=single(m.d)/1000;%convert from mm to m

%topofile
[topo.slope,topo.hdr]=GetTopography(topofile,'slope');
topo.aspect=GetTopography(topofile,'aspect');
topo.dem=GetTopography(topofile,'elevation');
topo.view=GetTopography(topofile,'viewfactor');
topo.topofile=topofile;

%reproject fsca if it doesn't match topo
if any(topo.hdr.RefMatrix~=fsca_hdr.RefMatrix,'all') 
    fsca=rasterReprojection(fsca,fsca_hdr.RefMatrix,fsca_hdr.ProjectionStructure,...
        topo.hdr.ProjectionStructure,'rasterref',topo.hdr.RasterReference);
end

%use midpoints
x=mean(topo.hdr.RasterReference.XWorldLimits);
y=mean(topo.hdr.RasterReference.YWorldLimits);
[~,lon]=minvtran(topo.hdr.ProjectionStructure,x,y);

%convert to - for west of PM and as fraction of 24 hr
tz=-timezone(lon)/24;
ldas=make_ldas_filelist(datevalsDay,ldas_dir,tz);

%LDAS variable list
% ldas.vars ={'TMP','PRES','WIND','SPFH'};
ldas.var={'Tair_f_inst','Psurf_f_inst',...
    'Wind_f_inst','Qair_f_inst'};
ldas_topo_names={'Z','aspect','slope'};
ldas_topo=load_coarse_topo(ldas_topofile,ldas_topo_names,topo);

%add in CERES data
if ~LDASOnlyFlag
%ceres eliminated the surface pressure variable in ed 4a
ceres.var={'adj_atmos_sw-down_all_surface_1h',...
        'adj_atmos_lw-down_all_surface_1h'};
%currently just fixing names in 4.1, see fix41names.m
% if dateval >= datenum([2018 10 01]) %dash becomes underscore in ed 4.1
%     ceres.var=strrep(ceres.var,'-','_');
% end
    
ceres.ceres_dir=ceres_dir;
ceres_topo=load_coarse_topo(ceres_topofile,ldas_topo_names,topo);
else
    ceres=[];
    ceres_topo=[];
end