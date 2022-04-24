function [FOREST, topo, ldas, ldas_topo, dateval, sFile, ceres, ...
    ceres_topo, merra, merra_topo, tz,outfile]=include_vars_melt(sFileDay,sFile,topofile,...
    landcoverfile,ldas_dir,ldas_dem_file,ceres_dir,ceres_topofile,merra_dir,merra_topofile,...
    outdir,LDASOnlyFlag)
% for each day or set of days
% assemble variables needed to downscale energy balance

%convert to Azure paths
[fsca_dir,fname,ext] = fileparts(sFile);
[fsca_dir,diaryFolder] = identifyFolders(fsca_dir);
sFile = fullfile(fsca_dir,[fname,ext]);
% get dates corresponding to sFileDay
try
    datevals=h5readatt(sFile,'/','MATLABdates');
catch
    error('no dates in %s',sFile);
end
dateval=floor(datevals(sFileDay)); %start at midnight

if isdeployed
    h = StartAzureDiary('downscale_energy_azure',diaryFolder,...
        datestr(dateval,'yyyymmdd')); %#ok<NASGU>
end

assert(exist(sFile,'file')==2,'Fsca file (%s) does not exist',sFile);

[topo_dir,fname,ext] = fileparts(topofile);
topo_dir = identifyFolders(topo_dir);
topofile = fullfile(topo_dir,[fname,ext]);
assert(exist(topofile,'file')==2,'Topo file (%s) does not exist',topofile);

[landcover_dir,fname,ext] = fileparts(landcoverfile);
landcover_dir = identifyFolders(landcover_dir);
landcoverfile = fullfile(landcover_dir,[fname,ext]);
assert(exist(landcoverfile,'file')==2,'Land Cover file (%s) does not exist',landcoverfile);

ldas_dir = identifyFolders(ldas_dir);

[ldas_dem_dir,fname, ext] = fileparts(ldas_dem_file);
ldas_dem_dir = identifyFolders(ldas_dem_dir);
ldas_dem_file = fullfile(ldas_dem_dir,[fname,ext]);
assert(exist(ldas_dem_file,'file')==2,'LDAS DEM file (%s) does not exist',ldas_dem_file);

if ~LDASOnlyFlag
    [ceres_topofile_dir,fname, ext] = fileparts(ceres_topofile);
    ceres_topofile_dir = identifyFolders(ceres_topofile_dir);
    ceres_topofile = fullfile(ceres_topofile_dir,[fname,ext]);
    assert(exist(ceres_topofile,'file')==2,'CERES topofile (%s) does not exist',ceres_topofile);

    [merra_topofile_dir,fname, ext] = fileparts(merra_topofile);
    merra_topofile_dir = identifyFolders(merra_topofile_dir);
    merra_topofile = fullfile(merra_topofile_dir,[fname,ext]);
    assert(exist(merra_topofile,'file')==2,'MERRA topofile (%s) does not exist',merra_topofile);
end

outdir = identifyFolders(outdir);
% error if energy file exists
outfile=fullfile(outdir,[datestr(dateval,'yyyymmdd'),'.mat']);

%topography stuff
[topo.slope,topo.hdr]=GetTopography(topofile,'slope');
topo.aspect=GetTopography(topofile,'aspect');
topo.dem=GetTopography(topofile,'elevation');
topo.view=GetTopography(topofile,'viewfactor');
topo.topofile=topofile;

% compute curvature
curve_len_scale=2000;%m
topo.curvature=curvature(topo.dem,topo.hdr.RasterReference,...
    curve_len_scale);

%landcover
[~,~,ext]=fileparts(landcoverfile);
switch ext
    case '.mat'
        m=load(landcoverfile);
        if isempty(m)
            error('could not read landcover file %s',landcoverfile);
        end
        fn=fieldnames(m);
        for i=1:length(fn)
            LandCover.(fn{i})=m.(fn{i});
        end
    case '.h5'
        LandCover.Z=h5read(landcoverfile,'/Grid/Z');
        LandCover.cc=h5read(landcoverfile,'/Grid/cc');
    otherwise
        error('landcover file %s not .h5 or .mat',landcoverfile);
end

FOREST = makeFOREST(LandCover);

%get timezone based on longitude

%use midpoints
x=mean(topo.hdr.RasterReference.XWorldLimits);
y=mean(topo.hdr.RasterReference.YWorldLimits);
[~,lon]=minvtran(topo.hdr.ProjectionStructure,x,y);

%convert to - for west of PM and as fraction of 24 hr
tz=-timezone(lon)/24;
ldas=make_ldas_filelist(dateval,ldas_dir,tz);

%LDAS variable list
if contains(ldas_dem_dir,'NLDAS')
    ldas.var ={'TMP','PRES','UGRD','VGRD','SPFH'};
elseif contains(ldas_dem_dir,'GLDAS')
    ldas.var={'Tair_f_inst','Psurf_f_inst',...
        'Wind_f_inst','Qair_f_inst','SWE_inst'};
    %add in radiation if LDAS only
    if LDASOnlyFlag
        ldas.var=['SWdown_f_tavg','LWdown_f_tavg',ldas.var];
    end
end

ldas_topo_names={'Z','aspect','slope'};
ldas_topo=load_coarse_topo(ldas_dem_file,ldas_topo_names,topo);

%add in CERES data
if ~LDASOnlyFlag
    %ceres eliminated the surface pressure variable in ed 4a
    ceres.var={'adj_atmos_sw_down_all_surface_1h',...
        'adj_atmos_lw_down_all_surface_1h'};
    %currently just fixing names in 4.1, see fix41names.m
    % if dateval >= datenum([2018 10 01]) %dash becomes underscore in ed 4.1
    %     ceres.var=strrep(ceres.var,'-','_');
    % end

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
