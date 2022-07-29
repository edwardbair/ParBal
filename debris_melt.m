function debris_melt(sFileDay,sFile,albedo,topofile,ldas_dir,...
    ldas_topo_file,ceres_dir,ceres_topo_file,merra_dir,merra_topo_file,...
    debris_depth_file,outfile,LDASOnlyFlag)
%calculate debris depth across a series of images from ebalance
%inputs
%sFileDay - which day to run
%sFile - h5 fsca file
%albedo, albedo of debris, scalar
%topofile, location of topo struct produced by TopoHorizons
%ldas_dir, path the GLDAS
%ldas_dem_file, path to GLDAS DEM and slope/aspect file struct
%ceres_dir, path to the CERES SW & LW
%ceres_topofile, path to CERES DEM and slope/aspect file struct
%merra_dir, path to MERRA wind forcings
%merra_topo_file, path to MERRA topofile

%debris_depth_file, mat file of debris depth w/ variable d in mm
%outfile, filename to write out mat files
%LDASOnlyFlag, if true use only LDAS (no CERES)
%slices
%calculate debris depth a one point in time for a directory of Tsfc images
% build topo structure,ldas filelist,and ldas topostruct

[datevalsDay,fsca,d,topo,gldas_filelist,gldas_topo,ceres,ceres_topo,...
    merra,merra_topo,tz]=...
    include_vars_debris_melt(sFileDay,sFile,debris_depth_file,...
    topofile,ldas_dir,ldas_topo_file,ceres_dir,...
    ceres_topo_file,merra_dir,merra_topo_file,LDASOnlyFlag);
%outfile=fullfile(outdir,[datestr(datevalsDay,'yyyymmdd'),'.mat']);
% need debris cover, but no snow cover, gets transfered to albedo and
% skipped by everything else when NaN

dmask=d>0 & fsca == 0;
albedo=single(dmask).*albedo;
albedo(albedo==0)=NaN;
%load interpolated forcings
[gldasInterp,ceresInterp,merraInterp]=makeInterp(gldas_filelist,...
    gldas_topo,topo,dmask,ceres,merra,tz);
dailyEnergy(topo,gldasInterp,gldas_topo,ceresInterp,...
    ceres_topo,merraInterp,merra_topo,false,false,'debris',outfile,albedo,d);