function downscale_energy(sFileDay,sFile,poolsize,topofile,landcoverfile,...
    ldas_dir,ldas_dem_file,ceres_dir,ceres_topofile,fast_flag,outdir)
% downscales energy for a day, specifically meant for running on individual
% workers on Azure. Contains a parfor loop for computing snow surface
% temperatures
% input:
% sFileDay - scalar value corresponding to day (slice) of sFile
% poolsize - size of parallel pool
% sFile - path to h5 smoothed fSCA cube that must contain: snow_fraction,grain_size,and
% deltavis
% topofile - path to h5 topofile created by TopoHorizons
% landcoverfile - path to h5 or mat land cover file
% ldas_dir - base path to ldas files, don't specify a year
% ldas_dem_file - path to h5 or mat ldas dem
% ceres_dir - base path to ceres files
% ceres_topofile - path to ceres topofile (mat - h5 not implemented yet)
% fast flag - true - only solve for M; false - solve for all outputs; only set for 'normal'
% outdir - path to write out files, must already be created

numarg=7;

if isdeployed && nargin~=numarg
    disp(['usage: ' mfilename ' sFileDay sFile topofile landcoverfile',...
    'ldas_dir ldas_dem_dir outdir']);
    return
end
p = inputParser;
validationFcn = @(x) isnumeric(x) || ischar(x);
addRequired(p,'sFileDay',validationFcn)
addRequired(p,'sFile',validationFcn)
addRequired(p,'poolsize',validationFcn)
addRequired(p,'topofile',validationFcn)
addRequired(p,'landcoverfile',validationFcn)
addRequired(p,'ldas_dir',validationFcn)
addRequired(p,'ldas_dem_file',validationFcn)
addRequired(p,'ceres_dir',validationFcn)
addRequired(p,'ceres_topofile',validationFcn)
addRequired(p,'outdir',validationFcn)

parse(p,sFileDay,sFile,poolsize,topofile,landcoverfile,ldas_dir,ldas_dem_file,...
    ceres_dir,ceres_topofile,outdir)

% variables
if ~isnumeric(p.Results.sFileDay)
    sFileDay = str2double(p.Results.sFileDay);
end

% get required variables and filenames
[FOREST, topo, ldas, ldas_topo, dateval, sFile, ceres,...
    ceres_topo,tz,outfile]=include_vars_melt(sFileDay,sFile,poolsize,topofile,...
    landcoverfile,ldas_dir,ldas_dem_file,ceres_dir,ceres_topofile,outdir);

%run downscaling
daily_melt(dateval,ldas,ceres,topo,ldas_topo,ceres_topo,...
    tz,FOREST,sFile,outfile,fast_flag)