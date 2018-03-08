function debris_depth(albedo,topofile,ldas_dir,ldas_topo_file,ceres_dir,...
    ceres_topo_file,Tsfc_dir,dmaskfile,outdir)
%calculate debris depth across a series of images from ebalance
%inputs
%albedo, albedo of debris, scalar
%topofile, location of topo struct produced by TopoHorizons
%ldas_dir, path the GLDAS
%ldas_dem_file, path to GLDAS DEM and slope/aspect file struct
%ceres_dir, path to the CERES SW & LW
%ceres_topofile, path to CERES DEM and slope/aspect file struct
%Tsfc_dir, directory of AST08 Tsfc TIFs from ASTER, will loop through
%dmaskfile, mat file of binary debris mask, same size as topofile struct
%outdir, directory to write out mat files
%calculate debris depth a one point in time for a directory of Tsfc images
d=dir(fullfile(Tsfc_dir,'*.tif'));

% build topo structure,ldas filelist,and ldas topostruct
m=matfile(dmaskfile);
dmask=m.dmask;
%start parallel pool
% parpool_check(poolsize);
for i=1:length(d)
    % for i=1:length(d)
    fname_base=d(i).name;
    Tsfc_filedate=datenum(fname_base(11:24),'mmddyyyyhhMMss');   
    [topo,gldas_filelist,gldas_topo,ceres,ceres_topo,tz]=include_vars_debris(...
        topofile,ldas_dir,ceres_dir,ceres_topo_file,ldas_topo_file,floor(Tsfc_filedate));
    fname=fullfile(d(i).folder,fname_base);
    %get cloud cover from line x of metfile
    linenum = 345;
    cloudthresh = 80;
    metname=[fname,'.met'];
    fid=fopen(metname);
    C = textscan(fid,'%s',1,'delimiter','\n', 'headerlines',linenum-1);
    fclose(fid);
    val=regexp(C{1},'\d*','match');
    val=str2double(cell2mat(val{1}));
    % flag bad pixels
    qaname=[fname(1:end-22),'QA_DataPlane.txt'];
    c=dlmread(qaname,' ');
    good=c==0;
    if val < cloudthresh
        [Tsfc,~,r,~]=geotiffread(fname);
        Tsfc=single(Tsfc);
        Tsfc(Tsfc==2000)=NaN;
        Tsfc(~good)=NaN;
        Tsfc=Tsfc*0.1;
        Tsfc(Tsfc<273.15)=NaN; % cannot map negative values
        info=geotiffinfo(fname);
        mstruct=geotiff2mstruct(info);
        %AST08 TIFFS are rotated, so need a geolocated input to reprojectRaster
        [x,y] = pixcenters(r,size(Tsfc));
        [lat,lon] = minvtran(mstruct,x,y);
        Tsfc=reprojectRaster(Tsfc,[],[],...
            topo.hdr.ProjectionStructure,'lat',lat,...
            'lon',lon,'rasterref',...
            topo.hdr.RasterReference);
        %if no values are non NaN (out of basin), skip
        if ~any(~isnan(Tsfc(:)))
            continue;
        end
        albedo_=single(dmask).*albedo;
        albedo_(albedo_==0)=NaN;
        %load interpolated forcings
        [gldasInterp,ceresInterp]=makeInterp(gldas_filelist,...
            gldas_topo,topo,dmask,ceres,tz);
        %now extract inputs for time of day when Tsfc image is from
        [~,ind]=min(abs(gldasInterp.datevalsUTC-Tsfc_filedate));
        gldasInterp=extractInd(gldasInterp,ind);
        %and for ceres
        ceresInterp=extractInd(ceresInterp,ind);
        % solve for debris depth
        outfile=fullfile(outdir,[fname_base(1:end-4),'.mat']);
        dailyEnergy(topo,gldasInterp,gldas_topo,ceresInterp,...
            ceres_topo,false,'debris depth',outfile,albedo_,Tsfc);
    end
end