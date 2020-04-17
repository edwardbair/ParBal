function reconstructSWE(poolsize,energy_dir,sFile,rFile,varargin)
% Reconstruct SWE  and write out cubes as mat or h5 files

% INPUT
% poolsize - parpool num workers
% energy_dir - Directory for energy cube files for yr, daily
% sFile - time/space smoothed SCA cube for yr, daily
% rFile - reconstruction ouputfile, mat or h5
% optional name-value pairs:
% maxswedatesfile - peak SWE dates to limit reconstruction; otherwise zeros
% watermaskfile - binary water mask where SWE is set to zero
% canopycoverfile - static canopy cover fraction file, h5 or mat file with
% "cc" or "/Grid/cc" variable
% note: if cc is supplied and fsca projection doesn't match, fsca will be
% reprojected to match
% matdates - list of matdates to process. overrides default behavoir which
% is to reconstruct all dates in sFile
tic;
numArgs = 4;
if nargin<numArgs
    disp(['minimum ' num2str(numArgs) '(>' num2str(nargin) ') arguments required'])
    if isdeployed
        disp(['usage: ' mfilename ' poolsize energy_dir sFile outdir [optional name/value pairs: ',...
            '[maxswefile maxswefilename],[watermaskfile watermaskfilename canopycoverfile matdates matdatesvector]']);
    else
        disp(['usage: ' mfilename '(poolsize,energy_dir,sFile,outdir [optional name/value pairs: ',...
            '[''maxswefile'' maxswefilename],[''watermaskfile'' watermaskfilename canopycoverfile matdates matdatesvector])']);
    end
    disp('see documentation about optional arguments')
    return
end

% parse inputs
defaultmaxswe = [];
defaultwatermask = [];
defaultcanopycover = [];
defaultmatdates = [];

p = inputParser;
addRequired(p,'poolsize',@isnumeric);
addRequired(p,'energy_dir',@ischar)
addRequired(p,'sFile',@ischar)
addRequired(p,'rFile',@ischar)
addParameter(p,'maxswefile',defaultmaxswe,@ischar);
addParameter(p,'watermaskfile',defaultwatermask,@ischar);
addParameter(p,'canopycoverfile',defaultcanopycover,@ischar);
addParameter(p,'matdates',defaultmatdates,@isnumeric);
parse(p,poolsize,energy_dir,sFile,rFile,varargin{:})

%get matdates from sFile or from name-value pair
if isempty(p.Results.matdates)
    try
        datevalsDay=h5readatt(sFile,'/','MATLABdates');
    catch
        error('could not find record MATLABdates in %s',sFile)
    end
else
    datevalsDay=p.Results.matdates;
end

[fsca_dir,fname,ext] = fileparts(sFile);
[fsca_dir,diaryFolder] = identifyFolders(fsca_dir);
sFile = fullfile(fsca_dir,[fname,ext]);
assert(exist(sFile,'file')==2,'fsca file (%s) does not exist',sFile);

[recon_dir,fname,ext] = fileparts(rFile);
recon_dir = identifyFolders(recon_dir);
rFile = fullfile(recon_dir,[fname,ext]);

if isdeployed
    h = StartAzureDiary(mfilename,diaryFolder,...
        datestr(datevalsDay(1),'yyyymmdd'),'-',datestr(datevalsDay(end))); %#ok<NASGU>
end

%read optional files
%canopycoverfile
canopycoverfile=p.Results.canopycoverfile;
ccflag=false;
if ~isempty(canopycoverfile)
    [~,~,ext]=fileparts(canopycoverfile);
    ccflag=true;
    switch ext
        case '.mat'
            list={'cc','RefMatrix','RasterReference','ProjectionStructure'};
            for j=1:length(list)
                m=load(canopycoverfile,'-regexp',['(?i)(',list{j},')']);
                if isempty(m)
                    error('could not read %s in canopycoverfile %s',list{j},...
                        canopycoverfile);
                else
                    fn=fieldnames(m(1));
                    target_fn=list{j};
                    switch target_fn
                        case 'cc'
                            cc.(target_fn)=m.(fn{1});
                        otherwise % hdr info
                            cc.hdr.(target_fn)=m.(fn{1});
                    end
                end
            end
        case '.h5'
            try
                cc.cc=h5read(canopycoverfile,'/Grid/cc');
                cc.hdr=GetCoordinateInfo(canopycoverfile,'/Grid/',size(cc.cc));
            catch
                error(['could not find dataset in ''/Grid''',...
                    ' in h5file:',canopycoverfile]);
            end
    end
    %use cc for size
    sz=[size(cc.cc) length(datevalsDay)];
else
    %use fsca for size
    info=h5info(sFile,'/Grid/MODIS_GRID_500m/');
    sz=info.Datasets(1).Dataspace.Size;
    cc.cc=zeros([sz(1) sz(2)]);
end

%max swe date mask
maxswefile=p.Results.maxswefile;
if ~isempty(maxswefile)
    [~,~,ext]=fileparts(maxswefile);
    switch ext
        case '.mat';
            m=load(maxswefile,'-regexp','.*swe.*dates.*');
            if isempty(m)
                error('could not read maxswefile %s',maxswefile);
            end
            fn=fieldnames(m(1));
            sweR.maxswedates=m.(fn{1});
        case '.h5';
            try
                sweR.maxswedates=h5read(maxswefile,'/Grid');
            catch
                error(['could not find dataset in ''/Grid''',...
                    ' in h5file:%s',maxswefile]);
            end
    end
else
    sweR.maxswedates=ones([sz(1) sz(2)]);
end

%water mask
watermaskfile=p.Results.watermaskfile;
if ~isempty(watermaskfile)
    [~,~,ext]=fileparts(watermaskfile);
    switch ext
        case '.mat';
            m=load(watermaskfile,'-regexp','.*mask.*');
            if isempty(m)
                error('could not read watermask %s',watermaskfile);
            end
            fn=fieldnames(m(1));
            watermask=m.(fn{1});
        case '.h5';
            try
                watermask=h5read(watermaskfile,'/Grid');
            catch
                error(['could not find dataset in ''/Grid''',...
                    ' in h5file:',watermaskfile]);
            end
    end
else
    watermask=false([sz(1) sz(2)]);
end

watermask=repmat(logical(watermask),[1 1 sz(3)]);


%start parallel pool
poolsize=p.Results.poolsize;
parpool_check(poolsize);

%constants
% Ci=1.9e3; %specific heat of ice, J/kg deg
% rho_ice=917; %density of ice, kg/m^3
% Tg=273; %temp of ground, K
% Tf=273; % freezing temp, K
Lf=3.34e5;  % latent heat of fusion, J kg^-1,
rho_water=1000;% density of water kg/m^3
sec_hr=3600; %seconds in an hour
mf=1/(rho_water*Lf)*sec_hr*1000; %melt factor, mm water/(W m^-2)
melt=zeros(sz,'single');

%compute pot melt
parfor d=1:length(datevalsDay)
    t=datevalsDay(d);
    fprintf('reconstructing %s\n',datestr(t));
    %load sca
    [rawsca,~,hdr]=GetEndmember(sFile,'snow',t);
    %if a canopy cover file was supplied and it has a different header,
    %reproject
    if ccflag && ~isequal(hdr,cc.hdr)
        rawsca=reprojectRaster(rawsca,hdr.RefMatrix,hdr.ProjectionStructure,...
            cc.hdr.ProjectionStructure,'rasterref',cc.hdr.RasterReference);
    end
    cctmp=cc.cc;
    ind=rawsca+cctmp > 1;
    cctmp(ind)=1-rawsca(ind);
    sca=rawsca./(1-cctmp);
    sca(rawsca==0)=0;
    %quick fix for binary fsca
%     sca(sca>0)=1;
    %load potential melt
    fname=fullfile(energy_dir,[datestr(t,'yyyymmdd'),'.mat']);
    %     [M,MATLABdates]=parload(fname,'M','MATLABdates');
    m=matfile(fname);
    M=m.M;
    melt(:,:,d)=single(M).*mf.*sca;
end
%% sum SWE
%allocate
rsize=[sz(1) sz(2)];
vecsize=rsize(1)*rsize(2);
% convert to uint16
sweR.melt=uint16(melt);
sweR.melt(watermask)=0;

[sca,~,hdr]=GetEndmember(sFile,'snow_fraction',datevalsDay);
% have to reproject again since last time was done in parfor loop
% if a canopy cover file was supplied and it has a different header,
%reproject
if ccflag && ~isequal(hdr,cc.hdr)
    sca=reprojectRaster(sca,hdr.RefMatrix,hdr.ProjectionStructure,...
        cc.hdr.ProjectionStructure,'rasterref',cc.hdr.RasterReference);
    hdr=cc.hdr;
end
swe=zeros([length(datevalsDay) vecsize],'single');

%reshape inputs
%transpose is req'd because reshape operates on columns
%rows=time, cols=image vector
melt=reshape(sweR.melt,vecsize,length(datevalsDay))';
sca=reshape(sca,vecsize,length(datevalsDay))';
maxswedates=reshape(sweR.maxswedates,1,vecsize);

%one pixel at every timestep
parfor k=1:vecsize;
    melt_vec=single(melt(:,k));
    sca_vec=sca(:,k);
    swe_temp=zeros(length(datevalsDay),1);
    if maxswedates(k) > 0 && any(sca_vec);
        % look for contiguous sca pds
        % eliminate spurious zeros in SCA
        %         sca_vec=slidefun(@max,7,sca_vec,'central');
        [start,finish] = contiguous(sca_vec);
        %get rid of contiguous pds w/ start and finish prior to peak
        %whats left are contiguous pds that start before or on peak
        ind=start < maxswedates(k) & finish < maxswedates(k);
        start(ind)=[];
        finish(ind)=[];
        for h=1:length(start);
            s=start(h);
            f=finish(h);
            % same as cumsum from s to f
            for i=s:f
                swe_temp(i)=sum(melt_vec(i:f));
            end
        end
        %set SWE to max for all dates prior to max
        swe_temp(1:maxswedates(k))=swe_temp(maxswedates(k));
    end
    swe(:,k)=swe_temp;
end

% go back to images/cubes
sweR.swe=reshape(swe',rsize(1),rsize(2),length(datevalsDay));
sweR.swe=uint16(sweR.swe);
%% write out files
%get extension
[~,~,ext]=fileparts(rFile);

fn=fieldnames(sweR);
units='mm';

switch ext
    case '.mat'
        save(rFile,'-struct','sweR','-v7.3');
        save(rFile','units','datevalsDay','hdr','-append');
    case '.h5'
        deflateLevel=9;
        if ~isempty(dir(rFile))
            delete(rFile);
        end
        for i=1:length(fn)
            sz=size(sweR.(fn{i}));
            if length(sz) > 2
                ChunkSize=[sz(1:2) 1];
            else
                ChunkSize=sz;
            end
            location=['/Grid/' fn{i}];
            h5create(rFile,location,sz,'Deflate',deflateLevel,'ChunkSize',...
                ChunkSize,'DataType',class(sweR.(fn{i})));
            h5write(rFile,location,sweR.(fn{i}));
            %write mm for units for melt and swe
            if strcmpi(fn{i},'melt') || strcmpi(fn{i},'swe')
                h5writeatt(rFile,location,'units',units);
            end
        end
        h5writeProjection(rFile,'/Grid',hdr.ProjectionStructure);
        h5writeatt(rFile,'/','MATLABdates',datevalsDay);
        %added yyyymmdd attribute for non-MATLAB users
        h5writeatt(rFile,'/','ISOdates',...
            strjoin(cellstr(datestr(datevalsDay,'yyyymmdd'))));
        h5writeatt(rFile,'/Grid','ReferencingMatrix',hdr.RefMatrix);
end
t=toc;
fprintf('SWE reconstructed and saved in %3.1f min\n',t/60);