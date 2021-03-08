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
% note: if any optional input is supplied and fsca projection doesn't match, fsca will be
% reprojected to match
% matdates - list of matdates to process. overrides default behavoir which
% is to reconstruct all dates in sFile
% binarymodevalue - transform fsca into binary using specified threshold. 
% For point validation. Default is [], which does not transform fsca.
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
defaultbinarymodevalue = []; 

p = inputParser;
addRequired(p,'poolsize',@isnumeric);
addRequired(p,'energy_dir',@ischar)
addRequired(p,'sFile',@ischar)
addRequired(p,'rFile',@ischar)
addParameter(p,'maxswefile',defaultmaxswe,@ischar);
addParameter(p,'watermaskfile',defaultwatermask,@ischar);
addParameter(p,'canopycoverfile',defaultcanopycover,@ischar);
addParameter(p,'matdates',defaultmatdates,@isnumeric);
addParameter(p,'binarymodevalue',defaultbinarymodevalue,@isnumeric);

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
        datestr(datevalsDay(1),'yyyymmdd'),'-',datestr(datevalsDay(end))); 
end

%set default target hdr from fsca, change if optional inputs don't match
gname='/Grid/MODIS_GRID_500m/';
info=h5info(sFile,gname);
sz=info.Datasets(1).Dataspace.Size;
sz(3)=length(datevalsDay);
target_hdr=GetCoordinateInfo(sFile,gname,[sz(1) sz(2)]);

%read optional files
%canopycoverfile
canopycoverfile=p.Results.canopycoverfile;

if ~isempty(canopycoverfile)
    [~,~,ext]=fileparts(canopycoverfile);
    switch ext
        case '.mat'
                m=load(canopycoverfile);
                if isempty(m)
                    error('could not read canopycoverfile %s',...
                        canopycoverfile);
                else
                    cc=m.cc;
                    target_hdr=m.hdr;
                end
        case '.h5'
            try
                cc=h5read(canopycoverfile,'/Grid/cc');
                target_hdr=GetCoordinateInfo(canopycoverfile,'/Grid/',size(cc));
            catch
                error(['could not find dataset in ''/Grid''',...
                    ' in h5file:',canopycoverfile]);
            end
    end
end

%max swe date mask
maxswefile=p.Results.maxswefile;
if ~isempty(maxswefile)
    [~,~,ext]=fileparts(maxswefile);
    switch ext
        case '.mat'
            m=load(maxswefile);
            if isempty(m)
                error('could not read maxswefile %s',maxswefile);
            end
            sweR.maxswedates=m.maxswedates;
            target_hdr=m.hdr;
            
        case '.h5'
            try
                sweR.maxswedates=h5read(maxswefile,'/Grid');
                target_hdr=GetCoordinateInfo(canopycoverfile,'/Grid/',...
                    sweR.maxswedates);
            catch
                error(['could not find dataset in ''/Grid''',...
                    ' in h5file:%s',maxswefile]);
            end
    end
end



%water mask
watermaskfile=p.Results.watermaskfile;
if ~isempty(watermaskfile)
    [~,~,ext]=fileparts(watermaskfile);
    switch ext
        case '.mat'
            m=load(watermaskfile);
            if isempty(m)
                error('could not read watermask %s',watermaskfile);
            end
            watermask=m.mask;
            target_hdr=m.hdr;
        case '.h5'
            try
                watermask=h5read(watermaskfile,'/Grid');
                target_hdr=GetCoordinateInfo(watermaskfile,'/Grid',...
                    size(watermask));
            catch
                error(['could not find dataset in ''/Grid''',...
                    ' in h5file:',watermaskfile]);
            end
    end
end
%now that target_hdr may have been changed from default,
%set all optional inputs that were not supplied to that size
if isempty(maxswefile)
   sweR.maxswedates=ones(target_hdr.RasterReference.RasterSize); 
end
if isempty(watermaskfile)
   watermask=zeros(target_hdr.RasterReference.RasterSize);
end
if isempty(canopycoverfile)
   cc=zeros(target_hdr.RasterReference.RasterSize); 
end

watermask=repmat(logical(watermask),[1 1 sz(3)]);

%binarymodevalue
binarymodevalue=p.Results.binarymodevalue;

%start parallel pool
poolsize=p.Results.poolsize;
parpool_check(poolsize);

%constants
Lf=3.34e5;  % latent heat of fusion, J kg^-1,
rho_water=1000;% density of water kg/m^3
sec_hr=3600; %seconds in an hour
mf=1/(rho_water*Lf)*sec_hr*1000; %melt factor, mm water/(W m^-2)

sz=target_hdr.RasterReference.RasterSize;
sz(3)=length(datevalsDay);

melt=zeros(sz,'single');
fsca=zeros(sz,'single');

%compute pot melt
parfor d=1:length(datevalsDay)
    t=datevalsDay(d);
    fprintf('reconstructing %s\n',datestr(t));
    %load sca
    [rawsca,~,hdr]=GetEndmember(sFile,'snow',t);
    %if a canopy cover file was supplied and it has a different header,
    %reproject
    if ~(isequal(hdr.RefMatrix,target_hdr.RefMatrix) &&...
            isequal(hdr.RasterReference.RasterSize,...
            target_hdr.RasterReference.RasterSize))
        rawsca=reprojectRaster(rawsca,hdr.RefMatrix,hdr.ProjectionStructure,...
            target_hdr.ProjectionStructure,'rasterref',target_hdr.RasterReference);
    end
    cctmp=cc;
    ind=rawsca+cctmp > 1;
    cctmp(ind)=1-rawsca(ind);
    rawsca=rawsca./(1-cctmp);
    rawsca(rawsca==0)=0;
    
    if ~isempty(binarymodevalue)
       bmask=rawsca>=binarymodevalue;
       rawsca(:)=0;
       rawsca(bmask)=1;
    end
    
    fsca(:,:,d)=rawsca;
    %load potential melt
    fname=fullfile(energy_dir,[datestr(t,'yyyymmdd'),'.mat']);
    %     [M,MATLABdates]=parload(fname,'M','MATLABdates');
    m=matfile(fname);
    M=m.M;
    melt(:,:,d)=single(M).*mf.*rawsca;
end
%% sum SWE
%allocate
rsize=[sz(1) sz(2)];
vecsize=rsize(1)*rsize(2);
% convert to uint16
sweR.melt=uint16(melt);
sweR.melt(watermask)=0;

% [sca,~,hdr]=GetEndmember(sFile,'snow_fraction',datevalsDay);
% % have to reproject again since last time was done in parfor loop
% % if a canopy cover file was supplied and it has a different header,
% %reproject
% if ~isequal(hdr.RefMatrix,target_hdr.RefMatrix)
%     sca=reprojectRaster(sca,hdr.RefMatrix,hdr.ProjectionStructure,...
%         cc_hdr.ProjectionStructure,'rasterref',cc_hdr.RasterReference);
%     hdr=cc_hdr;
% end

swe=zeros([length(datevalsDay) vecsize],'single');

%reshape inputs
%transpose is req'd because reshape operates on columns
%rows=time, cols=image vector
melt=reshape(sweR.melt,vecsize,length(datevalsDay))';
fsca=reshape(fsca,vecsize,length(datevalsDay))';
maxswedates=reshape(sweR.maxswedates,1,vecsize);

%one pixel at every timestep
parfor k=1:vecsize
    melt_vec=single(melt(:,k));
    sca_vec=fsca(:,k);
    swe_temp=zeros(length(datevalsDay),1);
    if maxswedates(k) > 0 && any(sca_vec)
        % look for contiguous sca pds
        % eliminate spurious zeros in SCA
        %         sca_vec=slidefun(@max,7,sca_vec,'central');
        [start,finish] = contiguous(sca_vec);
        %get rid of contiguous pds w/ start and finish prior to peak
        %whats left are contiguous pds that start before or on peak
        ind=start < maxswedates(k) & finish < maxswedates(k);
        start(ind)=[];
        finish(ind)=[];
        for h=1:length(start)
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
            
            FillValue=0;
            
            if isa(sweR.(fn{i}),'uint16')
                FillValue=intmax('uint16');
                sweR.(fn{i})(watermask)=FillValue;
            end
            
            h5create(rFile,location,sz,'Deflate',deflateLevel,'ChunkSize',...
                ChunkSize,'DataType',class(sweR.(fn{i})),'FillValue',FillValue);
            h5write(rFile,location,sweR.(fn{i}));
            %write mm for units for melt and swe
            if strcmpi(fn{i},'melt') || strcmpi(fn{i},'swe')
                h5writeatt(rFile,location,'units',units);
            end
        end
        h5writeProjection(rFile,'/Grid',target_hdr.ProjectionStructure);
        h5writeatt(rFile,'/','MATLABdates',datevalsDay);
        %added yyyymmdd attribute for non-MATLAB users
        h5writeatt(rFile,'/','ISOdates',...
            strjoin(cellstr(datestr(datevalsDay,'yyyymmdd'))));
        h5writeatt(rFile,'/Grid','ReferencingMatrix',target_hdr.RefMatrix);
end
t=toc;
fprintf('SWE reconstructed and saved in %3.1f min\n',t/60);