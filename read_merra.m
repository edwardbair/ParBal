function [C,hdr,varnames]=read_merra(merra_dir,datevalsUTC,tz,var)
% function to read merra to go with ldas data
% input
% merra_dir, where the merra nc files live
% datevalsUTC - matdates to retrieve

% var, variables names to retrieve, as a cell

% output
% C, struct containing merra info
% hdr, coordinate info
% varnames, names of variables specified by numbers

datetimes=dateshift(datetime(datevalsUTC,'ConvertFrom','datenum'),...
    'start','minute','nearest');
offset=0.5/24; %times start at 00:30

d=dir(fullfile(merra_dir,'*.nc'));

varnames=var;
hdr.gridtype='geographic';

umatdates=unique(floor(datevalsUTC));

for h=1:length(umatdates)
    match=false;
    i=1;

    while ~match && i <=length(d)
        fname=fullfile(merra_dir,d(i).name);
        [~,name]=fileparts(fname);
        matdateFromFile=datenum(name(end-11:end-4),'yyyymmdd');
        if umatdates(h) == matdateFromFile
            match=true; % found a file, for a given day
           
            timerange=double(ncread(fname,'time'));
            timerange=dateshift(datetime(matdateFromFile+offset+timerange/1440,...
                'ConvertFrom','datenum'),'start','minute','nearest');
            [~,id_dt,id_tr] = intersect(datetimes,timerange);
            for j=1:length(var)

                x=ncread(fname,var{j});
                mv=ncreadatt(fname,var{j},'missing_value');
                x(x==mv)=NaN;
                x=single(x);
                x=pagetranspose(x);
                x=flipud(x);

                if h==1 % allocate for first run
                    sz=size(x);
                    C.(var{j})=zeros([sz(1) sz(2) length(datetimes)]);
                    C.datevalsUTC=datevalsUTC;
                    C.datevalsLocal=datevalsUTC+tz;
                    lon=ncread(fname,'lon');
                    lat=ncread(fname,'lat');
                    hdr.RefMatrix=makerefmat(lon(1),lat(end),lon(2)-lon(1),...
                        -(lat(end)-lat(end-1)));
%                     hdr.RasterReference=refmatToGeoRasterReference(hdr.RefMatrix,[sz(2)...
%                         sz(1)]);
                end
                C.(var{j})(:,:,id_dt)=x(:,:,id_tr);
            end
        end
    i=i+1;
    end
end