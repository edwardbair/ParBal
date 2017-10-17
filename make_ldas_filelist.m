function ldas=make_ldas_filelist(datevalsDay,ldasdir,utc_offset,varargin)
% creates ldas filelist
% input
% datevalsDay, localtime, matdates for dates in LDAS
% ldasdir, string ldas directory, e.g 'C:\raid\scratch\reconstruction\NLDAS\';
% utc offset , scalar in fraction of a day, e.g. -8/24 = PST
% optional:
% nldas secondary forcing dir, which makes an additional field where
% secondary forcing data is located
% output
% ldas filelist structure

if isempty(varargin);
    secondary_flag=false;
else
    secondary_flag=true;
    secondary_dir=varargin{1};
end

if ~isempty(regexpi(ldasdir,'NLDAS'));
    ldastype='NLDAS';
    timestep=1/24; % hourly
elseif ~isempty(regexpi(ldasdir,'GLDAS'));
    ldastype='GLDAS';
    timestep=3/24; % 3 hr
else
    error('could not determine which LDAS (NLDAS or GLDAS) you are using')
end

ldas=struct('datevalsUTC',[],'datevalsLocal',[],'filenames',[]);

%utc dates
datevalsDayUTC=datevalsDay-utc_offset;

for i=1:length(datevalsDay);
    %datevals within day i, local
    daily_datevals=datevalsDay(i):timestep:datevalsDay(i)+1;%-timestep;
    %datevals within day i, utc
    daily_datevalsUTC=daily_datevals-utc_offset;
    %datevals within day that conform to ldas sampling (eg. 00,03,06,...)
    tmp1=floor(datevalsDayUTC(i)):timestep:floor(datevalsDayUTC(i))+1;%-timestep;
    tmp2 = abs(daily_datevalsUTC(1)-tmp1);
    [~,idx] = min(tmp2); %index of closest value
    daily_datevalsUTCldasInterval=tmp1(idx):timestep:tmp1(idx)+1;%-timestep;
    switch ldastype
        case 'NLDAS'
            %need 00:00 to 23:00
            daily_hour_length=length(daily_datevals)-1;
        case 'GLDAS'
            %need 00:00 to 00:00 (+1), include overlap with next day for
            %interpolation
            daily_hour_length=length(daily_datevals);
    end
    for j=1:daily_hour_length
        %set up all datevals to conform to available gldas intervals
        dv_utc=datevec(daily_datevalsUTCldasInterval(j));
        ldas(i).datevalsUTC(j)=daily_datevalsUTC(j);
        ldas(i).datevalsLocal(j)=daily_datevalsUTC(j)+utc_offset;
        doy=datenum(dv_utc(1:3))-datenum([dv_utc(1) 1 1])+1;
        subdir=fullfile(ldasdir,num2str(dv_utc(1),'%04i'),...
            num2str(doy,'%03i'));
        switch ldastype
            case 'NLDAS'
                fname=fullfile(subdir,sprintf('*%s*.grb',datestr(dv_utc,...
                'yyyymmdd.HHMM')));
            case 'GLDAS'
                 fname=fullfile(subdir,sprintf('*%s%03i%s*.grb',...
                datestr(daily_datevalsUTCldasInterval(j),'yyyy'),doy,...
                datestr(daily_datevalsUTCldasInterval(j),'.HHMM')));
        end
        f=dir(fname);
        if ~isempty(f)
            ldas(i).filenames{j}=fullfile(subdir,f.name);
        else
            error('%s\n does not exist\n',fname);
        end
        if secondary_flag
           subdir=fullfile(secondary_dir,num2str(dv_utc(1),'%04i'),num2str(doy,'%03i'));
           fname=fullfile(subdir,sprintf('*%s*.grb',datestr(dv_utc,'yyyymmdd.HHMM')));
           f=dir(fname);
           if ~isempty(f)
           ldas(i).secondary_filenames{j}=fullfile(subdir,f.name);
           else
           error('%s\n does not exist\n',fname);
           end
        end 
    end
end