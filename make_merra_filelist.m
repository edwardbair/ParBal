function ldas=make_merra_filelist(datevalsDay,merradir,utc_offset)
% creates ldas filelist
% input
% datevalsDay, localtime, matdates for dates in MERRA
% merradir, string ldas directory
% utc offset , scalar in fraction of a day, e.g. -8/24 = PST

% output
% merra filelist structure

timestep=1/24; % hourly
merra=struct('datevalsUTC',[],'datevalsLocal',[],'filenames',[]);

%utc dates
datevalsDayUTC=datevalsDay-utc_offset;

for i=1:length(datevalsDay)
    %datevals within day i, local
    daily_datevals=datevalsDay(i):timestep:datevalsDay(i)+1;
    %datevals within day i, utc
    daily_datevalsUTC=daily_datevals-utc_offset;
    udaily_datevalsUTC=unique(floor(daily_datevalsUTC));
      
    for j=1:udaily_datevalsUTC
        %set up all datevals to conform to available merra intervals
        dv_utc=datevec(daily_datevalsUTCmerraInterval(j));
        merra(i).datevalsUTC(j)=daily_datevalsUTC(j);
        merra(i).datevalsLocal(j)=daily_datevalsUTC(j)+utc_offset;
        
        fname=fullfile(merradir,sprintf('*%s*.nc',datestr(dv_utc,...
        'yyyymmdd')));
           
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