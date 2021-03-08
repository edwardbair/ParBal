function [out,hdr]=mean_swe(m,d,yr,recondir)
%create cube of mean swe for a particular month, day, and vector of years
%m - month
%d - day
%yr - vector of years
%recondir - where the h5 cubes live
%output:
% out - h5 cube of that date for each year
% hdr - map header info
for i=1:length(yr)
    fname=sprintf('*%i.h5',yr(i));
    dd=dir(fullfile(recondir,fname));
%     fname=sprintf(...
%         'reconstruction_indus_CY%i.h5',yr(i));
if ~isempty(dd)
    fname=fullfile(recondir,dd.name);
else
    error('could not find %s',fullfile(recondir,fname));
end
    [swe,hdr]=getMelt(fname,'swe',...
        datenum([yr(i) m d]));
    if i==1
        out=zeros(...
            size(swe,1),size(swe,2),length(yr));
    end
    out(:,:,i)=swe;     
end
    
    