function [ldasInterp,ceresInterp]=makeInterp(ldas_filelist,ldas_topo,topo,mask,ceres,tz,...
    varargin)
%make interpolated gldas and ceres data
%input:
%ldas_filelist - list of ldas files to load
%ldas_topo - ldas topo struct
%topo - finescale topo struct
%mask - mask for place to skip interpolation, e.g. from fsca or debris
%ceres - struct w/ ceres filelist
%tz - timezone
%LDASOnlyFlag - use LDAS only inputs, optional
%ouput
% LDASInterp - interp. LDAS struc
% ceresInterp - interp CERES struc (empty if LDASOnlyFlag)
% reproject and interpolate ldas for specified datevals
LDASOnlyFlag=false;
gldas_flag=false;
if contains(ldas_filelist.filenames{1},'GLDAS')
    gldas_flag=true;
end
if nargin==7
    LDASOnlyFlag=varargin{1};
end
ldas_stacked = stackGLDAS(ldas_filelist);
ldas_subset = subsetGLDAS(ldas_stacked,ldas_topo.CoarseRefMatrix,...
    topo.hdr,ldas_filelist.var,mask,true);
clear ldas_stacked
%only do hourly interpolation for gldas
if gldas_flag
    ldasInterp=interpGLDAS2hourly(ldas_subset,ldas_filelist.var);
else
    ldasInterp=ldas_subset;
end
clear ldas_subset
%for n/gldas, we now have local times 00:00-23:00 in 1 hr increments
% update 2018-6-28, now using ceres 4A 1 hr,so CERES is now hourly, but
% offset by 30 min from GLDAS.

% old*** do the same for ceres data, which is offset from GLDAS by +1.5 hr,
%and therefore needs the 2230 local time from the previous day to the
%130 local time in order interpolate from midnight to midnight ***old
if ~LDASOnlyFlag
    %     getCERESdatevalsUTC=(gldas_stacked.datevalsUTC(1)-...
    %     0.5/24:1/24:gldas_stacked.datevalsUTC(end)+0.5/24);
    %     equates to local vals of 1130 on day-1 to 0:30 on day +1, every hour
    getCERESdatevalsUTC=(ldasInterp.datevalsUTC(1)-0.5/24):1/24:...
        (ldasInterp.datevalsUTC(end)+1.5/24);
    [ceres_stacked,ceres_hdr,ceres_varnames]=read_ceres(ceres.ceres_dir,...
        getCERESdatevalsUTC,tz,ceres.var);
    ceres_subset = subsetGLDAS(ceres_stacked,ceres_hdr.RefMatrix,...
        topo.hdr,ceres_varnames,mask,true);
    clear ceres_stacked
    %still need to interpolate to shift to 00:00-23:00
    ceresInterp = interpGLDAS2hourly(ceres_subset,ceres_varnames);
    clear ceres_subset
else
    ceresInterp=[];
end
end