function [gldasInterp,ceresInterp]=makeInterp(gldas_filelist,gldas_topo,topo,mask,ceres,tz)
%make interpolated gldas and ceres data
%input:
%gldas_filelist - list of gldas files to load
%gldas_topo - gldas topo struct
%topo - finescale topo struct
%mask - mask for place to skip interpolation, e.g. from fsca or debris
%ceres - struct w/ ceres filelist
% reproject and interpolate ldas for specified datevals
    gldas_stacked = stackGLDAS(gldas_filelist);
    gldas_subset = subsetGLDAS(gldas_stacked,gldas_topo.CoarseRefMatrix,...
        topo.hdr,gldas_filelist.vars,mask,true);
    gldasInterp=interpGLDAS2hourly(gldas_subset,gldas_filelist.vars);
    % do the same for ceres data, which is offset from GLDAS by +1.5 hr,
    %and therefore needs the 2230 local time from the previous day to the
    %130 local time in order interpolate from midnight to midnight
getCERESdatevalsUTC=(gldas_stacked.datevalsUTC(1)-...
    1.5/24:3/24:gldas_stacked.datevalsUTC(end)+1.5/24);
     [ceres_stacked,ceres_hdr,ceres_varnames]=read_ceres(ceres.ceres_dir,...
        getCERESdatevalsUTC,tz,ceres.var_num);
    ceres_subset = subsetGLDAS(ceres_stacked,ceres_hdr.RefMatrix,...
        topo.hdr,ceres_varnames,mask,true);
    ceresInterp = interpGLDAS2hourly(ceres_subset,ceres_varnames);
end