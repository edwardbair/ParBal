function daily_melt(dateval,gldas_filelist,ceres,topo,gldas_topo,ceres_topo,...
    tz,FOREST,sFile,outfile,fast_flag,LDASOnlyFlag,metvars_flag)
% loads,stacks,subsets,and interpolates GLDAS data for a single day
%then calls daily_energy and saves energy outputs
%INPUT
% dateval -  dateval for single day
% gldas_filelist - struct of gldas filenames
% ceres - struct of ceres filenames and other info
% topo - fine scale topo struct (DEM,SlopeAspect,ViewFactor,Horizons)
% gldas_topo - coarse gldas topo struct
% ceres_topo - coarse ceres topo struct
% tz - timezone for target area
% FOREST - forest data structure
% sFile - fsca h5 file
% outfile - output file
% fast flag - true - only solve for M; false - solve for all outputs; 
% only set for 'normal'
% LDASOnlyFlag - run w/ only LDAS inputs
% metvars_flag - true, save met vars; false - dont save met vars

[fsca,~,hdr]=GetEndmember(sFile,'snow',dateval);
%if the topo hdr doesn't match the fsca hdr, reproject
if ~isequal(topo.hdr,hdr)
    fsca=reprojectRaster(fsca,hdr.RefMatrix,hdr.ProjectionStructure,...
        topo.hdr.ProjectionStructure,'rasterref',topo.hdr.RasterReference);
end

dateS=datestr(dateval);
mask=fsca>0;
t1=clock;

[gldasInterp,ceresInterp]=makeInterp(gldas_filelist,...
        gldas_topo,topo,mask,ceres,tz,LDASOnlyFlag);
t2=clock;
disp(['-->> prepped forcings in ' num2str(round(etime(t2,t1))) ...
    ' seconds for ',dateS])
t1=clock;
dailyEnergy(topo,gldasInterp,gldas_topo,ceresInterp,...
        ceres_topo,fast_flag,metvars_flag,'normal',outfile,FOREST,sFile);
t2=clock;
disp(['dailyEnergy in ' num2str(roundn(etime(t2,t1)/60,-1)) ' minutes for ',dateS])