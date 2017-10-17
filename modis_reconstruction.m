function modis_reconstruction(FOREST, topo, ldas, ldas_topo, ...
    ldas_vars, ldas_vars_locs, datevalsDay, eFiledir, sFile, ...
    reconfile, maxswedates, poolsize, mask)

t1=clock;
parpool_check(poolsize);
parfor d=1:length(datevalsDay)
% for d=73;
    ldas_filelist = ldas(d);
    ldas_filelist.vars=ldas_vars;
    ldas_filelist.vars_locs=ldas_vars_locs;
    energy_modis(datevalsDay(d),topo,FOREST,ldas_filelist,ldas_topo,...
        sFile,eFiledir);
end
t2=clock;
disp(['energy downscaled in ' num2str(roundn(etime(t2,t1)/60,-1)) ' minutes']);

t1=clock;
melt = reconstructSWE_modis(datevalsDay,eFileDir,sFile,true);
% melt = reconstructSWE_modis(datevalsDay,eFiledir,sFile,method,...
%     sw_bias_correct,veg_adj);
t2=clock;
disp(['SWE reconstructed in ' num2str(roundn(etime(t2,t1)/60,-1)) ' minutes']);

t1=clock;
make_recon_cube_modis(topo.hdr,melt,datevalsDay,reconfile,sFile,...
    maxswedates,veg_adj,mask)
t2=clock;
disp(['SWE cube saved in ' num2str(roundn(etime(t2,t1)/60,-1)) ' minutes']);