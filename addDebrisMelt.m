function addDebrisMelt(ReconCubeOutDir,DebrisMeltOutDir)
% add debris cover melt (under ice) to HDF5 reconstruction cubes
%input:
%ReconCubeOutDir - where the recon cubes live
%DebrisMeltOutDir - where the daily debris melt mat files live
d=dir(fullfile(ReconCubeOutDir,'*.h5'));
mf=0.2592; % mm W m^-2 day^-1
parfor i=1:length(d)
    rname=fullfile(d(i).folder,d(i).name);
    matdates=h5readatt(rname,'/','MATLABdates')
    fprintf('reading debris melt for %s\n',rname);
    for j=1:length(matdates)
        dname=fullfile(DebrisMeltOutDir,sprintf('%s.mat',datestr(matdates(j),'yyyymmdd')));
        M=matfile(dname);
        if j==1 %preallocate
            DebrisMelt=zeros([size(M.G,1) size(M.G,2) length(matdates)],'single');
        end
        G=M.G;
        C=class(G);
        G=single(G);
        FillValue=intmax(C);
        G(G==FillValue)=NaN;
        DebrisMelt(:,:,j)=G*mf;
    end
    fprintf('finished reading debris melt for %s\n',rname);
    %output variables
    h5create(rname,'/Grid/DebrisMelt',size(DebrisMelt),...
        'Deflate',9,...
        'ChunkSize',[size(DebrisMelt,1) size(DebrisMelt,2) 1],...
        'DataType',C,...
        'FillValue',FillValue);
    DebrisMelt(isnan(DebrisMelt))=FillValue;
    DebrisMelt=cast(DebrisMelt,C);
    h5write(rname,'/Grid/DebrisMelt',DebrisMelt);
    h5writeatt(rname,'/Grid/DebrisMelt','units','mm');
    fprintf('finished writing debris melt for %s\n',rname);
end
end