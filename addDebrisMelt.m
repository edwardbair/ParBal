function addDebrisMelt(Rname,DebrisMeltOutDir)
% add debris cover melt (under ice) to HDF5 reconstruction cubes
%input:
%Rname- recon h5 cube (can be remote location)
%DebrisMeltOutDir - where the daily debris melt mat files live

mf=0.2592; % mm W m^-2 day^-1

%copy to temp location (addresses h5 remote write issues)
[~,name,ext]=fileparts(Rname);
%local name
rname=fullfile(tempdir,[name ext]);
copyfile(Rname,rname);

%check if debris cover is already addded
info=h5info(rname,'/Grid');
skip=false;
k=1;
while ~skip && (k <= length(info.Datasets))
    skip=strcmp('DebrisMelt',info.Datasets(k).Name);
    k=k+1;
end
if skip
    fprintf('skipping %s because it already has a debris entry\n',rname);
else
    matdates=h5readatt(rname,'/','MATLABdates');
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
        G(G<0)=0; %ignore negative melt values for a day
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
%move from temp back to original location
movefile(rname,Rname);
end