function add_ice_mask(fscadir,rdir,wy)

% add an ice mask to reconstructed swe
% input:
% fscadir, fsca directory, e.g. 'C:\raid\scratch\nbair\modscag\h23v05';
% rdir, reconstructed swe directory, e.g. 'C:\raid\data\nbair\datasets\reconstructions\h23v05';
% wy vector, e.g. 2000:2015;
for i=1:length(wy);
    d=dir(fullfile(fscadir,sprintf('*WY%i.h5',wy(i))));
    fsca=GetEndmember(fullfile(fscadir,d.name),'snow');
    ice=min(fsca,[],3);
    icemask=false(size(ice));
    icemask(ice>0)=true;
    d2=dir(fullfile(rdir,sprintf('*WY%i.h5',wy(i))));
    h5name=fullfile(rdir,d2.name);
    h5create(h5name,'/Grid/ice',size(icemask),'DataType','uint8',...
        'ChunkSize',size(icemask),'Deflate',9);
    h5write(h5name,'/Grid/ice',uint8(icemask));
end
end