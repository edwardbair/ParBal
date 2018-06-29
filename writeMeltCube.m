function writeMeltCube(fsca_dir,melt_dir,out_dir,yr,nameprefix,varargin)
%write melt cube (no reconstruction) from daily inputs
%note: no canopy adjustment or watermask option
%input
%fsca_dir - where fsca h5 cubes live
%melt_dir - where melt daily mat files live
%yr - yrs for each cube
%hdr - map/geog projection info
%nameprefix -cube name prefix, yr will be addded
%optional:
%debris_dir - where debris melt files live
debris_flag=false;
debris_dir=[];
if ~isempty(varargin)
    debris_dir=varargin{1};
    debris_flag=true;
end
parfor i=1:length(yr)
   % matdates=datenum([yr(i) 1 1]):datenum([yr(i) 12 31]);
    d=dir(fullfile(fsca_dir,sprintf('*_%i.h5',yr(i))));
    [fsca,matdates,hdr]=GetEndmember(fullfile(d.folder,d.name),'snow');
    for j=1:length(matdates)
        ds=datestr(matdates(j),'yyyymmdd');
        d=dir(fullfile(melt_dir,[ds,'.mat']));
        m=matfile(fullfile(d.folder,d.name));
        if j==1
               SnowIceMelt=zeros([size(m.M,1),size(m.M,2),length(matdates)],'int16');
            if debris_flag
               DebrisMelt=zeros([size(m.M,1),size(m.M,2),length(matdates)],'int16');
            end
        end
        smelt=single(m.M).*0.0108.*fsca(:,:,j);
        SnowIceMelt(:,:,j)=int16(smelt);
        if debris_flag
            d=dir(fullfile(debris_dir,[ds,'.mat']));
            m2=matfile(fullfile(d.folder,d.name));
            Gmelt=sum(single(m2.G),3).*0.0108;
            %no negative melt and no melt w/ snow cover
            Gmelt(Gmelt<0 | fsca(:,:,j) > 0) = 0;
            DebrisMelt(:,:,j)=int16(Gmelt);
        end
        fprintf('done w/%s\n',ds);
    end

        fname=fullfile(out_dir,sprintf('%s%i.mat',nameprefix,yr(i)));
        m3=matfile(fname,'Writable',true);
        m3.SnowIceMelt=SnowIceMelt;
        m3.SnowIceMeltUnits='mm';
        if debris_flag
            m3.DebrisMelt=DebrisMelt;
            m3.DebrisMeltUnits='mm';
        end
        m3.hdr=hdr;
        m3.matdates=matdates;
end

    
    

