function [Z,cc,hdr]=mosaicMCD12Q1(ccfile,mcd12dir,tilesv,tilesh)
%mosaic and reproject MCD12Q1 into landcover structure
%input
%ccfile - canopy cover file w/ header
% eg '/raid/sandbox/snowhydro/nbair/datasets/Landcover/Indus_cc_463m_sinusoidal.mat'
%mcd12dir - location of MCD12 files
% eg     '/raid/sandbox/snowhydro/scratch/MODIS/mcd12q1'
%tilesv - vertical tile(s) eg {'v05'};
%tilesh - horizontal tiles(s), eg tilesh={'h23','h24','h25'};

CC=load(ccfile);
cc=CC.cc;
hdr=CC.hdr;

Z=false(2400*length(tilesv),2400*length(tilesh));

%1-deciduous
%2-coniferous
% 
% tiles={'h23v05','h24v05','h25v05'};
% tilesv={'v05'};
% tilesh={'h23','h24','h25'};

% pxy=[1:2400];
% pxx=[1:2400; 2401:4800; 4801:7200];

for i=1:2 %two veg types
    for j=1:length(tilesv)
        for k=1:length(tilesh)
            if i==1 && j==1 && k==1 %get refMatrix
                [R,mstruct]=sinusoidProjMODtile([tilesh{k} tilesv{j}]);
            end
            fname=['MCD12Q1*' tilesh{k} tilesv{j} '*.hdf'];
            d=dir(fullfile(mcd12dir,fname));
            x=hdfread(fullfile(mcd12dir,d.name),'LC_Type1');
            if i==1 % deciduous (IGBP class
                %3 - decid needleleaf forest > 60% tree cover
                %4 - decid broadleaf forests > 60% tree cover
                x=x==3 | x==4;
            elseif i==2 %coniferous
                % 1 - evergreen needleleaf forests > 60% tree cover
                % 2 - evergreen broadleaf forests > 60% tree cover
                % 8 - woody savanas, 30-60 % tree cover (canopy > 2m)
                % 9 - savannas, 10-30% tree cover (canopy > 2m) 
                % 5 - mixed decid. evergreen , tree cover > 60% (canopy > 2
                % m)
                x=x==1 | x==2 | x==8 | x==9 | x==5;
            end
            pxy=(1:2400)+2400*(j-1);
            pxx=(1:2400)+2400*(k-1);
            
            Z(pxy,pxx,i)=x;
        end
    end
end

Z=reprojectRaster(Z,R.RefMatrix_500m,mstruct,CC.hdr.ProjectionStructure,...
    'rasterref',CC.hdr.RasterReference,'method','nearest');