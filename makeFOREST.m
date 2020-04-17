function FOREST = makeFOREST(LandCover)

% Link and Marks (1999) 
% Distributed simulation of snowcover mass- 
% and energy-balance in the boreal forest
% Hydrological Processes 2439-2452
%INPUT
% FOREST LandCover structure called 
% fields:
% Z - m x n x L, with each layer of L corresponding to the 
% fractional cover of type L, 0-1
% currently only two veg types are implemented:
% 1 - deciduous
% 2 - coniferous
% cc - m x n, static canopy cover estimate, 0-1

%OUTPUT
% FOREST structure
%       .tau - canopy transmissivity coefficient
%       .u - canopy extinction coefficent
%       .cc - canopy cover fraction
%       .h - canopy height, m
%       .shd - snow holding depth, m
%       .type - canopy type structure w/:
%       .type.list - {1=deciduous,2=coniferous}
%       .type.num_val - map of dominant type 

% Allocate the transmissivity, extinction, 
%fraction canopy and height grids
siz=size(LandCover.Z);
FOREST.tau=ones([siz(1) siz(2)],'single');
FOREST.u=zeros([siz(1) siz(2)],'single');
FOREST.cc=FOREST.u;
FOREST.h=FOREST.u;

deciduous_fraction=squeeze(LandCover.Z(:,:,1));
coniferous_fraction=squeeze(LandCover.Z(:,:,2));
% Set values for transmissivity
FOREST.tau(deciduous_fraction > 0 & ...
    deciduous_fraction >= coniferous_fraction)=0.60;%Deciduous
FOREST.tau(coniferous_fraction > 0 & ...
    coniferous_fraction >= deciduous_fraction)=0.30;%Coniferous

%classification map
FOREST.type.list={'deciduous','coniferous'};
FOREST.type.num_val=zeros(size(FOREST.tau));
FOREST.type.num_val(FOREST.tau==0.60)=1;
FOREST.type.num_val(FOREST.tau==0.30)=2;

%vegetation snow holding depth (for winds)
FOREST.shd=zeros(size(FOREST.tau));
FOREST.shd(FOREST.type.num_val==1)=12;
FOREST.shd(FOREST.type.num_val==2)=15;
% Set values for forest height
% assume veg snow holding depth = height;
FOREST.h=FOREST.shd;
%bare ground
FOREST.shd(FOREST.type.num_val~=1 & FOREST.type.num_val~=2)=0.01; 

% Set values for extinction
FOREST.u(deciduous_fraction > 0 & ...
    deciduous_fraction >= coniferous_fraction)=0.016;%Deciduous
FOREST.u(coniferous_fraction > 0 & ...
    coniferous_fraction >= deciduous_fraction)=0.033;%Coniferous

FOREST.cc=LandCover.cc;
FOREST.tau(FOREST.cc==0)=1;
FOREST.u(FOREST.cc==0)=0;