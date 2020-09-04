function gldasInterp=interpGLDAS2hourly(gldas,vars)
% Interpolate (3hr) GLDAS data to 1 hr
% INPUTS: 
% gldas - stucture of gldas data from subset_GLDAS
% vars - list of GLDAS variables you want to interpolate
% OUTPUT:
% gldasInterp - hourly interpolated gldas data for one day

gldasInterp=gldas;
% 24 hr of data
%compute timezone
tz=gldas.datevalsLocal(1)-gldas.datevalsUTC(1);
%note the use of index 2 for cases to avoid starting from the previous day
%for CERES data
gldasInterp.datevalsLocal=floor(gldas.datevalsLocal(2)):1/24:...
    floor(gldas.datevalsLocal(2))+1-1/24;
gldasInterp.datevalsUTC=gldasInterp.datevalsLocal-tz;
sg=size(gldas.(vars{1}));
sizea=[sg(1) sg(2) length(gldasInterp.datevalsUTC)];

[X,Y,Z]=meshgrid(1:sizea(2),1:sizea(1),gldas.datevalsUTC);
[xq,yq,zq]=meshgrid(1:sizea(2),1:sizea(1),gldasInterp.datevalsUTC);
    
for i=1:length(vars)
    gldasInterp.(vars{i})= interp3(X,Y,Z,gldas.(vars{i}),xq,yq,zq);
end