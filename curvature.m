function c=curvature(dem,RasterReference,curve_len_scale)
% input dem and associated rasterrefereence
% curvature length scale, m, 1/2 wavelength of topographic features
% output scaled (-0.5:0.5) curvature, size of dem
% from Glen Liston's Snowmodel - Micromet Fortran code
% Liston et al. (2007) Simulating complex snow distrubutions 
% in windy envirnoments.
% Journal of Glaciology
 c=zeros(size(dem));
%  curve_len_scale=600; %m
 
 rows=RasterReference.RasterSize(1);
 cols=RasterReference.RasterSize(2);
 deltaxy = mean([RasterReference.CellExtentInWorldX,...
     RasterReference.CellExtentInWorldY]);
 inc = max(1,round(curve_len_scale/deltaxy));
 
%  nx=rows; snowtran code is backwards, 
%  ny=cols; 
 %compute curvature
 for i=1:rows
     for j=1:cols
     c(i,j) = (4.0 * dem(i,j) -...
     dem(max(1,i-inc),max(1,j-inc)) - ...
     dem(min(rows,i+inc),min(cols,j+inc)) - ...
     dem(min(rows,i+inc),max(1,j-inc)) - ...
     dem(max(1,i-inc),min(cols,j+inc))) / ...
     (sqrt(2.0) * 16.0 * real(inc) * deltaxy) + ...
     (4.0 * dem(i,j) - ...
     dem(min(rows,i+inc),j) - dem(max(1,i-inc),j) - ...
     dem(i,min(cols,j+inc)) - dem(i,max(1,j-inc))) / ...
     (16.0 * real(inc) * deltaxy);
     end
 end
 
% Scale the curvature such that the max abs(curvature) has a value
% of abs(0.5).  Include a 1 mm curvature in curve_max to prevent
% divisions by zero in flat terrain where the curvature is zero.
curve_max = 0.0 + 0.001;
curve_max = max(curve_max,max(abs(c(:))));
c=c./(2.*curve_max);
 
 