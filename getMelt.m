function [x,hdr,h5mdates]=getMelt(h5file,meltvariable,varargin)
%reads different melt varibles from h5 files output by
%reconstructSWE
%input:
%h5file - HDF 5 file, string
%meltvariable - variable to read, string, choices are:
% 'swe' - daily reconstructed swe, mm
% 'melt' - daily melt, mm
% 'maxswedates' - date of max swe
% optional 3rd input is scalar or vector of matdates to read, if not
% supplied, whole cube is read
%output:
% x - image or cube of requested variable, corresponding to dates
% requested
% hdr - geographic info
% h5mdates - matdates

if nargin==2
    dflag=false;
elseif nargin==3
    dflag=true;
    matdates=varargin{1};
end

varloc=['/Grid/' meltvariable];
i=h5info(h5file,varloc);
sz=i.Dataspace.Size;
hdr=GetCoordinateInfo(h5file,'/Grid',sz(1:2));
h5mdates=h5readatt(h5file,'/','MATLABdates');

if ~dflag
    x=h5read(h5file,varloc);
else
    
    if length(matdates)==1 %1 day
       n=find(matdates==h5mdates);
       x=h5read(h5file,varloc,[1 1 n],[sz(1) sz(2) 1]);
    else
       x=zeros([sz(1) sz(2) length(matdates)]);
        for i=1:length(matdates) %2 or more days
            n=find(matdates(i)==h5mdates);
            x(:,:,i)=h5read(h5file,varloc,[1 1 n],[sz(1) sz(2) 1]);
        end
    end     
end

