function [out,MATLABdates]=read_energy(filename,varname,row,col,varargin)
% input: filename - mat filename for daily recon outputs,
% row,col - pixel row/col to load
% optional, hr, 1-24 to load; otheriwse all 24 hr loaded
% varname, one of these
% M
% Tsfc 
% LinZ 
% Lout
% Ta 
% diffuseZ 
% directZ 
% albedo 
% presZ
% windspd
% output: var transformed to single precision original
% values
    m=matfile(filename);

if isempty(varargin)
    out=squeeze(m.(varname)(row,col,:));
%     MATLABdates=m.MATLABdates;
    MATLABdates=datenum(filename(end-11:end-4),'yyyymmdd');
    MATLABdates=MATLABdates:1/24:MATLABdates+23/24;
else
    ind=varargin{1};
    out=squeeze(m.(varname)(row,col,ind));
    out=out(row,col);
    MATLABdates=datenum(filename(end-11:end-4),'yyyymmdd')+ind/24;
%     MATLABdates=m.MATLABdates(1,ind);
end
% nullval=intmin(class(out));
fn=fieldnames(m);
out=single(out);
if any(strcmp('divisors',fn))
divisors=m.divisors;
divisor=table2array(divisors(varname,1));
out=out./divisor;
end
% out(out==nullval)=NaN;