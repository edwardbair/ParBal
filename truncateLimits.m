function [ Y ] = truncateLimits( X,minvalue,varargin )
%truncate X to stay between limits
% Y = truncateLimits( X,limits )
% Y = truncateLimits(X,minvalue,maxvalue)
%
% Input
%   X - input data, vector or matrix or N-D
%   To specifiy limits, either
%       limits - vector of length 2: min and max of output
%       or minvalue,maxvalue as separate arguments
%
% Output
%   Y - truncated X with min and max of limits
%       (NaNs in input are NaN in output)

assert(isfloat(X),'input X must be floating point')
numargs = 2;
if nargin==numargs
    limits = minvalue;
    assert(length(limits)==2,...
        'if just one minmax argument specified, must be vector of length 2')
else
    limits = [minvalue varargin{1}];
end
Y = X;
tmin = X<min(limits);
tmax = X>max(limits);
Y(tmin) = min(limits);
Y(tmax) = max(limits);
Y(isnan(X)) = NaN; % just to make sure

end