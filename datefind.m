function ind = datefind(sub,super,tol)
%DATEFIND Indices of date numbers in matrix.
%   IND = DATEFIND(SUB,SUPER,TOL) returns a vector of indices to the
%   date numbers in SUPER that are present in SUB, plus or minus the
%   tolerance TOL. SUPER is the superset matrix of non-repeating date
%   numbers whose elements are sought. SUB is a subset matrix of date
%   numbers used to find matching date numbers in SUPER. TOL is a
%   positive integer, representing the tolerance. The default is TOL = 0.
%
%   Notes:
%     If no date numbers match, IND = [].
%
%     The elements of SUB must be contained in SUPER, without
%     repetition. This function is designed to work with
%     non-repeating sequences of dates.
%
%
%   For example, given the following data:
%
%      Super = datenum(1997, 7, 1:31);
%      Sub = [datenum(1997, 7, 10); datenum(1997, 7, 20)];
%
%      ind = datefind(Sub,Super,1) returns
%
%      ind =  9
%            10
%            11
%            19
%            20
%            21
%
%   See also DATENUM.

%   Author(s): M. Reyes-Kattar, 06-04-97
%   Copyright 1995-2006 The MathWorks, Inc.
%   $Revision: 1.9.2.6 $   $Date: 2010/10/08 16:42:36 $

if nargin < 3
   tol = 0;
end

if nargin < 2
   error(message('finance:datefind:missingInputs'));
end

sub = sub(:);
super = super(:);

% Find the occurrences of sub in super and take into account tolerance
ind = [];
for idx = 1:length(sub)
   d = abs(super - sub(idx));
   i = find(d <= tol);

   ind = [ind; i]; %#ok
end

ind = sort(ind);


% [EOF]
