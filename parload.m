function varargout = parload(fname,varargin)
%varargout = parload(fname,varargin)
% Loads variables (varargin) from MATLAB file "fname"
%INPUT
%fname - filename with variables to be loaded
%varargin - cell array (list) of variables to be loaded
%OUTPUT
%varargout - output variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2011-03-18 Karl Rittger rittger@nsidc.org 303-735-3433 
% National Snow & Ice Data Center
% Copyright 2011-2011 Regents of the University of California
%
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
%
% http://www.apache.org/licenses/LICENSE-2.0
%
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Number of input and optional arguements
o=size(varargin,2);
%disp([' Number of optional input arguements = ' num2str(o) ' to parload'])
for j=1:o
    %disp(['Loading ' varargin{j}])
    varargout{j} = cell2mat(struct2cell(load(fname,varargin{j})));
end