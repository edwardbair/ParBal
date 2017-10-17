% Correct NLDAS beam for elevation at ~100m
function BinZ = beamAtElevation(S0,tauZ,airmassZ)
%INPUT
% S0 - exoatmospheric direct solar irradiance
% muZ - cosine solar zenith angle (a grid!)
% tauZ - optical depth
% airmassZ - airmass
%OUTPUT
% BinZ - beam irradiance at elevation, normal to the sun
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2010-08-22 Karl Rittger rittger@nsidc.org 303-735-3433 
% National Snow & Ice Data Center
% Copyright 2010-2014 Regents of the University of California
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
% note that S0 has been scaled for rv, so it is not included in the tau
% calculation
BinZ = ( S0 .* exp(-tauZ.*airmassZ));

% Note that IPW requires normal to sun so older version used muIPW as input
%BinZ = (muZ .* S0 .* exp(-tauZ.*airmassZ))/muIPW;
