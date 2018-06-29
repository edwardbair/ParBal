function [S] = SnowCloudLimits()
%limits of snow and cloud properties used in various functions
%
%Output
%   S - structure of 2-element vectors with limits for the following
%       variables
%   snowRadius, mum
%   iceCloudRadius, mum
%   waterCloudRadius, mum
%   iceCloudWE, mm
%   waterCloudWE, mm
%   dust (in snow) mass fraction
%   soot (in snow) mass 
%   dust radius
%   soot radius

S.snowRadius = [30 1500];
S.waterInSnowRadius = S.snowRadius;
S.iceCloudRadius = [5 50];
S.waterCloudRadius = [1 25];
S.dustRadius = [.1 10];
S.sootRadius = [0.01 5];
S.defaultDustRadius = 1;
S.defaultSootRadius = .1;
S.iceCloudWE = [.01 3];
S.waterCloudWE = [.01 10];
S.dust = [0 1000]/1e6;
S.soot = [0 5000]/1e9;
S.wetSnow = [0 .15];
S.unitsSize = 'mum';
S.unitsWE = 'mm';
end