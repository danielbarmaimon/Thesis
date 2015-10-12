function [ noisyPoints ] = addNoise( points, a)
%   Add Gaussian noise to a banch of 3D points, given the points and sigma
%   value.
%   Inputs:
%           a : Sigma value. It is the same for x, y and z dimensions
%           points : Original 3D points over which the noise is going to be added 

% Defining the Gaussian distribution 
MU = [0 0 0];
SIGMA = a*eye(3);
obj = gmdistribution(MU,SIGMA);

% Generating random points that belongs to the distribution
rng(1); % For reproducibility (the same results would be generated with this line)
numPoints = size(points,1);
noise = random(obj,numPoints);

% Adding the noise to the points
noisyPoints = points + noise;
end

