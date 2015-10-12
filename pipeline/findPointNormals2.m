function [ normals, curvature ] = findPointNormals2(points, point, viewPoint, dirLargest)
%FINDPOINTNORMALS Estimates the normals of a sparse set of n 3d points by
% using a set of the closest neighbours to approximate a plane.
%
%   Required Inputs:
%   points- nx3 set of 3d points (x,y,z)
%
%   Optional Inputs: (will give default values on empty array [])
%   numNeighbours- number of neighbouring points to use in plane fitting
%       (default 9)
%   viewPoint- location all normals will point towards (default [0,0,0])
%   dirLargest- use only the largest component of the normal in determining
%       its direction wrt the viewPoint (generally provides a more stable
%       estimation of planes near the viewPoint, default true)
%
%   Outputs:
%   normals- nx3 set of normals (nx,ny,nz)
%   curvature- nx1 set giving the curvature
%
%   References-
%   The implementation closely follows the method given at
%   http://pointclouds.org/documentation/tutorials/normal_estimation.php
%   This code was used in generating the results for the journal paper
%   Multi-modal sensor calibration using a gradient orientation measure 
%   http://www.zjtaylor.com/welcome/download_pdf?pdf=JFR2013.pdf
%
%   This code was written by Zachary Taylor
%   zacharyjeremytaylor@gmail.com
%   http://www.zjtaylor.com

%% check inputs
validateattributes(points, {'numeric'},{'ncols',3});

if(nargin < 3)
    viewPoint = [];
end
if(isempty(viewPoint))
    viewPoint = [0,0,0];
else
    validateattributes(viewPoint, {'numeric'},{'size',[1,3]});
end

if(nargin < 4)
    dirLargest = [];
end
if(isempty(dirLargest))
    dirLargest = true;
else
    validateattributes(dirLargest, {'logical'},{'scalar'});
end

%% setup

%ensure inputs of correct type
points = double(points);
viewPoint = double(viewPoint);

%remove self
points = points(2:end,:);
nNeigh = size(points,1);
%find difference in position from neighbouring points
p = repmat(point,nNeigh,1) - points;
%p = reshape(p, size(point,1),numNeighbours,3);

%calculate values for covariance matrix
C = zeros(size(point,1),6);
C(1) = sum(p(:,1).*p(:,1),1);
C(2) = sum(p(:,1).*p(:,2),1);
C(3) = sum(p(:,1).*p(:,3),1);
C(4) = sum(p(:,2).*p(:,2),1);
C(5) = sum(p(:,2).*p(:,3),1);
C(6) = sum(p(:,3).*p(:,3),1);
C = C ./ nNeigh;

%% normals and curvature calculation

normals = zeros(size(point));
%curvature = zeros(size(point,1),1);
%form covariance matrix
    Cmat = [C(1) C(2) C(3);...
        C(2) C(4) C(5);...
        C(3) C(5) C(6)];  
    
    %get eigen values and vectors
    [v,d] = eig(Cmat);
    d = diag(d);
    [lambda,k] = min(d);
    
    %store normals
    normals(:) = v(:,k)';
    
    %store curvature
    curvature = lambda / sum(d);

%% flipping normals

%ensure normals point towards viewPoint
point = point - repmat(viewPoint,size(point,1),1);
if(dirLargest)
    [~,idx] = max(abs(normals),[],2);
    idx = (1:size(normals,1))' + (idx-1)*size(normals,1);
    dir = normals(idx).*point(idx) > 0;
else
    dir = sum(normals.*point,2) > 0;
end

normals(dir,:) = -normals(dir,:);

%%%% My line %%%% not sure
%curvature(dir) = - curvature(dir);
%%%%%%%%%%

end