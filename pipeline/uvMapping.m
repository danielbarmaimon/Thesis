function [ intensity3Dpoints ] = uvMapping( points3Dface, uvMap, img)
%   Returns the values of the intensity (gray scale), given: 
%       points3Dface =  vertices, 
%       uvMap =         mapping over the 2D image
%       img =           2D gray image to get the intensity
[height, width] = size(img);
npixels = size(points3Dface,1);
a = repmat([height width], npixels, 1); % Repetition of the size of the 2D image
uvMapPoints2D = uvMap.*a;                % Conversion of the normalized 2D points

% Creating the interpolation space in 2D ---> Map values and 2D coordinates
total2Dpoints = height*width;
ind = (1:total2Dpoints);
s = [height, width];
[I,J] = ind2sub(s,ind);
F = scatteredInterpolant(I',J',double(img(:)));
intensity3Dpoints = F(uvMapPoints2D(:,1),uvMapPoints2D(:,2));
% For plotting
% grayMap = normalize(intensity3Dpoints, 1,0);
% RGBint = repmat(grayMap,1,3);
% scatter3(points3Dface(:,1),points3Dface(:,2),points3Dface(:,3),25,RGBint,'filled');

end

