function [ Xplane,Yplane, Zplane] = createPlaneGridReduced( Zprincipal, xFace, yFace, xNose3D, yNose3D, zNose3D, pointsGridX, pointsGridY, reduction)
% Creates a gridded plane over the plane in front of the face. It uses the
% nose location in 3D and the principal components of the point cloud for
% the plane generation. The amount of points in X and Y directions should
% be given (usually is the same to have square images later on). The last  
% parameter regulates the reduction of the face.
% Inputs:
%           Zprincipal : Principal direction that matches with the normal to the plane
%           xFace, yFace : List of x and y coordinates for the 3D facial points 
%           xNose3D, yNose3D, zNose3D : Location of the nose
%           pointsGridX, pointsGridY : Number of points in x and y for the grid over the plane
%           reduction : Plane is limited by (minX, maxX) and (minY, maxY). This parameter will limit these for values by multiplying them by a factor between (0, 1)
%
% Outputs:
%           Xplane, Yplane, Zplane : Points of the grid in 3D
%
% NOTE: The function is considering that the face is aligned with two of the
% principal components in X and Y axis (or almost). This is very convinient
% due to the function 'meshgrid' that allows the creation of a grid in x, y 
% standard coordinates.
%
% TO DO (improvement):
% For a more precise result a homography should be calculated 
% align the facial points with respect to x, y standard coordinate system.
% After alignment grid should be created. And as a last step, all the
% points in the grid should be back to the xprincipal, yprincipal
% coordinate system using the inverse of the homography.

% Plane representation Ax + By + Cz + D = 0
A = Zprincipal(1);
B = Zprincipal(2);
C = Zprincipal(3);
D = -dot(Zprincipal, [xNose3D, yNose3D, zNose3D]);
xlimMax = reduction * max(max(xFace));
ylimMax = reduction * max(max(yFace));
xlimMin = reduction * min(min(xFace));
ylimMin = reduction * min(min(yFace));
xlimMax = xlimMax + (xNose3D - (xlimMax+xlimMin)/2);    % Centering in the nose
xlimMin = xlimMin + (xNose3D - (xlimMax+xlimMin)/2);    % Centering in the nose
ylimMax = ylimMax + (yNose3D - (ylimMax+ylimMin)/2);    % Centering in the nose
ylimMin = ylimMin + (yNose3D - (ylimMax+ylimMin)/2);    % Centering in the nose
xx = linspace(xlimMin, xlimMax,pointsGridX);
yy = linspace(ylimMin, ylimMax,pointsGridY);
[Xplane,Yplane] = meshgrid(xx,yy);
Zplane = (A*Xplane + B*Yplane + D)/(-C);

end

