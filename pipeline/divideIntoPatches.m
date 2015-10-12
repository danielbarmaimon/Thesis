function [ patches ] = divideIntoPatches( img, sideOfPatches )

% Inputs:
%           img = Input image to divide into patches (imgSize = nxn)
%           sideOfPatches = number of rows and columns of the patches image
% Output:
%           patches = Cell of square patches (sorted by row)

% Checking the size of the image
imSize = size(img);
if (mod(imSize(1),2)==0) % if it is even
    imCenter = imSize(1)/2;    
else % if it is odd
    imCenter = imSize(1)/2;    
    %imCenter = idivide(imSize(1),uint8(2)) + 1; 
end
% Calculating the size of each patch
sideOfPatch = idivide(uint32(imCenter), (sideOfPatches/2));
% Cropping the image
rowsToRemove = imCenter - double(sideOfPatch)*(sideOfPatches/2);
if (rowsToRemove<0)
    sideOfPatch = sideOfPatch-1;
    rowsToRemove = imCenter - double(sideOfPatch)*(sideOfPatches/2);
end
if(mod(rowsToRemove,2)==0.5)||(mod(rowsToRemove,2)==1.5)     %==0.5 % If # of rows to remove is even
        img(:,end:end)=[]; % Button
        img(end:end,:)=[]; % Right
        rowsToRemove = rowsToRemove - 0.5;
end
if( rowsToRemove > 0)               % Always remove from button to top and from right to left
    rowsToRemove = uint16(rowsToRemove);
    img(:,end-(rowsToRemove-1):end)=[]; % Button
    img(end-(rowsToRemove-1):end,:)=[]; % Right
    img(:,1:rowsToRemove)=[]; % Left
    img(1:rowsToRemove,:)=[]; % Top
end
patches = imsplit(img, [sideOfPatches sideOfPatches], 1:2);
end