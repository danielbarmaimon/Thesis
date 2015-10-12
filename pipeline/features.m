function [ features ] = features( fileFullPath, nPatchesSide, lower,upper,nbins, normalize)
% The function retuns the features of a sequence of curvature frames. This
% features will represent the histogram of different patches of the
% curvature image. In this way, there will be, per sequence:
% # features = nPatchesSide*nPatchesSide*nbins*nFrames 
%
%   Inputs: fileFullPath = path of the current sequence to process
%           nPatchesSide = number of divisions for the image across X or Y
%           nbins        = number of bins of the histogram for each patch
%           lower        = lower value for the histogram computation
%           upper        = upper value for the histogram computation
%           normalize    = normalization of the histogram (should be true)
%   Outputs: features    = matrix of features [nFrames x (nbins*nPatches)]

%% Load the file with the sequence of curvature and size of the frame
load(fileFullPath);

%% Define the histogram object
his = vision.Histogram(lower,upper,nbins,'Normalize',normalize);

%% Get the features for a sequence of frames
nFrames = size(curvatureSeq,2);
nPatches = nPatchesSide^2;
curvatureFrameSeq = reshape(curvatureSeq, [sizeFrame, nFrames]);
features = zeros(nPatches*nbins, nFrames);
for i=1:nFrames
    patches = divideIntoPatches( curvatureFrameSeq(:,:,i), nPatchesSide );
    frame = [];
    for j=1:nPatches
        patchFeat = step(his, patches{j});
        frame = [frame; patchFeat];
    end
    features(:,i)=frame;
end
release(his);
end

