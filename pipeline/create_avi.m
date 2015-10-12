function create_avi( seq,lower,upper, name)
% This process creates a video sequence given the frames.
% Inputs:
%           seq : Sequence of images as a matrix with dimensions M x N x nFrames
%           lower / upper : Limits for the representation (maintaining constant values for the colors)
%           name : name of the video file that will be created

% Create the video writer and open it
writerObj = VideoWriter(name);
open(writerObj);

% Define the figure frame
fig = figure;

% Display frame by frame
for j =  1 : size(seq,3)
   imagesc(imrotate(seq(:,:,j),180));
   caxis([lower  upper]);
   title(j);
   colormap('jet');
   axis tight
   axis xy
   axis off
   frame = getframe ( fig );
   writeVideo(writerObj,frame);
end
close(writerObj);
end

