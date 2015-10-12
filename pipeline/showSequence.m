function showSequence(seq)
% Shows a sequence stored in a matrix 
%   Inputs:
%           seq : Sequence to be shown in format (m x n x nframes)

name = 'video.avi';
% Create the video writer and open it
writerObj = VideoWriter(name);
open(writerObj);
% Define the figure frame
fig = figure;
% Set limits
topMax = max(max(max(seq)));
topMin = min(min(min(seq)));
for j =  1 : size(seq,3)
   imagesc(imrotate(seq(:,:,j),180));
   caxis([topMin  topMax]);
   title(j);
   colormap('jet');
   axis tight
   axis xy
   axis off
   frame = getframe ( fig );
   writeVideo(writerObj,frame);
end
close(writerObj);
delete(name);
end
