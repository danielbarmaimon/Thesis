
%% Parameters
% Folders
path = 'C:\Users\Daniel Barmaimon\Desktop\test1';
% Parameters of the frame
nPatchesSide = 5; % Number of the side of the square of patches
% Parameters of the histogram
lower = 0;
upper = 0.3; % TO CHECK
nbins = 10;  % 
normalize = true; % Assure that histogram values for bins add to 1

% Parameters of the window (if needed)
winSize = 10;       % Number of frames to consider an emotion

%% Directories and files management
curvatureResultsFullPath = fullfile(path, 'curvature_results');
if (~exist(curvatureResultsFullPath,'file'))
    error('classif:Folder', 'Folder with results does not exist');
end
path_result = fullfile(curvatureResultsFullPath, 'features');
% Create a directory for the features results if it doesn't exist
if (~exist(path_result,'file'))
    mkdir(path_result);
end
% Get all the folders for the same category
d = dir(curvatureResultsFullPath);
% Subject of the database
isubd = [d(:).isdir]; % logical vector
subd = {d(isubd).name}';
subd(ismember(subd,{'.','..','curvature_results','features'}))=[]; % remove . and ..
nSubjects = numel(subd);
for subject = 1:nSubjects
    subjectName = subd{subject};
    dd = dir(fullfile(curvatureResultsFullPath, subjectName,'*.mat'));
    nTasks = numel(dd);
    %% Get the features for a sequence of frames
    for task = 1:nTasks 
        fileNameSplit = strsplit(dd(task).name, '.');
        fileName = fileNameSplit{1};
        fileNameResult = strcat(fileName,'_feat.mat');
        path_result = fullfile(curvatureResultsFullPath, 'features', fileNameResult);
        file = fullfile(curvatureResultsFullPath,subjectName,dd(task).name);
        featuresResult = features(file, nPatchesSide, lower, upper, nbins, normalize);
        save(path_result, 'featuresResult');
    end
end


% nFrames = size(curvatureSeq,2);
% curvatureFrameSeq = reshape(curvatureSeq, [sizeFrame, nFrames]);
% features = zeros(nPatches*nbins, nFrames);
% for i=1:nFrames
%     patches = divideIntoPatches( curvatureFrameSeq(:,:,i), nPatchesSide );
%     frame = [];
%     for j=1:nPatches
%         patchFeat = step(his, patches{j});
%         frame = [frame; patchFeat];
%     end
%     features(:,i)=frame;
% end
% % Get the main emotion for each window in time
% emotionSum = zeros(size(emotions(:,:,1)));
% emotionWindows = zeros(1,nFrames-(winSize - 1));
% for i = 1:nFrames
%     if(i<winSize)
%         emotionSum = emotionSum + emotions(:,:,i);
%     end
%     emotionSum = emotionSum + emotions(i);
%     [val, emotionWindows(i)] = max(emotionSum);
%     emotionSum = emotionSum - emotions(i - (winSize - 1));
% end