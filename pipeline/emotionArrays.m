function emotionArrays(mainPath)
% This function will collect all the Action Units from the '.csv' files and 
% arrange them into arrays with different values for each emotion. The
% emotions will be saved in different variables depending the subject and
% the task. All of them will be storaged in a file called 'emotions.mat'

% Get all the names of the files
d = dir(fullfile(mainPath,'AU_OCC','*.csv'));
idx = [d(:).isdir]; % logical vector
filesNames = {d(~idx).name}';
nFiles = numel(filesNames);
pathResult = fullfile(mainPath,'AU_OCC','emotions.mat');
% Read each of the '.csv'
for i=1:nFiles
    name = filesNames{i};
    AUsFile = csvread(fullfile(mainPath,'AU_OCC',name));
    nameSplit = strsplit(name,'.');
    name = nameSplit(1);
    codes = AUsFile(2:end, 2:end);  % get codes for all action units
    nFrames = size(AUsFile,1)-1;
    emotions = zeros(1,nFrames);
    for j = 1:nFrames
        AUs = codes(j,:);     % get all AU for frame i
        emotions(j) = getEmotionFromAU(AUs);
    end
    eval([name{:} '=emotions;']);
    if(exist(pathResult, 'file'))
        save(pathResult, name{:}, '-append');
    else
        save(pathResult, name{:});
    end
end
end