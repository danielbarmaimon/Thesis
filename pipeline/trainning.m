%% Classification and training
% Folders
path = 'C:\Users\Iram\Desktop\BU4DSpontatneous\F006\test';

% Parameters
winSize = 10;       % Number of frames to consider an emotion
frameSide = 200;     % Size of the side of the square window to normalize the features

% Normalize the curvatures in size
normCurvSeqMat = zeros(frameSide,frameSide,nFrames);
for i = 1:nFrames
    normCurvSeqMat(:,:,i)= imresize(curvSeqMat(:,:,i), [frameSide,frameSide], 'bilinear');
end

% Get the emotion for each frame
AUsFile = csvread(dir(fullfile(path,'*T1.csv')));
frno = AUsFile(2:end,1);        % get all the frame numbers
codes = AUsFile(2:end, 2:end);  % get codes for all action units
emotions = zeros(1,7,nFrames);
for i = 1:nFrames
    AUs = codes(frno==i,:);     % get all AU for frame i
    emotions(:,:,i) = getEmtionFromAU(AUs);
end

% Get the main emotion for each window in time
emotionSum = zeros(size(emotions(:,:,1)));
emotionWindows = zeros(1,nFrames-(winSize - 1));
for i = 1:nFrames
    if(i<winSize)
        emotionSum = emotionSum + emotions(:,:,i);
    end
    emotionSum = emotionSum + emotions(i);
    [val, emotionWindows(i)] = max(emotionSum);
    emotionSum = emotionSum - emotions(i - (winSize - 1));
end
