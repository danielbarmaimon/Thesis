function [ emotion ] = getEmotionFromAU(AUs)
% The definition for the emotions with respect to AUs are given in the
% paper 'BP4D Spontaneous: High resolution 3D dynamical database'
emotion = zeros(1,7);
if (AUs(:,12))      
    emotion = 1;          % Happiness
    return
elseif ((AUs(:,1) && AUs(:,2) && AUs(:,4)) ||(AUs(:,1) && AUs(:,2) && AUs(:,5)))
    emotion = 4;          % Fear
    return
elseif (AUs(:,4))  
    emotion = 2;          % Sadness
    return
elseif (AUs(:,5))
    emotion = 3;          % Surprise or startle
    return
elseif (AUs(:,23))
    emotion = 5;          % Anger or upset
    return
elseif (AUs(:,9) && AUs(:,10))
    emotion = 6;          % Disgust
    return
else
    emotion = 7;          % Neutral
    return
end

% emotion = zeros(1,7);
% if (AUs(:,6) && AUs(:,12))      
%     emotion(1,1) = 1;          % Happiness
%     return
% elseif (AUs(:,1) && AUs(:,4) && AUs(:,15))
%     emotion(1,2) = 1;          % Sadness
%     return
% elseif (AUs(:,1) && AUs(:,2) && AUs(:,5) && AUs(:,26))
%     emotion(1,3) = 1;          % Surprise
%     return
% elseif (AUs(:,1) && AUs(:,2) && AUs(:,4) && AUs(:,5) && AUs(:,7) && AUs(:,20) && AUs(:,26))
%     emotion(1,4) = 1;          % Fear
%     return
% elseif (AUs(:,4) && AUs(:,5) && AUs(:,7) && AUs(:,23))
%     emotion(1,5) = 1;          % Anger
%     return
% elseif (AUs(:,9) && AUs(:,15) && AUs(:,16))
%     emotion(1,6) = 1;          % Disgust
%     return
% else
%     emotion(1,7) = 1;          % Neutral
%     return
% end
end

