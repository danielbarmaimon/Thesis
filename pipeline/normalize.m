function normalized = normalize( array, maxRange, minRange )
%NORMALIZE Summary of this function goes here

% Normalization [0 , 1]
minArray = min(array);
maxArray = max(array);
initRange = maxArray - minArray;
initNorm = (array - minArray)./ initRange;

% Adjust to desired range
finalRange = maxRange - minRange;
normalized = (initNorm * finalRange) + minRange;
end

