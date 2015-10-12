function [ SVMModels ] = train(mainPath, task)
% The function will train a classifier.
% Inputs:   mainPath = path to the main folder with all the data.
%                    file with the emotions labels ('emotions.mat')
%                    should be into the AU_OCC folder after running
%                    the 'emotionArrays' script
%           task =  the classifier will be trained for one subset of
%                   experiments, depending on the task that was performed
%                   This will be a string of the form 'T1' (1-8)
% Outputs:  classifier = the trained classifier (SVM) will be returned

% Establish the search parameter
searchParameter = strcat('*', task, '*');

% Get the labels
AUsPath = fullfile(mainPath, 'AU_OCC');
if (~exist(AUsPath,'file'))
    error('train:Folder', 'Folder with labels does not exist');
end
load(fullfile(AUsPath,'emotions.mat'));
allFiles = who('*T*');
list = who(searchParameter);
allFiles(ismember(allFiles,list))=[]; % Deleting labels of files not needed
clear(allFiles{:});

% Get the data and labels needed to train
featuresPath = fullfile(mainPath, 'curvature_results', 'features');
if (~exist(featuresPath,'file'))
    error('train:Folder', 'Folder with features does not exist');
end
d = dir(fullfile(featuresPath, searchParameter)); % Get the names for files 
nSubjects = numel(d);
% load(fullfile(featuresPath, d(1).name));        % Just for allocation
% featuresInFrame = size(featuresResult,1);       % of the features
labels = [];
values = [];
for i=1:nSubjects
    load(fullfile(featuresPath, d(i).name)); % Data
    values = [values, featuresResult];
    nameSplit = strsplit(d(i).name,'_');    % Correspondance data-label
    labelsVarName = strcat(nameSplit(1),'_', nameSplit(2));
    labels = [labels, eval(labelsVarName{:})]; % Labels
end
%% ---------- Support Vector Machines ----------------
fprintf('-----------Support Vector Machines ------------\n');
% ONLY 2 CLASSES  - Define the classifier parameters, train it and return
%%% SVMModels = fitcsvm(values',labels', 'KernelFunction', 'linear');
fprintf('--- 1. Automatic - Matlab does 10-fold cross validation over the model\n');
% MORE THAN 2 CLASSES
classes = unique(labels);
SVMModels = cell(numel(classes),1);
CVSVMModels = cell(numel(classes),1); % Just to check with original data
classesLoss = cell(numel(classes),1); % Just to check with original data
rng(1); % for reproducibility
for j = 1:numel(classes)
    indx =(labels==classes(j));
    SVMModels{j}= fitcsvm(values',indx','ClassNames', [false true], ...
        'Standardize', true, 'KernelFunction', 'linear');
    fprintf('Class %d, with %d elements in total', j, sum(indx));
    CVSVMModels{j} = crossval(SVMModels{j});
    classesLoss{j} = kfoldLoss(CVSVMModels{j})
end
% fprintf('--- 2. Manual - Ramdomly divide into 2 parts and perform the cross validation 10 times\n');
% % Setting the parameters
% nTimes = 10;
% 
% % Perform the validation of the model
% classes = unique(labels);
% trueLabels = zeros(1, numel(classes));
% falseLabels = zeros(1, numel(classes));
% trainData = [];
% testData = [];
% trainLabels = [];
% testLabels = [];
% for i = 1:nTimes
%     for j = 1:numel(classes) % Separate half of the data for each of the class
%         idxClass = find(labels==classes(j));
%         randIdxClass = randperm(numel(idxClass)); % Permutation of the indexes for a class
%         middle = idivide(numel(randIdxClass),int32(2),'round');
%         trainData = [trainData, values(:, randIdxClass(1:middle))];
%         testData = [testData, values(:, randIdxClass(middle+1:end))];
%         trainLabels = [trainLabels, labels(:, randIdxClass(1:middle))];
%         testLabels = [testLabels, labels(:, randIdxClass(middle+1:end))];
%     end
%     correctLabels = zeros(1,numel(classes));
%     incorrectLabels = zeros(1,numel(classes));    
%     for j = 1:numel(classes)
%         idx =(trainLabels==classes(j));  % Train a specific class
%         SVMModels{j}= fitcsvm(trainData',idx','ClassNames', [false true], ...
%         'Standardize', true, 'KernelFunction', 'linear');
%         idxTest = (testLabels==classes(j));
%         [labelsPredicted,scores] = predict(SVMModels{j},testData');
%         correctLabels(j)=sum( labelsPredicted == idxTest');
%         incorrectLabels(j)= numel(idxTest)-correctLabels(j);
%     end    
%     trueLabels = trueLabels + correctLabels;
%     falseLabels = falseLabels + incorrectLabels;
% end
% trueLabelsRatios = trueLabels./(trueLabels+falseLabels);
% falseLabelsRatios = falseLabels./(trueLabels+falseLabels)

%% ---------- Classification tree  ----------------
fprintf('-----------Classification Tree ---------------\n');
ctree = fitctree(values', labels');
resuberror = resubLoss(ctree)
cvrtree = crossval(ctree);
cvloss = kfoldLoss(cvrtree)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

