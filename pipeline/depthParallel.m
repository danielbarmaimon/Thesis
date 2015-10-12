%% This script generates the depth sequences as a PMD camera would do (sorted grid).
% The parameters to set are listed in Configuration Parameters section

%% Initialization of the window
close all;  
clc; 
clear all;

tic;
%% Configuration parameters
% Setup flags and counters
savingFlag = true;         % True: save the sequence

% Process parameters
alpha = 1.00;       % Resolution --> Min = 1 (maximum resolution).
reduction = 0.6;    % Cropping of the face

%% Folder and files
addpath (genpath('./toolbox_graph'));
addpath (genpath('./icp'));
% Main folder of the same category
mainPath = 'C:\Users\Daniel Barmaimon\Desktop\test2\';
curvatureResultsFullPath = fullfile(mainPath, 'curvature_results');
if (~exist(curvatureResultsFullPath,'file'))
    mkdir(curvatureResultsFullPath);
end

%% Processing and getting the curvatures
% Get all the folders for the same category
d = dir(mainPath);
% Subject of the database
isubd = [d(:).isdir]; % logical vector
subd = {d(isubd).name}';
subd(ismember(subd,{'.','..','curvature_results'}))=[]; % remove . and ..
nSubjects = numel(subd);
for subject = 1:nSubjects
% Task on the database
   subjectName = subd{subject};
   dd = dir(fullfile(mainPath, subjectName));
   isubdd = [dd(:).isdir]; % logical vector
   subdd = {dd(isubdd).name}';
   subdd(ismember(subdd,{'.','..','curvature_results'}))=[]; % remove . and ..
   nTasks = numel(subdd);
   path = repmat({''},nTasks,1);
   path_result = repmat({''},nTasks,1);
   nameFile = repmat({''},nTasks,1);
   curvatureSeq = cell(1,nTasks);
   sizeFrame = cell(1,nTasks);
   % Working in parallel for each task of the same subject
   for (p = 1:nTasks)    % for all the tasks % IT GETS OUT OF MEMORY FOR rdist = 0.5 ---> change 'parfor' by 'for'
        % Full path of the sequence
        path{p} = fullfile(mainPath, subjectName,subdd{p});
        % Name for the resulting file and folder
        alphaStr = strrep(num2str(alpha,'%1.2f'),'.','_');
        nameFile{p} = strcat(subjectName,'_',subdd{p},'_a_',alphaStr,'.mat');
        path_result{p} = fullfile(mainPath, 'depth_results', subjectName);
        if exist(path{p}, 'file')
            %%%%%%%%%%% CODE TO PERFORM %%%%%%%%%%%%%%%
            path{p}
            % 1.- Process and save the curvature sequence
            if (~exist(path_result{p},'file'))
                mkdir(path_result{p});
            end
            path_result{p} = fullfile(mainPath, 'depth_results', subjectName,nameFile{p});
            depthPoints( path{p}, path_result{p}, alpha, reduction, true);
        end
    end
end