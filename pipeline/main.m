%% Script to compute the curvature maps after preprocessing
% This script is only to test changes in the preprocessing and compare
% them. The script to compute the final curvature maps is 'mainParallel'
% that calls the function 'processCurvature'.

%% Initialization of the window
% close all;  
% clc; 
% clear all;

tic;
%% Folder, files and configuration parameters
addpath (genpath('./toolbox_graph'));
addpath (genpath('./icp'));
%path = '';
%path = '';
%path = '';
%path = 'C:\Users\Daniel Barmaimon\Desktop\test\F010\T1';
%path = 'C:\Users\Daniel Barmaimon\Desktop\test5';
path ='C:\Users\Daniel Barmaimon\Desktop\test1\F001\T1';

% Setup flags and counters
savingFlag = true;         % True: save the video
mapping2DFlag = false;     % True: method to compute Viola-Jones and map
% Name of the file where we save our results and the video
results_main_folder = 'C:\Users\Iram\Desktop\BU4DSpontatneous\Sequences(2D+3D)\M001\test2\resultsDist0_50';

% Alpha = ResolutionROI / ResolutionGrid  ---> alpha >= 1
alpha = 1.00;
alphaStr = strrep(num2str(alpha,2),'.','_');

results_file_name = strcat(results_main_folder,alphaStr,'.mat');
videoName = strcat(results_main_folder,alphaStr,'.avi');

% Cropping the face once that is found, just for analysis
reduction = 0.6; % Amount of the plane over the face to consider

% Distance of the radius to consider for curvature approximation
rdist = 0.50;

%% Process the files
% List the documents needed
files = dir(fullfile(path,'*.obj'));
files = {files.name}';
nFrames = numel(files);
% Check if the sequence has been read or not
if exist(fullfile(path,'data.mat'),'file')
%    fprintf('Loading the storaged sequence...');    
    load(fullfile(path,'data.mat'));
%    fprintf(' done\n');
else        % Read the sequence and storage vertex and faces in '.mat' file
    if mapping2DFlag
        [ points3DfaceSeq, triSeq, uvSeq ] = readObjSeqUV( path, files, nFrames );
         %%%%% NEED TO READ AND STORAGE THE 2D TEXTURE IMAGES %%%%%%%%%%%%
    else
        [ points3DfaceSeq, triSeq] = readObjSeq( path, files, nFrames );
    end
end

% Getting the same amount of points for the faces in all the frames
%fprintf('Cropping the faces and getting same amount of points... ');    
points3DfaceSeqCrop = cell(1,nFrames);
if mapping2DFlag
    uvMapSeqCrop = cell(1,nFrames); %%% NEW
end
for i = 1:nFrames  % Check the maximum point in Y and remove all Z behind that
    points3Dface = points3DfaceSeq{i}';
    [val, idx] = max(points3Dface(:,2));
    listPixels = (points3Dface(:,3) > points3Dface(idx,3));
    points3DfaceSeqCrop{i} = points3Dface(listPixels,:)';
    if mapping2DFlag    
        uvMap = uvSeq{i}; %%%%%%%%%%%% NEW
        uvMapSeqCrop{i}=uvMap(:,listPixels); %%%%%%% NEW
    end
end
%fprintf(' done\n');
minPoints = min(cellfun('length',points3DfaceSeqCrop));
%fprintf('Checking minimum number of points...');    
%fprintf(' done\n');
%fprintf('Removing points assuring same number of points for each frame... ');
for i = 1:nFrames
    points3Dface = points3DfaceSeqCrop{i}';
    [copy, idx] = sort(points3Dface(:,3),'descend'); 
    numberOfPoints = minPoints;
    points3DfaceSeqCrop{i} = points3Dface(idx(1:numberOfPoints),:)';
    if mapping2DFlag
        uvMapSeqCrop{i}=uvMap(:,idx(1:numberOfPoints));  %%%%%%% NEW
    end
end
%fprintf('done\n');

% For each document, subsample, process and store
for i = 1:nFrames
%    fprintf('\n --------------------------\n');
%    fprintf('FRAME %d\n', i);
    % Get the data for the file
%    fprintf('Acquiring the frame... ');
    points3Dface = points3DfaceSeqCrop{i}';
    tri = triSeq{i};
%    fprintf('done\n');
    if mapping2DFlag
        uvMap = uvMapSeqCrop{i}';    %%%%%%% NEW
        img = images{i};
    end
    
    % Find the nose tip -------------- Careful (Check min or max) !!!!!!!!!
%    fprintf('Finding nose tip...');
    if mapping2DFlag
        [noseR, noseC, val] = noseFinding( points3Dface, uvMap, img )
    end
    %[val, idx] = min(points3Dface(:,3));
    [val, idx] = max(points3Dface(:,3));
    noseR = points3Dface(idx,1);
    noseC = points3Dface(idx,2); 
    %scatter3(noseR, noseC, val, 20, '*k');   
    
    % Assure the same amount of points
%     [idx, D] = knnsearch(points3Dface, [noseR, noseC,val], 'K', 24999);
%     points3Dface = [points3Dface(idx,:); [noseR, noseC,val]];
%     numberOfPoints = size(points3Dface,1);
%     
    % Data normalization
%    fprintf('Data normalization ... ');
%     points3Dface(:,1)= normalize(points3Dface(:,1), -1, 1);
%     points3Dface(:,2)= normalize(points3Dface(:,2), -1, 1);
%     points3Dface(:,3)= normalize(points3Dface(:,3), -1, 1);
%    fprintf('done\n');

[faceCroppedIdx,~] = rangesearch(points3Dface,[noseR, noseC, val],90);
%% X = points3Dface(faceCroppedIdx{:},1);
%% Y = points3Dface(faceCroppedIdx{:},2);
%% Z = points3Dface(faceCroppedIdx{:},3);
%% Xn = normalize(X, 1, -1);
%% Yn = normalize(Y, 1, -1);
%% Zn = normalize(Z, 1, -1);
%% principalDirections = pca([Xn,Yn,Zn]);
%% minY = min(Yn);  % Finding the nose approximately
%% maxY = max(Yn);
%% indices = (Yn(:)<0.4*maxY)&(Yn(:)>0.4*minY);
%% noseR = Xn(idx);
%% noseC = Yn(idx);
%% val = Zn(idx);

    
    % ICP
    if (i>1)
%        fprintf('Performing ICP ... ');
        [Ricp, Ticp] = icp(points3DfaceLast',points3Dface', 5, 'Matching', 'kDtree');
        points3Dface = Ricp * points3Dface' + repmat(Ticp, 1, numberOfPoints);
        points3Dface = points3Dface';
%        fprintf('done\n');
    end
    points3DfaceLast = points3Dface;
    % Grid parameters for curvature measurment
%    fprintf('Calculation of the grid... ');
    pointsGridX = round(sqrt(size(points3Dface,1)/alpha));
    pointsGridY = pointsGridX;
    pointsGridTotal = pointsGridX*pointsGridY;
%    fprintf('done\n');
    xFace = points3Dface(:,1);
    yFace = points3Dface(:,2);
    zFace = points3Dface(:,3);
    %scatter3(xFace,yFace,zFace,5,zFace,'filled');
    %hold on;    
    % Find the nose tip -------------- Careful (Check min or max) !!!!!!!!!
%    fprintf('Finding nose tip...');
    %%%%%%%%
    maxY = max(points3Dface(:,2));
    minY = min(points3Dface(:,2));
    % finding area to look for the nosetip (kind of cropping for hair problems)
    indices = (points3Dface(:,2)<0.4*maxY)&(points3Dface(:,2)>0.4*minY);
    [val, idx] = min(points3Dface(:,3).*indices);
    %%%%%%%%
    %[val, idx] = min(points3Dface(:,3));
    %%[val, idx] = max(points3Dface(:,3));
    noseR = points3Dface(idx,1);
    noseC = points3Dface(idx,2); 
    %scatter3(noseR, noseC, val, 20, '*k');
%    fprintf('done\n');
    % Plotting for the mesh and surface
       %plot_mesh(points3Dface, tri);
    % Calculating the princpal components to set the plane (PCA)
%    fprintf('Applying PCA...');
    principalDirections = pca([xFace,yFace,zFace]);
    Xprincipal = principalDirections(:,1);
    Yprincipal = principalDirections(:,2);
    Zprincipal = principalDirections(:,3);
%    fprintf('done\n');
    % Setting the coordinate system on the nose tip
   plot3([noseR,noseR + Xprincipal(1)], [noseC, noseC + Xprincipal(2)],[val,val + Xprincipal(3)], 'LineWidth', 5, 'Color', 'r'); % 3D coordinate axes
   plot3([noseR,noseR + Yprincipal(1)], [noseC, noseC + Yprincipal(2)],[val,val + Yprincipal(3)], 'LineWidth', 5, 'Color', 'b');
   plot3([noseR,noseR - Zprincipal(1)], [noseC, noseC - Zprincipal(2)],[val,val - Zprincipal(3)], 'LineWidth', 5, 'Color', 'g');
%    hold on;

    % Plane representation Ax + By + Cz + D = 0
%    fprintf('Creation of the plane with grid over nose tip...');
    [ Xplane,Yplane, Zplane] = createPlaneGridReduced( Zprincipal, xFace, yFace, noseR, noseC, val, pointsGridX, pointsGridY, reduction);
%    S = surf(Xplane, Yplane, Zplane, 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'FaceColor', [0 0 1]);
    % Representing the grid
    XplaneLine = Xplane(:);
    YplaneLine = Yplane(:);
    ZplaneLine = Zplane(:);
%    scatter3(XplaneLine,YplaneLine,ZplaneLine,15,'red','filled');
%    fprintf('done\n');
    % Interpolating over the surface, around the grid projection
%    fprintf('Interpolation over surface with grid...');
    F = scatteredInterpolant(xFace,yFace,zFace);  %%%% , 'nearest'PROBLEM HERE .....// Points go to the back the data get corrupted
    Zq = F(Xplane,Yplane);
    bigger = (Zq > 1);          % limit values out of the range
    lower = (Zq < -1);
    Zq(bigger)=1;
    Zq(lower)=-1;        
    %a = plot3(Xplane,Yplane,Zq,'o', 'MarkerSize', 5);     
%    fprintf('done\n');
    % Finding the nearest neighbours to each of the estimated points on surf.
%    fprintf('Calculation of the curvature...');
    curvPoints = [XplaneLine, YplaneLine, Zq(:)];
    initialPoints = [xFace, yFace, zFace];
    %[K,H,Pmax,Pmin] = surfature(Xplane,Yplane,Zq);
    %imagesc(K)
    %%%%%%%%%%%%%%
    [curvNeigh, distNeigh] = rangesearch(initialPoints,curvPoints,rdist);
    nNeigh =round(mean(cellfun(@numel,curvNeigh)));
    %%%%%%%%%
    [ normalsGrid, curvatureGrid] = findPointNormals(curvPoints, nNeigh, [noseR, noseC, 500]);
    %[ normalsGrid, curvatureGrid] = findPointNormals(curvPoints, 9, [noseR, noseC, 500]);
    curvatureFrame = reshape(curvatureGrid,[pointsGridY,pointsGridX]); % CAREFUL with X and Y
%    fprintf('done\n');
    % Removing outliers
%    fprintf('Removing outliers... ');
    %curvatureFrame = removeOutliers(curvatureFrame, 97, 200);
    % Storing the curvature frame in the sequence (vectorized)
    curvatureSeq(:,i) = curvatureFrame(:);
%    fprintf('done\n');
    % Approach 2 - Neighbours, filter by euclidian distance, filter by adj
    %[idx, D] = knnsearch(initialPoints, curvPoints, 'K', 50);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     [curvNeigh, distNeigh] = rangesearch(initialPoints,curvPoints,rdist);
%     % Curvature and normals for all the face
%     [ normalsFace, curvatureFace ] = findPointNormals(initialPoints, 8, [noseR, noseC, -500], false);   % normals and curvature over all face [xnose3D, ynose3D, val]
%     % Curvature and normals for the grid
%     [ normalsGrid, curvatureGrid ] = interpolateNormals( normalsFace, curvatureFace, pointsGridTotal, rdist, curvNeigh, distNeigh);
    
    % Estimated normals over interpolated points
%     hold on;
%     quiver3(Xplane(:),Yplane(:),Zq(:),...
%         normalsGrid(:,1),normalsGrid(:,2),normalsGrid(:,3),'r');
%     axis equal;
    
    %All the normals
%     figure;
%     %surf(xFace,yFace,repmat(zFace,1,size(xFace,1)));%,curvature
%     hold on;
%     quiver3(xFace(:),yFace(:),zFace(:),...
%         normalsFace(:,1),normalsFace(:,2),normalsFace(:,3),'r');
%     axis equal;
    
    % Asigning each normal and curvature matrix to a frame of the seq
%    curvatureSeq(:,i) = curvatureGrid;
%    normalSeq(:,:,i) = normalsGrid;
end
time = toc

% Curvature representation -- From dots to interpolated image
%fprintf('Preparing the representation and video...\n');
curvSeqMatrix = reshape(curvatureSeq, [pointsGridX, pointsGridY, nFrames]);
% 2D -- Bicubic interpolation
%scaleInterp = 4;
%curvRep = zeros(scaleInterp*pointsGridX, scaleInterp*pointsGridY, 2);
% for i = 1:nFrames
%     %curvRep(:,:, i) = imresize(curvSeqMatrix(:,:,i), scaleInterp, 'bicubic');
%     curvSeqMatrix(:,:,i)=imrotate(curvSeqMatrix(:,:,i),90);
% end
%%%%%%%%%%%%%%% filtering the output %%%%%%%%%
% for i =1:nFrames
%     curvSeqMatrix(:,:,i) = medfilt2(curvSeqMatrix(:,:,i));
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (savingFlag)
    % Save in a video format
    topMax = max(max(max(curvatureSeq)));
     %create_avi(curvRep(:,:,1:nFrames),15,0,0.6,videoName); % if bicubic
     create_avi(curvSeqMatrix(:,:,1:nFrames),15,0,topMax,videoName); % no interpolation   
     %save(results_file_name, 'curvatureSeq', 'curvSeqMatrix', 'normalSeq');
     save(results_file_name, 'curvatureSeq', 'curvSeqMatrix');
end
%fprintf('done\n');
% 3D -- Surface representation
%fprintf('Showing a 3D curvature representation...\n');
frmsToShow = 5;
frmsStep = nFrames/uint8(frmsToShow); %nFrames
dirName = strcat(results_main_folder,alphaStr,'\');
mkdir(dirName);
for i = 1:frmsToShow
    a = figure;
    surf(curvSeqMatrix(:,:,i));
    caxis([0 topMax]);
    title(['Frame ' num2str(i+(i-1)*frmsStep)]);
    savefig(a, strcat(dirName,num2str(i+(i-1)*frmsStep),'.fig'));
end
%fprintf('done\n');