function depthPoints( path, path_result, alpha, reduction, savingFlag  )
% Get the sequence of sorted depth points frames for computer
% super-resolution over them.

%% Process the files
% List the documents needed
files = dir(fullfile(path,'*.obj'));
files = {files.name}';
nFrames = numel(files);
% Check if the sequence has been read or not
if exist(fullfile(path,'data.mat'),'file')
%    fprintf('Loading the storaged sequence...');    
    load(fullfile(path,'data.mat'));
    nFrames = size(points3DfaceSeq,2);
%    fprintf(' done\n');
else        % Read the sequence and storage vertex and faces in '.mat' file
    [ points3DfaceSeq, triSeq] = readObjSeq( path, files, nFrames );
    %[ points3DfaceSeq, triSeq, UVmap ] = readObjSeqUV( path, files, nFrames );
    dataName = fullfile(path,'data.mat');
    save(dataName,'points3DfaceSeq', 'triSeq');    
end

% Getting the same amount of points for the faces in all the frames
%fprintf('Cropping the faces and getting same amount of points... ');    
points3DfaceSeqCrop = cell(1,nFrames);
for i = 1:nFrames  % Check the maximum point in Y and remove all Z behind that
    points3Dface = points3DfaceSeq{i}';
    [val, idx] = max(points3Dface(:,2));
    listPixels = (points3Dface(:,3) > points3Dface(idx,3));
    points3DfaceSeqCrop{i} = points3Dface(listPixels,:)';
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
    
    % Find the nose tip -------------- Careful (Check min or max) !!!!!!!!!
%%    fprintf('Finding nose tip...');
%     %[val, idx] = min(points3Dface(:,3));
%     [val, idx] = max(points3Dface(:,3));
%     noseR = points3Dface(idx,1);
%     noseC = points3Dface(idx,2); 
%     %scatter3(noseR, noseC, val, 20, '*k');   

    % Data normalization
%    fprintf('Data normalization ... ');
    points3Dface(:,1)= normalize(points3Dface(:,1), -1, 1);
    points3Dface(:,2)= normalize(points3Dface(:,2), -1, 1);
    points3Dface(:,3)= normalize(points3Dface(:,3), -1, 1);
%    fprintf('done\n');
    
    % ICP
    if (i>1)
%        fprintf('Performing ICP ... ');
        [Ricp, Ticp] = icp(points3DfaceLast',points3Dface', 9, 'Matching', 'kDtree');
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
    %[val, idx] = max(points3Dface(:,3));
    noseR = points3Dface(idx,1);
    noseC = points3Dface(idx,2); 
    %scatter3(noseR, noseC, val, 20, '*k');
%    fprintf('done\n');
    % Plotting for the mesh and surface
       %plot_mesh(points3Dface, tri);
    % Calculating the princpal components to set the plane (PCA)
%    fprintf('Applying PCA...');
    principalDirections = pca([xFace,yFace,zFace]);
    Yprincipal = principalDirections(:,1);
    Xprincipal = principalDirections(:,2);
    Zprincipal = principalDirections(:,3);
%    fprintf('done\n');
    % Setting the coordinate system on the nose tip
%    plot3([noseR,noseR + Xprincipal(1)], [noseC, noseC + Xprincipal(2)],[val,val + Xprincipal(3)], 'LineWidth', 1, 'Color', 'r'); % 3D coordinate axes
%    plot3([noseR,noseR + Yprincipal(1)], [noseC, noseC + Yprincipal(2)],[val,val + Yprincipal(3)], 'LineWidth', 1, 'Color', 'b');
%    plot3([noseR,noseR + Zprincipal(1)], [noseC, noseC + Zprincipal(2)],[val,val + Zprincipal(3)], 'LineWidth', 1, 'Color', 'g');
%    hold on;

    % Plane representation Ax + By + Cz + D = 0
%    fprintf('Creation of the plane with grid over nose tip...');
    [ Xplane,Yplane, Zplane] = createPlaneGridReduced( Zprincipal, xFace, yFace, noseR, noseC, val, pointsGridX, pointsGridY, reduction);
%    S = surf(Xplane, Yplane, Zplane, 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'FaceColor', [0 0 1]);
    % Representing the grid
    XplaneLine = Xplane(:);
    YplaneLine = Yplane(:);
    ZplaneLine = Zplane(:);
%    scatter3(XplaneLine,YplaneLine,ZplaneLine,0.5,'red','filled');
%    fprintf('done\n');
    % Interpolating over the surface, around the grid projection
%    fprintf('Interpolation over surface with grid...');
    F = scatteredInterpolant(xFace,yFace,zFace);  %%%% , 'nearest'PROBLEM HERE .....// Points go to the back the data get corrupted
    Zq = F(Xplane,Yplane);
    bigger = (Zq > 1);          % limit values out of the range
    lower = (Zq < -1);
    Zq(bigger)=1;
    Zq(lower)=-1;        
    %plot3(Xplane,Yplane,Zq,'o', 'MarkerSize', 1);     
%    fprintf('done\n');  
    curvPoints = [XplaneLine, YplaneLine, Zq(:)];
    curvPointsSeq(:,:,:,i) = curvPoints;  
    sizeFrame = [pointsGridY pointsGridX];
%    fprintf('done\n');
end

% Storing the results
%fprintf('Saving the curvature values...');
if (savingFlag)  % NOT WITH PARALEL
    save(path_result, 'nFrames', 'sizeFrame','curvPointsSeq');   
end
%fprintf('done\n');
end

