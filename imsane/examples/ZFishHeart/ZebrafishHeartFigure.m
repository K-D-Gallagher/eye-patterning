%% ImSAnE: Drosophila melanogaster embryo tutorial
%
% In this tutorial we detect the embryos apical surface during
% cellularization, fit a spherelike surface, generate a surface of interest
% and pullback the image data to various charts.
%
% This an example of how to detect and fit what we call spherelike surfaces,
% where the surface can be described as a slowly varying function of a
% preferred axis, which we call the z-axis. It also demonstrates using
% batch processing of multiple time points with jitter correction. Moreover
% it shows how to choose from two different detectors and describes the 
% options of loading an externally provided point cloud.
%
% Here is how to access the documentation of detectors and fitters that are
% introduced in this tutorial

%%
% 
% doc surfaceDetection.fastCylinderDetector
% doc surfaceDetection.radialEdgeDetector
%
% doc surfaceFitting.spherelikeFitter


%% Initialize ImSAnE project
%
% We start by clearing the memory and closing all figures.
%
clear all; close all;
%
% Setup a working directory for the project, where extracted surfaces,
% metadata and debugging output will be stored. Also specify the directory
% containing the data. 
%
[scriptPath,~,ext] = fileparts(matlab.desktop.editor.getActiveFilename);

%dataDir = fullfile(scriptPath, 'rawData');
dataDir = '/Users/idse/Dropbox/flows/flows_shared/data/zebrafish/huisken/video9/';
%'/Users/idse/Dropbox/flows/flows_shared/codes/ImSAnE_NatureMethods/examples/DrosoEmbryo/rawData/';
%projectDir = fullfile(scriptPath, 'projectFiles');
projectDir = dataDir;
%
% Start by creating an experiment object, optionally pass on the project
% directory (otherwise it will ask), and change into the directory of the
% data. This serves as a frontend for data loading, detection, fitting etc.
%
xp = project.Experiment(projectDir, dataDir);

% optionally switch working directory into the data directory.
%
cd(dataDir)

% Set file and experiment metadata
%
% Set required additional information on the files
% 
% We assume one individual image stack for each time point, labeled by time.
% To be able to load the stack, we need to tell the project where the data
% is, what convention is assumed for the file names, available time points
% and the stack resolution. Options to modules in ImSAnE are organised in
% matlab structures, that is a pair of field name and value are provided
% for each option. 
%
% The following file metadata information is required   
%
% * 'directory'       , the project directory (full path)
% * 'dataDir'         , the data directory (full path)
% * 'filenameFormat'  , fprintf type format spec of file name
% * 'timePoints'      , list of times available stored as a vector
% * 'stackResolution' , stackresuolution in microns, e.g. [.25 .25 1]
%
%

%
% The following file metadata information is optional
%
% * 'imageSpace'      , bit depth of image, such as uint16 etc, defined in Stack class
% * 'stackSize'       , size of stack in pixels per dimension [xSize ySize zSize]
% * 'swapZT'          , for some datasets time is the third dimension, and z 
% the fourth. Swap = 1 if this is the case. 
%
% This tutorial uses a SPIM dataset, that is 2 fold downsampled, for
% speedup and memory requirements. There is also a full size data file
% included, however, it is strongly recommended to have at least 6 GB of
% free memory when trying this dataset. Moreover, some fileMeta and
% detector options should be changed, where indicated. 
fileMeta                 = struct();
fileMeta.dataDir         = dataDir;
fileMeta.filenameFormat  = 'VSilent_T%d_G.tif';%'Time%06d_bin2.tif'; % for full data sample use Time000000.tif
fileMeta.timePoints      = 30:69; % for full data sample use 0;
fileMeta.stackResolution = [.45 .45 1];%[.5 .5 .5]; 
fileMeta.swapZT          = 0; % for full data sample use 1;

%
% Set required additional information on the experiment. A verbal data set
% description, Jitter correct by translating the sample, which time point to 
% use for fitting, etc.
%
% The following project metadata information is required 
%
% * 'channelsUsed'   , the channels used, e.g. [1 3] for RGB
% * 'channelColor'   , mapping from element in channels used to RGB = 123
% * 'dynamicSurface' , Not implmented yet, future plan: boolean, false: static surface
% * 'detectorType'   , name of detector class, e.g. radialEdgeDetector
%                        (user thresholded), fastCylinderDetector
% * 'fitterType'     , name of fitter class, here spherelikeFitter
% * 'fitTime'        , time point used for fit
%
% The following project metadata information is optional 
%
% * 'description'     , string describing the data set set experiments metadata, 
%                                such as a description, and if the surface is dynamic,
%                                or requires drift correction of the sample.
% * 'jitterCorrection', Boolean, false: No fft based jitter correction 
%
expMeta                  = struct();
expMeta.channelsUsed     = 1;
expMeta.channelColor     = 1;
expMeta.description      = 'Beating Zebrafish Heart';
expMeta.dynamicSurface   = 1;
expMeta.jitterCorrection = 0; % 1: Correct for sample translation
expMeta.fitTime          = 60; 
expMeta.detectorType     = 'surfaceDetection.fastCylinderDetector';%
%expMeta.fitterType       = 'surfaceFitting.cylinderMeshWrapper';
expMeta.fitterType       = 'surfaceFitting.meshWrapper'; 

% Now set the meta data in the experiment.
%
xp.setFileMeta(fileMeta);
xp.setExpMeta(expMeta);

%
% For a new project call initNew()
% This reads the stack size from the first available time point, then
% initialize fitter and detector and create fitOptions and detectOptions
% based on their defaults.
%
xp.initNew();

%% Load data for surface detection and rescale to unit aspect ratio
%
% Load a timepoint from the data.
% loadTime sets currentTime, loads the stack and resets detector and fitter
% with the appropriate options for that time. Optionally rescale the stack 
% to unit aspect ratio if desired.
%
xp.loadTime(60); %fileMeta.timePoints(1)

%
% The stack used in this tutorial already has unit aspect ratio,
% therefore rescaleStackToUnitAspect will do nothing. However, most
% datasets won't, which will cause detectors to work less well. Then 
% rescaling axes to unit aspect ratio will be useful.
%
xp.rescaleStackToUnitAspect();

%% Detect the surface
%
% Before detecting the surface we set the options of the used detector.
% Changing the detectOptions resets the detector so that one cannot have
% for example a detected pointcloud that was detected with other parameters
% than the current one. 
%
% Default mode of this tutorial is the fast cylinder detector. Below are
% the options to set for the radialEdgeDetector. Don't forget to switch the
% detector type in the expMeta structure, when trying this out.
%
%% 
%
% options for fast cylinder detector
% Standard filters are applied to the images for point cloud detection. The
% options used are 
%
% * 'channel'             , vector of integers, specify which channel(s) to use for detection 
% * 'sigma'               , standard deviation of a gaussian filter in pixels 
% * 'ssfactor'            , integer, specifying the degree of data subsampling. 
% * 'nBins'               , number of radial bins to determine the point cloud in
% * 'rmRadialOutliers'    , remove radial outliers, low means stringent, high sloppy, 0 is off.
% * 'rmIntensityOutliers' , remove outliers based on intensity. 
% * 'zDim'                , specify long axis in data.
%
%%
%
% Comment below, when using radial edge detector 
% for full data sample change ssfactor to 4;
myDetectOpts = struct('channel', 1, 'sigma', 1, 'ssfactor', 2, ...
    'nBins', 120,'rmRadialOutliers', 2, 'rmIntensityOutliers',2,...
    'zDim', 2);  
%%
%
% options for radial edge detector
%
% * 'thresh' , determine cutoff for foreground detection
% * 'amin'   , minimal number of connected pixels to be foreground
% * 'bgdisc' , diameter of rolling ball filter for background estimation
% * 'dildisc', diameter of dilation disc for morphological opening and closing.
% * 'sp'     , permute stack dimensions before detection.
%
%%
%
%  Uncomment below, when using radial edge detector
%
%myDetectOpts = struct('channel', 1, 'sigma',1, 'ssfactor', 4,...
%   'rmRadialOutliers', 1.2, 'rmIntensityOutliers',2, 'thresh',75,...
%   'amin',5,'bgdisc',0,'dildisc',20,'sp',[1 3 2]);  

%%
%
% set the detect options in the project
%
xp.setDetectOptions(myDetectOpts);
%%
%
% calling detectSurface runs the surface detector and creates
% the detector.pointCloud object
%
xp.detectSurface();

%% Inspect the point cloud in a cross section
%
% inspect point cloud over a cross section inthe  data. Dimensions 
% are x,y or z and the value has to be within the corresponding axis range.
%
inspectOptions= struct('dimension', 'z', 'value', 200, 'pointCloud', 'b');
xp.detector.inspectQuality(inspectOptions, xp.stack);

%% optionally save the quality inspection to disc
%
%storeQualityOptions = struct('inspectOptions',inspectOptions,...
% 'range',1:2:200,'outName',fullfile(projectDir,...
% 'debugOutput/InspectionDetector.tif'),'export','true','closeFig','true');
%xp.detector.storeQualityInspection(storeQualityOptions, xp.stack)

%% Inspect pointcloud in 3d
%
% Plot the points on the surface as a three dimensional point cloud.
% The subsampling factor reduces the number of points shown. 
%
ssfactor = 6;
xp.detector.pointCloud.inspect(ssfactor);

%% Save pointCloud to obj format

clear OBJ
OBJ.vertices = xp.detector.pointCloud.unalignedPoints;
OBJ.objects(1).type='f';
OBJ.objects(1).data.vertices=[];
write_wobj(OBJ, fullfile(projectDir, 'pointCloud.obj'));

%% Copy meshlab command to clipboard

% if meshlabserver is in the bash path, just paste in the terminal and
% press enter

% add to path by making a symbolic link
% sudo ln -s /Applications/meshlab.app/Contents/MacOS/meshlabserver /usr/bin/meshlabserver

PCfile = fullfile(projectDir, 'obj3', ['heartPC',num2str(xp.currentTime),'.obj']);
outputMesh = fullfile(projectDir, 'obj3',['heartPC',num2str(xp.currentTime),'_mesh.obj']);
meshlabScript = fullfile(projectDir, 'heartTest2.mlx');
clipboard('copy', ['meshlabserver -i ' PCfile ' -o ' outputMesh ' -s ' meshlabScript ' -om vn']);

%% 
%PCfile = fullfile(projectDir, 'obj', ['heartPC',num2str(69),'.obj']);
%outputMesh = fullfile(projectDir, 'obj',['heartPC',num2str(69),'_mesh.obj']);
%meshlabScript = fullfile(projectDir, 'heartTest2.mlx');
%clipboard('copy', ['meshlabserver -i ' PCfile ' -o ' outputMesh ' -s ' meshlabScript ' -om vn']);

%% Read surface mesh produced by meshlab
%outputMesh = fullfile(projectDir, 'obj','heartPC95_mesh.obj');
tic
mesh = readObj(outputMesh);
mesh.f = mesh.f.v;
toc

% dumb stuff
mesh.v = mesh.v(:,[3,1,2]);
mesh.vn = mesh.vn(:,[3,1,2]);

%% create seeds for the centers of the charts and load into fitter

% % random
% n = 3;
% seeds = floor(rand([n 1])*size(mesh.v,1))+1;

% random
n = 5;
seeds = floor(rand([n 1])*size(mesh.v,1))+1;
newseeds = seeds;

% improve seeds using lloyd algorithm
nIter = 10;
for i = 1:nIter
    newseeds = perform_lloyd_mesh(mesh.v, mesh.f, newseeds);
end

%%

seedsX = [332 256 250; 250 350 230; 139 280 268];
% 250 256 332 ZYX
% 119 381 273 ZYX
% 268 280 139 
seeds = zeros([1 size(seedsX,1)]);

for i = 1:size(seedsX,1)
    [~,idx] = min(sqrt((mesh.v(:,1) - seedsX(i,1)).^2 + (mesh.v(:,2) - seedsX(i,2)).^2 + (mesh.v(:,3) - seedsX(i,3)).^2));
    seeds(i) = idx;
end

%%
% the fitter takes the mesh and partions it into overlapping submeshes
% (patches) based on the seeds with zero transitionWidth these patches
% would be geodesic Voronoi cells

fitOptions = struct('chartSeeds', seeds, 'transitionWidth', 0);

%fitOptions = struct('chartSeeds', seeds, 'transitionWidth', 100, 'fixAxis',1);
xp.setFitOptions(fitOptions);
xp.fitSurface(mesh);

noOverlap = xp.fitter.fittedParam;
%%
fitOptions = struct('chartSeeds', seeds, 'transitionWidth', 20,...
                                'diskSeeds', seeds(2), 'diskRadius', 105);

%fitOptions = struct('chartSeeds', seeds, 'transitionWidth', 100, 'fixAxis',1);
xp.setFitOptions(fitOptions);
xp.fitSurface(mesh);

%% visualize 

% whole mesh
xp.fitter.inspectMesh(1:4);

%%
% whole heartedly 

this = xp.fitter;

% color different regions of the surface
v = this.fittedParam.mesh.v;
C = 0*v(:,1);
for i = 1:numel(seeds)
    C(noOverlap.submVidx{i} & ~this.fittedParam.submVidx{4}) = i;
end
C(this.fittedParam.submVidx{4}) = 4;

trisurf(this.fittedParam.mesh.f, v(:,1), v(:,2), v(:,3), C, 'CDataMapping','direct')
shading interp
cmap = hsv(numel(this.fittedParam.submeshes));
colormap(cmap);
axis equal 

% add the chart boundaries 
hold on
for i = 1:numel(this.fittedParam.submeshes)
    v = this.fittedParam.submeshes{i}.v;
    bidx = this.fittedParam.submeshes{i}.b{1};
    plot3(v(bidx,1),v(bidx,2),v(bidx,3), '--', 'Color', cmap(i,:), 'LineWidth', 2);
end
hold off
camlight;
axis off;
viewdir = [1 0 1];
view(viewdir);

%%

% broken heart
cmap = hsv(numel(this.fittedParam.submeshes));
separation = 30;

% for radial separation, find CM and seed positions
PC = surfaceDetection.PointCloud(this.fittedParam.mesh.v);
[centroid, ~] = PC.getMoments();

seedIdx = xp.fitter.fitOptions.chartSeeds;
V = this.fittedParam.mesh.v;

% voronoi only colors
Cvor = 0*v(:,1);
for i = 1:numel(seeds)
    Cvor(noOverlap.submVidx{i}) = i;
end

clf;
hold on;

% the Voronoi cells
for i = 1:numel(seeds)

    fittedParam = xp.fitter.fittedParam; % noOverlap
    submVidx = fittedParam.submVidx{i};
    v = fittedParam.submeshes{i}.v;
    f = fittedParam.submeshes{i}.f;

    %d = V(seedIdx(i),:) - centroid;
    d = mean(fittedParam.submeshes{i}.vn);
    d = d./norm(d);
    d = d*separation;

    for j = 1:3
        v(:,j) = v(:,j) +  d(j);
    end    
    trisurf(f, v(:,1), v(:,2),v(:,3), 0*v(:,1) + i, 'CDataMapping','direct')
end

% the stitches
for i = 1:numel(seeds)

    v = xp.fitter.fittedParam.submeshes{i}.v;
    f = xp.fitter.fittedParam.submeshes{i}.f;

    d = mean(xp.fitter.fittedParam.submeshes{i}.vn);
    d = d./norm(d);
    d = d*separation;

    for j = 1:3
        v(:,j) = v(:,j) +  d(j);
    end    

    bidx = xp.fitter.fittedParam.submeshes{i}.b{1};
    plot3(v(bidx,1),v(bidx,2),v(bidx,3), '--', 'Color', cmap(i,:), 'LineWidth', 2);
end

% add the disk
i = 4;
seedIdx = xp.fitter.fitOptions.diskSeeds(i - numel(seeds));
submVidx = xp.fitter.fittedParam.submVidx{i};
v = xp.fitter.fittedParam.submeshes{i}.v;
f = xp.fitter.fittedParam.submeshes{i}.f;
%d = V(seedIdx,:) - centroid;
%d = xp.fitter.fittedParam.mesh.vn(seedIdx,:);
d = mean(xp.fitter.fittedParam.submeshes{4}.vn);
d = d./norm(d);
d = 3*d*separation;
for j = 1:3
    v(:,j) = v(:,j) +  d(j);
end    
trisurf(f, v(:,1), v(:,2),v(:,3), C(submVidx,:), 'CDataMapping','direct')


axis equal;
shading interp
colormap(cmap)
camlight
view(viewdir);
axis off;
hold off;

%%
% submesh
xp.fitter.inspectMesh();
view([0 0 1]);

%% Inspect the fit in a cross section
%
% inspect fit and point cloud over a cross section inthe  data. Dimensions 
% are x,y or z and the value has to be within the corresponding axis range.
%
zval = 100;
inspectOptions = struct('dimension','x','value',zval,'pointCloud','b', 'noalign', true);
xp.fitter.inspectQuality(inspectOptions, xp.detector, xp.stack);

%% normally evolve

shift = 2;
xp.normallyEvolve(shift);

%%
xp.fitter.setDesiredChart('conformal', true);
xp.fitter.setDesiredChart('exponential', false);
xp.fitter.setFitDomain(xp.stack.image.domain);
xp.generateSOI();

%% Pullback the stack to the desired charts
%
% Pass the region of interest and the current time to pull back the stack
% in the desired charts. This generates the data fields containing the
% pullback.
%
% 21, 2, 20
onionOpts = struct('nLayers', 3, 'layerDistance', 2, 'sigma', 1, 'makeMIP', true);
xp.SOI.pullbackStack(xp.stack, [], xp.currentTime, onionOpts);

%%
% Now we extract the data field from the surface of interest at the current
% time, which is the time of the fit. 
%
data = xp.SOI.getField('data_MIP');
data = data(xp.tIdx(xp.currentTime));

for i = 1:4;

    addpath('/Users/idse/Dropbox/Ken/lattice_analysis/lattice_analysis_Aug28_2014/functions/')
    subplot_tight(2,2,i);
    
    %type = 'cylinder';
    %type = 'exponential';
    type = 'conformal';
    patchName     = [type '_' num2str(i) '_index'];
    transformName = [type '_' num2str(i)];
    pb = data.getPatch(patchName).getTransform(transformName).apply{1};
    pb(pb == 0) = max(pb(:));
    imshow(imadjust(pb'),[],'InitialMagnification',66)

    % the stitches
    tidx = xp.tIdx(xp.currentTime);
    u = xp.fitter.fittedParam.submeshes{i}.u{1};
    f = xp.fitter.fittedParam.submeshes{i}.f;
    bidx = xp.fitter.fittedParam.submeshes{i}.b{1};
    curChart = xp.SOI.atlas(tidx).getChart(transformName);
    ub = u(bidx,1);
    vb = u(bidx,2);
    umin = curChart.image.boundary{1}(1);
    vmin = curChart.image.boundary{2}(1);
    ub = (ub-umin)./curChart.image.stepSize(1);
    vb = (vb-vmin)./curChart.image.stepSize(2);
    hold on
    plot(vb,ub, '--', 'Color', cmap(i,:), 'LineWidth', 2);
    hold off
end

%%
% heart with texture map

viewdir = ([1 1 1]);
emb = xp.SOI.embedding(xp.tIdx(xp.currentTime));
data = xp.SOI.data(xp.tIdx(xp.currentTime));

separation = 50*[1 1 1 2];

% for radial separation, find CM and seed positions
seedIdx = xp.fitter.fitOptions.chartSeeds;
V = xp.fitter.fittedParam.mesh.v; 
PC = surfaceDetection.PointCloud(V);
[centroid, ~] = PC.getMoments();

figure
hold on;

MIP = xp.SOI.getField('data_MIP');

% texture 
C = {};
Cc = {};
for i = 1:numel(seedIdx)+1
    C{i} = MIP(31).patches{i}.apply{1};
    if i == 1
        Cmin = min(C{i}(:));
        Cmax = max(C{i}(:));
    end
    C{i} = mat2gray(C{i}, double([Cmin 0.5*Cmax]));
    
    % make it color
    Cc{i} = zeros([size(C{i}) 3]);
    for c = 1:3
        Cc{i}(:,:,c) = cmap(i,c)*C{i};
    end
end

% the pieces
for i = 1:numel(seedIdx)+1
    X = emb.patches{i}.apply;
    d = mean(xp.fitter.fittedParam.submeshes{i}.vn);
    d = d./norm(d);
    d = separation(i)*d;
    surf(X{1} + d(1),X{2} + d(2),X{3} + d(3), C{i}, 'FaceColor','texturemap');
end

axis equal;
shading flat
colormap gray
%camlight
view(viewdir);
axis off;

% % seeds
% for i = 1:numel(seedIdx)
%     scatter3(V(seedIdx(i),1), V(seedIdx(i),2), V(seedIdx(i),3), '.r', 'filled');
% end

% the stitches
for i = 1:numel(seeds)+1

    v = xp.fitter.fittedParam.submeshes{i}.v;
    f = xp.fitter.fittedParam.submeshes{i}.f;

    d = mean(xp.fitter.fittedParam.submeshes{i}.vn);
    d = d./norm(d);
    d = d*separation(i);

    for j = 1:3
        v(:,j) = v(:,j) +  d(j);
    end    

    bidx = xp.fitter.fittedParam.submeshes{i}.b{1}(1:5:end);
    plot3(v(bidx,1),v(bidx,2),v(bidx,3), '--', 'Color', cmap(i,:), 'LineWidth', 2);
end

hold off;

%% Save the surface of interest to disc
%
% Here we save the SOI using SOI.save. We set the following options:
%
% * dir:            The directory to save the SOI to.
% * imwriteOptions: Pullbacks are saved to image files using imwrite, we
% can pass options to change file format, compression etc. For example we
% could change this option to
% imwriteOptions = {'jp2', 'Mode', 'lossless'}; 
% * make8bit:       Often absolute intensities don't matter and 8 bit offers
% a large enough dynamic range. This options rescales the lookup table and
% converts to 8 bit before saving.
%
imwriteOptions = {'tif'};
saveDir = fullfile(projectDir, 'heart');

options = struct('dir',saveDir,'imwriteOptions',{imwriteOptions},...
                    'make8bit',false);
xp.SOI.save(options)


