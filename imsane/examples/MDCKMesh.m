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
clear all; %close all;
%
% Setup a working directory for the project, where extracted surfaces,
% metadata and debugging output will be stored. Also specify the directory
% containing the data. 
%
[scriptPath,~,ext] = fileparts(matlab.desktop.editor.getActiveFilename);

%dataDir = fullfile(scriptPath, 'rawData');
dataDir = '/Users/idse/Dropbox/MDCK';
projectDir = fullfile(dataDir, 'projectFiles');
%projectDir = cd;

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
fileMeta.filenameFormat  = '92_ctrlCK666_642_48xz9CropScale2.tif';
fileMeta.timePoints      = [60]; % for full data sample use 0;
fileMeta.stackResolution = [1 1 1];%[.5 .5 .5]; 
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
expMeta.description      = 'MDCK cyst';
expMeta.dynamicSurface   = 1;
expMeta.jitterCorrection = 0; % 1: Correct for sample translation
expMeta.fitTime          = fileMeta.timePoints(1); 
expMeta.detectorType     = 'surfaceDetection.fastCylinderDetector';%
expMeta.fitterType       = 'surfaceFitting.cylinderMeshWrapper';%meshWrapper'; 

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

% FOR MESH RELATED SCRIPTS -> SHOULD GO IN FILEMETA
[~,barefName,~] = fileparts(fileMeta.filenameFormat);
PCfile = fullfile(projectDir, ['pointCloud_' barefName '.obj']);

%% Load data for surface detection and rescale to unit aspect ratio
%
% Load a timepoint from the data.
% loadTime sets currentTime, loads the stack and resets detector and fitter
% with the appropriate options for that time. Optionally rescale the stack 
% to unit aspect ratio if desired.
%
xp.loadTime(fileMeta.timePoints(1));

%%
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
OBJ.vertices = PC.unalignedPoints;
OBJ.objects(1).type='f';
OBJ.objects(1).data.vertices=[];
write_wobj(OBJ, PCfile);

%% Copy meshlab command to clipboard

% if meshlabserver is in the bash path, just paste in the terminal and
% press enter

outputMesh = fullfile(projectDir, ['poissonMesh_' barefName '.obj']);
meshlabScript = fullfile(scriptPath, 'makePoissonRefined.mlx');
clipboard('copy', ['meshlabserver -i ' PCfile ' -o ' outputMesh ' -s ' meshlabScript ' -om vn']);

%% Read surface mesh produced by meshlab
%outputMesh = fullfile(projectDir, 'obj','heartPC95_mesh.obj');

tic
mesh = readObj(outputMesh);
mesh.f = mesh.f.v;
toc

%% create seeds for the centers of the charts and load into fitter

% % random
% n = 3;
% seeds = floor(rand([n 1])*size(mesh.v,1))+1;

% the tips
zdir = 1;
% seeds(1) = find(mesh.v(:,zdir) == min(mesh.v(:,zdir)));
% seeds(2) = find(mesh.v(:,zdir) == max(mesh.v(:,zdir)));
seeds = [];

% the fitter takes the mesh and partions it into overlapping submeshes
% (patches) based on the seeds with zero transitionWidth these patches
% would be geodesic Voronoi cells

fitOptions = struct('chartSeeds', seeds, 'transitionWidth', 100);
xp.setFitOptions(fitOptions);
xp.fitSurface(mesh);

%xp.fitter.smoothMesh('distance', 100);
shift = -45;
xp.normallyEvolve(shift);
%xp.fitter.smoothMesh('distance', 10);

%% visualize 

% whole mesh
xp.fitter.inspectMesh();

%%
% submesh
xp.fitter.inspectMesh(1);
view([0 0 1]);

%% Inspect the fit in a cross section
%
% inspect fit and point cloud over a cross section inthe  data. Dimensions 
% are x,y or z and the value has to be within the corresponding axis range.
%
zval = 550;
inspectOptions = struct('dimension','x','value',zval,'pointCloud','b', 'noalign', true);
figure, xp.fitter.inspectQuality(inspectOptions, [], xp.stack);

%%
xp.fitter.setDesiredChart('cylinder1', true);
xp.fitter.setDesiredChart('cylinder2', true);
xp.fitter.setDesiredChart('cylinder1_proper', true);
xp.fitter.setDesiredChart('cylinder2_proper', true);
xp.generateSOI();

% Pullback the stack to the desired charts
%
% Pass the region of interest and the current time to pull back the stack
% in the desired charts. This generates the data fields containing the
% pullback.
%
%%

onionOpts = struct('nLayers', 1, 'layerDistance', 3, 'sigma', 5);
xp.SOI.pullbackStack(xp.stack, [], xp.currentTime, onionOpts);

%%
% Now we extract the data field from the surface of interest at the current
% time, which is the time of the fit.
%
data = xp.SOI.getField('data');
data = data(xp.tIdx(xp.currentTime));

patchName     = 'cylinder1_index';%[type '_' num2str(i) '_index'];
transformName = 'cylinder1';%[type '_' num2str(i)];
%transformName = 'cylinder2';%[type '_' num2str(i)];

pb = data.getPatch(patchName).getTransform(transformName).apply{1};
pb = mat2gray(pb', double([min(pb(pb>0)) 0.6*max(pb(:))]));

%pb = mat2gray(adapthisteq(pb));
figure, imshow(pb,'InitialMagnification',66)

%% make an onion stack

patchName     = 'cylinder2_index';%[type '_' num2str(i) '_index'];
transformName = 'cylinder2';%_proper';%[type '_' num2str(i)];

nLayers = onionOpts.nLayers;
halfLayers = (nLayers - 1)/2;

onion = zeros([data.getPatch(patchName).getTransform(transformName).domain.gridSize nLayers],'uint16');

for li = 1:nLayers

    idx = li - halfLayers - 1;

    % convert index to field name
    if idx < 0
        fieldName = ['data_layer_m' num2str(-idx)];
    elseif idx > 0
        fieldName = ['data_layer_p' num2str(idx)];
    else
        fieldName = 'data';
    end


    data = xp.SOI.getField(fieldName);
    data = data(xp.tIdx(xp.currentTime));

    pb = data.getPatch(patchName).getTransform(transformName).apply{1};
    pb = mat2gray(pb', double([min(pb(pb>0)) 0.6*max(pb(:))]));

    onion(:,:,li) = pb;

end

%% write to file

imwrite(onion(:,:,1), 'MDCKonion2.tif')
for li = 2:nLayers
    imwrite(onion(:,:,li), 'MDCKonion2.tif', 'writemode', 'append');
end

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
saveDir = fullfile(projectDir, 'heart');

options = struct('dir',saveDir, 'make8bit',false, 'saveMetric', false);
xp.SOI.save(options)

%% 

SOIdir = saveDir;
loadedSOI = surfaceAnalysis.SurfaceOfInterest(SOIdir);

%%
loadedSOI.NCalcInducedMetric(xp.currentTime);

%% curvature

addpath(genpath('/Users/idse/Dropbox/heart_mesh/toolbox_graph/'));

options = struct();
options.curvature_smoothing = 5;
[Umin,Umax,Cmin,Cmax,Cmean,Cgauss,Normal] = compute_curvature(mesh.v, mesh.f, options);
1/max(Cmax)

subm = xp.fitter.fittedParam.mesh;
fv.Vertices = [subm.v(:,1), subm.v(:,2), subm.v(:,3)];
fv.Faces = subm.f;
fv.FaceVertexCData = [1 1 1];

trisurf(subm.f, subm.v(:,1), subm.v(:,2), subm.v(:,3), Cmin)
axis equal;
colormap hsv;
view([1 0 0]);

%% mesh fairing
% try heat kernel

subm = xp.fitter.fittedParam.mesh;

laplacian_type = 'distance';
%laplacian_type = 'conformal';

vertex = subm.v';
face = subm.f;

% Laplacian
if not(strcmp(laplacian_type,'conformal'))
    options.symmetrize = 1;
    options.normalize = 0;
    L = compute_mesh_laplacian(vertex,face,laplacian_type,options);
else
    options.symmetrize = 0;
    options.normalize = 1;
    L = compute_mesh_laplacian(vertex,face,laplacian_type,options);
end

%% heat diffusion flow

Tlist = [0 10 40 500];
options.dt = 0.3;
clf;
for i=1:length(Tlist)
    options.Tmax = Tlist(i);
    vertex1 = perform_mesh_heat_diffusion(vertex,face,L,options);
    % display
    subplot(1,length(Tlist),i);
    plot_mesh(vertex1,face,options);
    shading interp; camlight; axis tight;
end
saveas(gcf, [rep name '-smoothing-' laplacian_type '.png'], 'png');
