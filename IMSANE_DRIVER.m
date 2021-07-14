%% ImSAnE Eye Disc Tutorial
%
% NOTE: prior to using imsane, you must install it by running setup.m,
% which is located within the imsane directory
%
% This an example of how to detect and fit what we call planar surfaces,
% where the surface can be described as a height for every x,y.
%
% Note that ImSAnE is fully documented and additional information about
% available properties, methods and options can be found using the matlab
% documentation system. 
% 
% For example type:
doc surfaceDetection.planarDetector

% This is also useful to double checking that imsane is installed correctly
% and within your matlab path. If the documentation does not pop up, this
% means that matlab is not finding the imsane files.

%% Initialize the project
%
% We start by creating an experiment object, which holds this metadata and 
% provides a frontend for a number of tasks such as loading the raw data.
% To construct the experiment object we need to pass dataDir, the directory 
% containing the raw data and projectDir, the directory where results of 
% the script will be saved.
%
% NOTE: data should be saved as individual z-stacks for each time point.
% Provide the path to the directory containing these files below.

clear all; close all;

[scriptPath,~,ext] = fileparts(matlab.desktop.editor.getActiveFilename);

dataDir = '/Users/kevin/Documents/ImagingData/WT1/surface_detect_prep/';
projectDir = fullfile(scriptPath, 'projectFiles');

% test
xp = project.Experiment(projectDir, dataDir);

%%
% Next we set the metadata pertaining to the raw data files in the structure
% fileMeta. ImSAnE assumes that timepoints are saved individually and that 
% filenames in a timeseries are identical up to an integer specifying the 
% timepoint. Therefore we have:
%
% * filenameFormat:               Filename, %u in the position of the integer
% * timePoints = [t1, t2, ..] :   List of times available. In this example 
%                                 we have just a single time point 0.
% * stackResolution :             Stack resolution in micron.

fileMeta = struct();
fileMeta.dataDir = dataDir;
fileMeta.filenameFormat = 'WT1_fullStack_T%d.tif';
fileMeta.swapZT = 0;
fileMeta.timePoints = 0:119;
fileMeta.stackResolution = [0.0505 0.0505 0.5];

%% 
% In the structure expMeta we set general parameters for the surface
% analysis we will do. 
%
% * channelsUsed:       Which channels do we need.
% * channelColor:       Assign color to the channels, RGB = [1 2 3].
%                       In this example the first channel is E-cadherin
% * dynamicSurface:     Does the surface shape change with time?  
%                       For a single time point this is false. True is not yet
%                       supported.
% * jitterCorrection:   Not needed here.
% * detectorType:       Which type of surface detector will be used.
%                       We will look at only one side of the eye, which is
%                       a planar surface so we use MIPDetector, which is 
%                       one of the planarDetector methods.
% * fitterType:         Which type of fitter will be used.
%                       We fit planar surfaces using Thin Plate Spline:
%                       tpsFitter.

expMeta = struct();
expMeta.description = 'eye disc, resonance scanner, Ecad';
expMeta.channelsUsed = [1];
expMeta.channelColor = [1];
expMeta.dynamicSurface = 1;
expMeta.jitterCorrection = 0;
expMeta.fitTime = fileMeta.timePoints(120);     % total number of time points
expMeta.detectorType = 'surfaceDetection.MIPDetector';
expMeta.fitterType = 'surfaceFitting.tpsFitter';

xp.setFileMeta(fileMeta);
xp.setExpMeta(expMeta);

%%
% Finally we call initNew(), which reads the stack size from the first 
% available time point, then initializes fitter and detector and creates 
% fitOptions and detectOptions based on their defaults.

xp.initNew();

%% Load a time point from the data
%
% Now that we have set up the project, we can load a time point.
% loadTime sets xp.currentTime, loads the stack into xp.stack 
% and resets detector and fitter with the default options for that time.
%
% We are loading in just one time point for now in order to determine the
% best parameters for finding the surface. Once we figure this out, we will
% batch process the remaining time points
%
% We rescale to unit aspect ratio to make detection work better later and
% visualize in the right proportions.

xp.loadTime(5);
xp.rescaleStackToUnitAspect();

%% check the image content; 

% display MIP of the loaded stack
 im = xp.stack.image.apply();
 imshow(max(im{1},[],3),[]);


%% look at single slice

% xp.stack is not an array but a Stack object.
% The easy way to look at a slice through the data is using getSlice.

imshow(xp.stack.getSlice('z', 110), []);

%% Detect the surface
%
% MIPDetector.detectSurface detects the surface as the position of the 
% brightest pixel in z for every xy coordinate
%
% A number of detection options directly affect detection:
%
% * sigma :     Width of the Gaussian z-derivative.
% * channels :  Channels (summed) to use for detection.
% * zdir :      Dimension corresponding to z, minus flips direction.
% Flipping the direction can sometimes improve detection.
%
% Then there are options which filter the result and can be modified
% without redetecting:
%
% * maxIthresh:     Throw out points with MIP dimmer than this.
% * summedIthresh:  Throw out points with SIP dimmer than this.
% * sigZoutliers:   Remove height outliers after all other masks.
% * scaleZoutliers: Spatial scale of outlier removal.
%
% scaleZoutliers is the linear size of a region over which the
% distribution of height is computed, sigZoutliers is then a cutoff in
% units of standard deviation of this distribution to remove misdetected
% points far above or below the other points in the region.

%customOptions = struct('dilSize', 0, 'erSize', 20,'areaOpenSize', 1000, ...
%                        'fillSize', 500);

detectOptions = struct('sigma', 5, 'channels', 1, 'zdir', -3,...
                        'maxIthresh', 0.05, 'summedIthresh', 0.05,...
                        'sigZoutliers', 1, 'scaleZoutliers', 50,...
                        'seedDistance', 20); %,'customOptions', customOptions); 

% Calling detectSurface runs the surface detector and creates the point
% cloud in detector.pointCloud.

xp.setDetectOptions(detectOptions);
xp.detectSurface();


%% visualize height map of detected surface

% Different from the other detectors, the detected surface is
% represented not only by a PointCloud object but also by an image
% surfaceMatrix, containing z values for each xy.
% Looking at this height map masked by the filters specified in
% detectOptions one can judge how well the surface was detected.

imshow(xp.detector.mask.*xp.detector.surfaceMatrix, [],...
                                            'InitialMagnification', 100);


%% visualize cross section from point cloud of detected surface

% We can also inspect a point cloud cross section over the data with
% detector.inspectQuality. In the pointCloud option, 'c' specifies the 
% color cyan.

inspectOptions= struct('dimension', 'y', 'value', 100, 'pointCloud', 'c');
xp.detector.inspectQuality(inspectOptions, xp.stack);

%% visualize point cloud of detected surface in 3D

% Or we can look at the point cloud in 3d, with some subsampling factor.
ssfactor = 200;
xp.detector.pointCloud.inspect(ssfactor);

%% Fit the surface for the disc proper cells
%
% By detecting the brightest pixel in z for each x,y coordinate in the
% E-cad channel and filtering out local outliers we have found the apical
% surface of the disc proper cells. We can now fit a smooth surface
% representation to that.
%
% tpsFitter fits the pointcloud using a thin plate spline fit. It has the
% following options:
%
% * gridSize:     Size of grid on which to generate fitted surface
%               default [50 50], full size takes long.
% * smoothing:    TPS smoothing parameter (default 1000).

fitmask = imfill(imdilate(xp.detector.mask, strel('disk',20)), 'holes');
fitmask = imerode(fitmask,strel('disk',10));

fitOptions = struct('smoothing', 1000, 'gridSize', [80 80], 'fitMask', fitmask);
fitOptions.shift = -10;     % TEST DIFFERENT VALUES FOR YOUR OWN DATA 
shift = xp.fitOptions().shift;

xp.setFitOptions(fitOptions);
xp.fitSurface();

%% visualize cross section of surface spline

% We can visualize the result on a cross section with
% fitter.inspectQuality.

inspectOptions= struct('dimension', 'y', 'value', 500, 'pointCloud', 'c');%1150
xp.fitter.inspectQuality(inspectOptions, xp.detector, xp.stack);

%% visualize spline surface in 3D

% visualize 3D
xp.fitter.inspectTPS()

%% genereate charts

% We now generate the Surface Of Interest. The charts to be generated are 
% specified in xp.fitter.charts. In this case there is only one, called
% 'xy'. 

xp.generateSOI();

%% Pull back the data to the surface
% 
% We pull back the data to the SOI using pullbackStack.

xp.SOI.pullbackStack(xp.stack, xp.currentROI, xp.currentTime);

%%
% To look at the pullback, we call the data field of the SOI at the right
% time, and get a particular patch from that with getPatch. A patch is a 
% part of a surface. In this case, there is only one called xy_index.
% Then we get the data in some patch in a particular coordinate system with
% getTransform. In this case there is only one coordinate system: xy.
% What we get is an object not only holding the image data but also
% metadata and methods to manipulate it. The actual data is obtained by
% calling the method apply. This returns a cell array with entries for each
% channel.

% xp.tIdx converts the time into an index in a list of time points
tidx = xp.tIdx(xp.currentTime);

% the second channel is Ecad
channel = 1;

discProperPatch = xp.SOI.data(tidx).getPatch('xy_index');

discProperImage = discProperPatch.getTransform('xy').apply{channel};

figure, imshow(discProperImage, xp.stack.Ilim{1}, 'InitialMagnification', 50);

%% Onion options

% Onion options (onionOpts) allows you extract a series of layers of defined
% distance from the surface we have found and create an MIP of these. This is
% useful for capturing all of the signal from the surface in regions where
% the surface may not have been fit the best. If you run this, when you
% save the data, the real surface, the MIP, and all the layers will be
% saved in different folders.

onionOpts = struct('nLayers', 9, 'layerDistance', 1, 'sigma', 1,'makeMIP','MIP');
xp.SOI.pullbackStack(xp.stack, xp.currentROI, xp.currentTime, onionOpts);


%% compare MIP and pullback

MIP = max(xp.stack.image.apply{1},[],3);
blaMIP = mat2gray(MIP);
blaSOI = mat2gray(discProperImage,double([min(MIP(:)) max(MIP(:))]));
color = cat(3, blaSOI, blaMIP, 0*blaMIP);
figure, imshow(color)

imwrite(blaSOI, fullfile(dataDir, 'analysis', 'SOI.tif'));
imwrite(blaMIP, fullfile(dataDir, 'analysis', 'MIP.tif'));
imwrite(color, fullfile(dataDir, 'analysis', 'SOIMIP.tif'));

%%
% All metadata is saved in SOI.xml. Pullbacks, geometry and any other data
% defined on the surface are saved to image files in subdirectories 
% reflecting the structure of patches and coordinates explained earlier in 
% this tutorial. We can reload a surface of interest with
% SOI.load(directory)

%% Batch proces all other times
% first set the fit options for all times to what we determined above

allTimes = true;

xp.setDetectOptions(detectOptions, allTimes);

fitOptions.shift = -10;
xp.setFitOptions(fitOptions, allTimes);

% then run batchProcess

timePoints = 0:119;
batchOptions = struct('unitAspect', true, 'seeded', false);
%seed = [];

%xp.generateSOI();

for t = 1:length(timePoints)

    debugMsg(1, ['Processing time point ' num2str(timePoints(t)) '\n']);

    % load stack and rescale
    xp.loadTime(timePoints(t));
    if batchOptions.unitAspect
        xp.rescaleStackToUnitAspect();
    end

    % surface detection
    xp.setDetectOptions(xp.detectOptions(t));
    xp.detectSurface();

    % surface fitting
    fitOptions = xp.fitOptions(t);
    fitmask = imfill(imdilate(xp.detector.mask, strel('disk',20)), 'holes');
    fitmask = imerode(fitmask,strel('disk',10));
    fitOptions.fitMask = fitmask;

    xp.setFitOptions(fitOptions);
    xp.fitSurface();

    shift = xp.fitOptions(t).shift;
    if shift > 0
        xp.zEvolve(shift);
    end
    
    % SOI and pullbacks
    xp.fitter.populateSOI(xp.SOI, xp.currentTime);
    xp.SOI.pullbackStack(xp.stack, xp.currentROI, xp.currentTime);
    
    onionOpts = struct('nLayers', 9, 'layerDistance', 1, 'sigma', 5,'makeMIP','MIP');
    xp.SOI.pullbackStack(xp.stack, xp.currentROI, xp.currentTime, onionOpts);

    if batchOptions.seeded 
        seed = xp.fitter.generateSeed();
    end 
end


%% Save the result
%
% Finally we save the SOI using SOI.save. We set the following options:
%
% * dir:            The directory to save the SOI to.
% * imwriteOptions: Pullbacks are saved to image files using imwrite, we
%                   can pass options to change file format, compression etc.
%                   For example we could change this option to
%                   imwriteOptions = {'jp2', 'Mode', 'lossless'}; 
% * make8bit:       Often absolute intensities don't matter and 8 bit offers
%                   a large enough dynamic range. This options rescales the 
%                   lookup table and converts to 8 bit before saving.

imwriteOptions = {'tif'};
savedir = fullfile(scriptPath, 'discProperApicalSOI');

options = struct(   'dir',              savedir,...
                    'imwriteOptions',   {imwriteOptions},...
                    'make8bit',         true);
xp.SOI.save(options)

