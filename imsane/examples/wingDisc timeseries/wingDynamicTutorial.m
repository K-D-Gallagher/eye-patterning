%% ImSAnE Wing Disc Tutorial
%
% In this tutorial we detect the apical surface of the disc proper cells
% and the peripodial cells on one side of a wing disc two hours after
% puparium formation. We read out actin and E-cadherin on these surfaces.
%
% This an example of how to detect and fit what we call planar surfaces,
% where the surface can be described as a height for every x,y.
% It also demonstrates using multiple channels and making multiple
% surfaces. Finally, because finding folds as regions of high curvature is
% part of detecting the peripodial surface, this example shows how to
% calculate the surface metric and curvature.
%
% Note that ImSAnE is fully documented and additional information about
% available properties, methods and options can be found using the matlab
% documentation system. 
% 
% For example type:

%%
doc surfaceDetection.planarDetector

%% Initialize the project
%
% We start by creating an experiment object, which holds this metadata and 
% provides a frontend for a number of tasks such as loading the raw data.
% To construct the experiment object we need to pass dataDir, the directory 
% containing the raw data and projectDir, the directory where results of 
% the script will be saved.

clear all; close all;

[scriptPath,~,ext] = fileparts(matlab.desktop.editor.getActiveFilename);

%dataDir = fullfile(scriptPath, 'rawData');
dataDir = '/Users/idse/BIG data/Ken/20140501Yki';
projectDir = fullfile(scriptPath, 'projectFiles');

% test
xp = project.Experiment(projectDir, dataDir);

%%
% Next we set the metadata pertaining to the raw data files in the structure
% fileMeta. ImSAnE assumes that timepoints are saved individually and that 
% filenames in a timeseries are identical up to an integer specifying the 
% timepoint. Therefore we have
%
% * filenameFormat:               Filename, %u in the position of the integer
% * timePoints = [t1, t2, ..] :   List of times available. In this example 
% we have just a single time point 0.
% * stackResolution :             Stack resolution in micron.

fileMeta = struct();
fileMeta.dataDir = dataDir;
fileMeta.filenameFormat = '20140501_Live_EcadGFP_AyGal_B_YkiS168A_T%u.ome.tif';
fileMeta.timePoints = 0:160;
fileMeta.swapZT = 0;
fileMeta.stackResolution = [0.1 0.1 0.7]; 

%% 
% In the structure expMeta we set general parameters for the surface
% analysis we will do. 
%
% * channelsUsed:   Which channels do we need.
% * channelColor:   Assign color to the channels, RGB = [1 2 3].
%                   In this example the first channel is E-cadherin, and 
%                   the second is actin. We want these in green and red,
%                   respectively.
% * dynamicSurface: Does the surface shape change with time?  
%                   For a single time point this is false. True is not yet
%                   supported.
% * jitterCorrection:   Not needed here.
% * detectorType:       Which type of surface detector will be used.
%                       We will look at only one side of the wing, which is
%                       a planar surface so we use planarDetector.
% * fitterType:         Which type of fitter will be used.
%                       We fit planar surfaces using Thin Plate Spline:
%                       tpsFitter.

expMeta = struct();
expMeta.description = 'Wing with Yki clones, Ecad';
expMeta.channelsUsed = [1 2];
expMeta.channelColor = [1 2];
expMeta.dynamicSurface = true;
expMeta.jitterCorrection = false;
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
% We rescale to unit aspect ratio to make detection work better later and
% visualize in the right proportions.

xp.loadTime(0);
xp.rescaleStackToUnitAspect();

%% 
% xp.stack is not an array but a Stack object.
% The easy way to look at a slice through the data is using getSlice.

imshow(xp.stack.getSlice('z', 70), []);

%% Detect the surface
%
% planarDetector.detectSurface detects the surface as the position of the 
% maximal Gaussian z-derivative in some direction, i.e. the position of the
% largest intensity jump along some direction and smoothened over some
% scale.
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

customOptions = struct('dilSize', 6, 'erSize', 20, 'areaOpenSize', 1000,...
                        'fillSize', 1000);

detectOptions = struct(  'sigma', 2, 'channels', 2, 'zdir', 3,...
                        'maxIthresh', 0.1, 'summedIthresh', 0,...
                        'sigZoutliers', 1, 'scaleZoutliers', 50,...
                        'seedDistance', 20, 'customOptions', customOptions); 

% Calling detectSurface runs the surface detector and creates the point
% cloud in detector.pointCloud.

xp.setDetectOptions(detectOptions);
xp.detectSurface();

% Different from the other detectors, the detected surface is
% represented not only by a PointCloud object but also by an image
% surfaceMatrix, containing z values for each xy.
% Looking at this height map masked by the filters specified in
% detectOptions one can judge how well the surface was detected.

imshow(xp.detector.mask.*xp.detector.surfaceMatrix, [],...
                                            'InitialMagnification', 40);

%% 
% One can then find better filter parameters without redetecting the
% surface by changing the second block of options in detectOptions and 
% calling resetMask and applyMasks. 

xp.detector.resetMask();

customOptions = struct('dilSize', 50, 'erSize', 30, 'areaOpenSize', 5000,...
                        'fillSize', 1000);

customOptions = struct('dilSize', 0, 'erSize', 0, 'areaOpenSize', 5000,...
                        'fillSize', 1000);

detectOptions = struct(  'sigma', 2, 'channels', 2, 'zdir', 3,...
                        'maxIthresh', 0.1, 'summedIthresh', 0,...
                        'sigZoutliers', 1, 'scaleZoutliers', 50,...
                        'seedDistance', 20, 'customOptions', customOptions); 

xp.detector.setOptions(detectOptions);    
xp.detector.applyMasks();

%%
imshow(xp.detector.mask.*xp.detector.surfaceMatrix, [],...
                                            'InitialMagnification', 40);
                                        
%%

%imshow(mat2gray(xp.detector.surfaceMatrix))

bla = imextendedmax(xp.detector.surfaceMatrix, 10);
imshow(bla, [])
%%
imshow(xp.detector.mask.*~bla.*xp.detector.surfaceMatrix, [])

%%
imshow(cat(3, mat2gray(xp.detector.surfaceMatrix), imextendedmax(xp.detector.surfaceMatrix, 10),0*imregionalmax(xp.detector.surfaceMatrix)), [])

%%
% We can also inspect a point cloud cross section over the data with
% detector.inspectQuality. In the pointCloud option, 'c' specifies the 
% color cyan.

inspectOptions= struct('dimension', 'y', 'value', 904, 'pointCloud', 'c');
xp.detector.inspectQuality(inspectOptions, xp.stack);

%% 
% Or we can look at the point cloud in 3d, with some subsampling factor.
ssfactor = 50;
xp.detector.pointCloud.inspect(ssfactor);

%% Fit the surface for the disc proper cells
%
% By detecting the largest intensity jump along z for each x,y in the
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

fitOptions = struct('smoothing', 100, 'gridSize', [80 80], 'fitMask', fitmask);

xp.setFitOptions(fitOptions);
xp.fitSurface();

%%
% We can visualize the result on a cross section with
% fitter.inspectQuality.

inspectOptions= struct('dimension', 'y', 'value', 400, 'pointCloud', 'c');%1150
xp.fitter.inspectQuality(inspectOptions, xp.detector, xp.stack);

%%
% visualize 3D
xp.fitter.inspectTPS()

%%
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
channel = 2;

discProperPatch = xp.SOI.data(tidx).getPatch('xy_index');

discProperImage = discProperPatch.getTransform('xy').apply{channel};

figure, imshow(discProperImage, xp.stack.Ilim{1}, 'InitialMagnification', 50);

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

fitOptions.shift = 0;
xp.setFitOptions(fitOptions, allTimes);

% then run batchProcess

timePoints = 11:50;
batchOptions = struct('unitAspect', true, 'seeded', true);
seed = [];

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
    xp.detectSurface(seed);

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

    if batchOptions.seeded 
        seed = xp.fitter.generateSeed();
    end 
end  

%%
t = 0;
% xp.tIdx converts the time into an index in a list of time points
tidx = xp.tIdx(t);

% the first channel is Ecad
channel = 2;

discProperPatch = xp.SOI.data(tidx).getPatch('xy_index');
discProperImage = discProperPatch.getTransform('xy').apply{channel};
imshow(discProperImage, xp.stack.Ilim{1}, 'InitialMagnification', 50);

%% Save the result
%
% Finally we save the SOI using SOI.save. We set the following options:
%
% * dir:            The directory to save the SOI to.
% * imwriteOptions: Pullbacks are saved to image files using imwrite, we
% can pass options to change file format, compression etc. For example we
% could change this option to
% imwriteOptions = {'jp2', 'Mode', 'lossless'}; 
% * make8bit:       Often absolute intensities don't matter and 8 bit offers
% a large enough dynamic range. This options rescales the lookup table and
% converts to 8 bit before saving.

imwriteOptions = {'tif'};
savedir = fullfile(scriptPath, 'discProperApicalSOI');

options = struct(   'dir',              savedir,...
                    'imwriteOptions',   {imwriteOptions},...
                    'make8bit',         true);
xp.SOI.save(options)


%% make a video of the shape

clear M;

nTimePts = 5;

minFrameH = 10^10;
minFrameW = 10^10;


for t = timePoints % 0:nTimePts;
    % xp.tIdx converts the time into an index in a list of time points
    tidx = xp.tIdx(t);

    grids = xp.SOI.embedding(tidx).getPatch('xy_index').apply;
    
    discProperPatch = xp.SOI.data(tidx).getPatch('xy_index');
    
    Ecad = discProperPatch.getTransform('xy').apply{2};
    Ilim = xp.stack.Ilim{2};
    Ilim(2) = (Ilim(2)-Ilim(1))*1.5 + Ilim(1);
    Ecad = mat2gray(Ecad,Ilim);
    
    clones = mat2gray(discProperPatch.getTransform('xy').apply{1}, xp.stack.Ilim{1});
    skin = cat(3, 0*Ecad', Ecad', clones');
    
    %surf(double(grids{1})',double(grids{2})',double(grids{3})', skin, 'FaceColor', 'texture');
    surf(double(grids{3})', skin, 'FaceColor', 'texture');
    
    %surf(double(grids{3})');
    %surf(xp.stack.imageSize(3)-double(grids{3}));

    view([1 0 1])
    whitebg([0 0 0]);
    axis off;
    shading interp;
    colormap gray;
    axis equal
    
    M(tidx) = getframe;
    
    frameH = size(M(tidx).cdata, 1);
    frameW = size(M(tidx).cdata, 2);
    
    if frameH < minFrameH,  minFrameH = frameH; end
    if frameW < minFrameW,  minFrameW = frameW; end
end

%%
writerObj = VideoWriter(fullfile(dataDir, '3DshapeTexture.avi'));
writerObj.Quality = 100;
writerObj.FrameRate = 2; % lower frame is slower video
open(writerObj);

for t = timePoints
    
    tidx = xp.tIdx(t);
    % make the frame all the same size 
    M(tidx).cdata = M(tidx).cdata(1:minFrameH, 1:minFrameW, :);
    
    % write to video
    writeVideo(writerObj, M(tidx));
end
    
close(writerObj);
