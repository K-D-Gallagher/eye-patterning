
% The goal of this driver file is to segment and track images. Unless you
% find a way to achieve perfect pixel classification, this will
% unfortunately involve some manual correction. This driver file brings you
% from initial pixel classification using an external program (some options
% are given) through tracking objects in your segmented images and using a
% GUI to discover and correct errors in your segmentation that lead to
% errors in cell tracking. Once you've segmentation errors and can track
% your cells-of-interest throughout the course of your movie, you can move
% on to analysis.

%--------------------------------------------------------------------------
% STEP 1: READ IN RAW IMAGES
%--------------------------------------------------------------------------

% While we technically will only be making measurements / doing analysis on
% a segmentation mask, it is helpful to view the segmentation mask as an
% overlay on top of the raw images. Therefore, we'll start by loading the
% raw images into MATLAB and storing them in a 3D tensor.

%--------------------------------------------------------------------------
% STEP 2: PIXEL CLASSIFICATION & SEGMENTATION
%--------------------------------------------------------------------------

% There are two ways of performing pixel classification. Either option is
% completed outside of MATLAB and then loaded in prior to detection of
% cells.

% OPTION 1: PIXEL CLASSIFICATION USING ILASTIK

% The first option for pixel classification is using the pixel
% classification workflow in Ilastik (ilastik.org), which transforms the
% image from 8-bit space (or whatever bit depth you're in), where pixel
% value represents fluorescence intensity to a new 8-bit space where pixel
% value represents the probability of being either a cell edge or not
% (where 0s represent 100% probability that these pixels are background and
% 255s represent 100% probability that these pixels are cell edges).
% https://www.ilastik.org/documentation/pixelclassification/pixelclassification


% OPTION 2: U-NET (or any other method of pixel classification that saves
% the result as a binary image)

% The second option for pixel classification is using U-NET or any other
% pixel classification workflow that creates binary images where 0s
% are pixels classified as background or cell interiors and 1s are pixels
% classified as being a cell edge.

%--------------------------------------------------------------------------
% STEP 3: DETECT CELLS - watershed transform & bwlabel
%--------------------------------------------------------------------------

% After transforming our images into a space where 0s represent background
% pixels or cell interior and 1s represent cell edges, we next need to
% detect the location of cells. To do this, we are going to use a watershed
% transform ( https://en.wikipedia.org/wiki/Watershed_(image_processing) , 
% https://www.mathworks.com/help/images/ref/watershed.html ) to define
% cells and clean up noise from the pixel classification, followed by a
% function called bwlabel that will assign identities to binary objects
% defined using a defined 2D connectivity.

%--------------------------------------------------------------------------
% STEP 4: TRACKING - hungarian (munkres) algorithm
%--------------------------------------------------------------------------

% Bwlabel gives cells a unique identify for every time point they exist. To
% track cells across time, we must create a map that connects cells between
% adject time points. We will be using the munkres assignment algorithm
% (https://en.wikipedia.org/wiki/Hungarian_algorithm , 
% https://www.mathworks.com/matlabcentral/fileexchange/20328-munkres-assignment-algorithm ). 

%--------------------------------------------------------------------------
% STEP 5: MANUAL CORRECTIONS - using the GUI
%--------------------------------------------------------------------------

% Try as we might, there is currently no methodology that can generate
% perfect segmentation. U-Net performed the best out of all methods we
% tested. However, it still had ~0.5% percent error in segmentation that,
% when tracked over 120 time points, compounded to over 10% error in
% tracking! Therefore, we developed a matlab GUI ('segmeter') that uses
% tracking errors to discover and correct the underlying segmentation
% errors


%%


%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

% step 1: read in raw images

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------


% it is useful to be able to view the pixel classification and watershed
% transform transform overlayed on top of the raw images, so we can make
% judgement calls about the quality of the watershed transform.  Therefore,
% we need to read in the raw images and store them in a matrix.


% set directory for data
image_dir = dir('/Users/kevin/Documents/ImagingData/WT_replicate2/raw/*tif');
folder = image_dir.folder;

% establish dimensionality of movie
test_image = imread(fullfile(folder,image_dir(1).name));
x_resolution = size(test_image,1);
y_resolution = size(test_image,2);
t_resolution = length(image_dir);

% initalize matrix to store raw images
raw_images = uint8(zeros(x_resolution,y_resolution,t_resolution));

       
% loop through time
for t = 1:t_resolution
    
    % read image and store
    raw_images(:,:,t) = imread(fullfile(folder,image_dir(t).name));
    
end

%% visualize

implay(raw_images)


%%


%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

% STEP 2 / OPTION 1: PIXEL classification using ilastik

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

% NOTE: it is important that 'Label 1' in Ilastik is the cell edges and
% 'Label 2' is the cell interior / background.

% Ilastik outputs images with probability values for each pixel in the
% h5 file format.  The following code is an h5 file reader that will read
% in this format and convert it into a matrix of doubles. It will load in
% all the .h5 files ending in '_Probabilities.h5' and will load them in the
% order of whatever iterator is in the file name (i.e. _T0075, _T0076, etc)

% read in h5 probability files from ilastik
ilastik_probabilities = read_ilastik_h5('/Users/kevin/Documents/ImagingData/WT_replicate2/raw/');


%% visualize

implay(ilastik_probabilities)


%% 

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

% STEP 3 / OPTION 1: watershed Ilastik pixel classification

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

% The 'ilastik_watershed' function uses a number of built in matlab
% functions, including bwlabel, to apply a watershed transform. Because the
% default output from bwlabel is called 'L', we'll be sticking with this
% convention and henceforth call our segmentation mask L.

% this can sometimes take awhile

L = ilastik_watershed(ilastik_probabilities);


%% visualize watershed as red skeleton over raw images

ilastik_seg_display = uint8(zeros(size(L,1),size(L,2),3,size(L,3)));

% iterate through time
for t = 1:size(L,3)
    t
    % make rgb version of current time frame
    ilastik_seg_display(:,:,:,t) = cat(3,raw_images(:,:,t),raw_images(:,:,t),raw_images(:,:,t));
    
    % find seg skeleton in L and overlay as redline (1,0,0)
    [x,y] = find(not(L(:,:,t)));
    for j = 1:length(x)
        ilastik_seg_display(x(j),y(j),1,t) = 255;
        ilastik_seg_display(x(j),y(j),2,t) = 0;
        ilastik_seg_display(x(j),y(j),3,t) = 0;
    end
    
end

implay(ilastik_seg_display)
    


%%

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

% STEP 2 / OPTION 2: U-NET (or alternative) pixel classification

% NOTE: pixel classification must be binary (can be 8-bit) with 0s
% representing background / cell interiors and 1s (or 255s if 8-bit)
% representing cell edges

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

% initalize matrix to store raw images - needs to be unit 8

% set directory
current_dir = dir('/Users/kevin/Documents/ImagingData/WT_replicate2/seg/*tif');
folder = current_dir.folder;

% establish dimensionality of movie
test_image = imread(fullfile(folder,current_dir(1).name));
x_resolution = size(test_image,1);
y_resolution = size(test_image,2);
t_resolution = length(current_dir);

seg_mask = zeros(x_resolution,y_resolution,t_resolution);
       
% loop through time
for t = 1:t_resolution
    
    % read image
    seg_mask(:,:,t) = uint8(imread(fullfile(current_dir(1).folder,current_dir(t).name)));
    
end

%% visualize

implay(seg_mask)


%% 

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

% STEP 3 / OPTION 2: watershed U-Net (or alternative) pixel classification

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

% Apply watershed transform and label objects using bwlabel. Because the
% default output from bwlabel is called 'L', we'll be sticking with this
% and call our segmentation mask L.

L = zeros(size(seg_mask));
for ii = 1:size(seg_mask,3)
    L(:,:,ii) = bwlabel(watershed(seg_mask(:,:,ii),8));
end

%% visualize watershed as red skeleton over raw images

unet_seg_display = uint8(zeros(size(L,1),size(L,2),3,size(L,3)));

% iterate through time
for t = 1:size(L,3)
    t
    % make rgb version of current time frame
    unet_seg_display(:,:,:,t) = cat(3,raw_images(:,:,t),raw_images(:,:,t),raw_images(:,:,t));
    
    % find seg skeleton in L and overlay as redline (1,0,0)
    [x,y] = find(not(L(:,:,t)));
    for j = 1:length(x)
        unet_seg_display(x(j),y(j),1,t) = 255;
        unet_seg_display(x(j),y(j),2,t) = 0;
        unet_seg_display(x(j),y(j),3,t) = 0;
    end
    
end

implay(unet_seg_display)

%%

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

% tracking

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

% Once we have found cells in the segmented images (we stored this as the 
% 'L' matrix), we can track cells across time. Bwlabel gives every cell a 
% unique ID. However, these IDs only refer to the time point where the
% object exists. Therefore, we must map these together across time. This is
% what we mean by 'tracking'.

%% tracking

% Parameter for maximum distance cells are permitted to move between time
% points. This does not just exist for the sake of increasing the
% efficiency of the munkres assignment algorithm. It also constrains the
% number of possible solutions - i.e. if the cutoff is too large, it will
% create strings of mis-tracked cells connected instances where a cell
% appears on one side of the FOV and disappears on the other side of the
% FOV in the same time point. If the cutoff is too small, cells tracks will
% be lost in regions where there is rapid cell expansion and contraction
% (near division events), which result in large cell movements over a
% single time frame.
cutoff = 30;    % WT1

% The manual correction GUI changes the segmentation mask and the matrix
% where we labeled cells ('L'). When you export info from the GUI back into
% your matlab workspace, its recorded as a separate variable (ending in
% '_GUI'). This line will check for whether 'L_GUI' exists and update 'L'.
% If you don't want to incorporate the changes you made using the GUI, then
% you should delete 'L_GUI' from your workspace.
if exist('L_GUI','var')
    L = L_GUI;
end

% record info about cell centroids and area in cell arrays
[cell_area, cell_centroid] = cellStats(L);

% This is the wrapper function for the munkres assignment algorithm. It
% outputs a matrix ('cellPairs') that contains timepoint-to-timepoint maps
% for cells across the movie. cellPairs gives you the index of cell 'c' at
% timepoint t-1. i.e. cellPairs{20}(100) gives you the corresponding index
% of cell 100 at timepoint 20 at timepoint 19.
[ cellPairs ] = hungarian( L, cutoff );

% The GUI for manually correcting mistakes in segmentation also allows for 
% manual correction in tracking. If any manual corrections were made, these
% recorded in the matrix 'track_links'. If this matrix exists, we need to
% use it to modify the timepoint-to-timepoint maps in 'cellPairs'
if exist('track_links','var')
    cellPairs = manualUpdateCellPairs(cellPairs, track_links);
end

% Instead of calling cellPairs{t}(c) to find the map of cell 'c' between
% time point 't' and 't-1', it is easier to assemble a matrix (called 
% 'tracks') that spans all of time ('columns') and contains the unique
% label from 'L' in every row.
[tracks, display_tracks] = generateTracksGUI(L, cell_area, cellPairs);

% Find cell births and deaths. Cell births are defined as cell tracks that
% begin later than the first time point and that do not originate along the
% boundary of the segmented FOV. Cell deaths are defined as cell tracks
% that end prior to the last time point and do not end along the boundary
% of the FOV. 'birth_death_display' will be used in the manual correction
% GUI in order to evaluate each birth/death to determine whether it is a
% real biological phenomena or the result of a segmentation error. This is
% how we discover mistakes in segmentation
[birth,death,total_born,total_dead,boundary_cells,birth_death_display] = ...
    birthDeath(L,tracks,cell_area);



%%

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

% visualize tracking mask as transparency over raw images

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

% display tracks as color mask 
for t = 1:t_resolution
    
    % concatenate 8-bit raw into rgb
    rgbImage = cat(3, raw_images(:,:,t), raw_images(:,:,t), raw_images(:,:,t));
    
    % set transparency for overlay
    transparency = 1 - im2double(raw_images(:,:,t)) - 0.4;
    %transparency = im2double(raw_images(:,:,t));
    
    % show raw
    imshow(rgbImage)
    hold on
    
    % show tracking color masks
    h = imshow(display_tracks(:,:,:,t));
    hold off
    
    % apply transparency
    set(h, 'AlphaData',transparency);
    
    % draw
    drawnow
    
end

%% display tracks as color mask 

implay(display_tracks)


%%

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

% manual correction of segmentation errors

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

% first, we need to different version of the data the GUI uses for
% visualization the raw images and segmentation masks

% real raw images
raw_images1 = raw_images;

% if using ilastik
raw_images2 = ilastik_probabilities;

% if using U-Net
% raw_images2 = seg_mask;

% open GUI
segmenter

% after you've finished correcting errors, you should export to your matlab
% workspace and then re-run the tracking code above.
