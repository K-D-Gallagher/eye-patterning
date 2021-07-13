# eye-morphogenesis [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

&nbsp;

## Table of Contents
- [PART 1: ANALYSIS (ANALYSIS_DRIVER.M)](#part-2-analysis-analysis_driverm)
    - [Cell tracking](#cell-tracking)
    - [Annotation of cell classes](#annotation-of-cell-classes)
    - [Morphogenetic furrow location](#morphogenetic-furrow-location)
    - [Cell centroid velocity field](#cell-centroid-velocity-field)
    - [Ommatidial lattice annotation](#ommatidial-lattice-annotation)
    - [Ommatidial lattice annotation - column annotation](#ommatidial-lattice-annotation---column-annotation)
    - [Ommatidial lattice annotation - row annotation](#ommatidial-lattice-annotation---row-annotation)
- [PART 2: SEGMENTATION AND TRACKING (SEG_TRACKING_DRIVER.M)](#part-1-segmentation-and-tracking-seg_tracking_driverm)
    - [STEP 1: READ IN RAW IMAGES](#step-1-read-in-raw-images)
    - [STEP 2: PIXEL CLASSIFICATION & SEGMENTATION](#step-2-pixel-classification--segmentation)
    - [STEP 3: DETECT CELLS - watershed transform & bwlabel](#step-3-detect-cells---watershed-transform--bwlabel)
    - [STEP 4: TRACKING - hungarian (munkres) algorithm](#step-4-tracking---hungarian-munkres-algorithm)
    - [STEP 5: MANUAL CORRECTIONS - using the GUI](#step-5-manual-corrections---using-the-gui)
- [Installation](#installation)
- [License](#license)

&nbsp;
&nbsp;
&nbsp;

# PART 1: ANALYSIS (ANALYSIS_DRIVER.M)

The data we are publishing is very rich and contains many more phenomena than the ones we've reported on. This driver file (ANALYSIS_DRIVER.m) contains code demonstrating how to access our data and the various ways we've annotation the data to extract measurements from subpopulations of cells and/or different regions of interest in the tissue. We hope this code helps facilitate further analysis from researchers interested in expanding on our work or, even better, asking new questions entirely!


## Cell tracking

After [segmenting cells](#part-1-segmentation-and-tracking-seg_tracking_driverm), we use the [munkres assignment algorithm](#step-4-tracking---hungarian-munkres-algorithm) to map cells between adjacent time points. This is prerequisite before any additional analysis can be completed. Here, we demonstrate how to display a representation of tracked cells via a color mask that uses a unique, random color to label each tracked object.

![tracking_example](github_media/tracking.gif)

&nbsp;

## Annotation of cell classes

With cells segmented and tracked, we were able to annotate the identities of photoreceptor cells (R-cells) onto our data manually. Cell classes were easily determined by their unique morphology once they've differentiated on the posterior side of the morphogenetic furrow. These identities were then propegated backwards in time, allowing us to extract measurements of position, topology, and morphology from these cells before they were committed to these fates. This code shows how to access these cell class annotations and color in the cell area corresponding to each instance of a R-cell throughout the duration of the movie.

![precluster_ID_example](github_media/preclusters.gif)

&nbsp;

## Morphogenetic furrow location

The morphogenetic furrow is the wavefront of differenetiation in this system. As it moves from the posterior to anterior margin of the eye imaginal disc, it triggers a wave of simultaneous cell differentiation and morphological changes. Two scales of patterning happen concurrently in the wake of the morphogenetic furrow: 1) locally, cells lined up in the dorsal-ventral direction on the posterior side of the morphogenetic furrow buckle into hairpin structures and close into multicellular rosettes; these are the first photoreceptor cells to differentiate. 2) each multicellular rosette emerges at the precise location to assimilate into the already established triangular lattice of rosettes on the posterior side of the tissue, thus growing the lattice in a process reminiscent of directional solidification in the crystal manufacturing process. Here, we demonstrate code that allows one to access and plot the location of the morphogenetic furrow at each time point. This data is available for all four datasets.

![morphogenetic_furrow_location_example](github_media/MF.gif)

&nbsp;

## Cell centroid velocity field

The cell centroid velocity field is calculated through central displacement of cell centroids over one hour of developmental time (the movie has a frame rate of 5 minutes, so this is 12 frames of the movie). One hour was found to be the timescale over which deterministic behavior emerged. The velocity field derivated from shorter time steps introduced additional noise to the organization of the velocity field, whereas timesteps greater than an hour did not substantial change the behavior/organization of the velocity field. Here, we demonstrate how to visualize the velocity field using matlab's quiver function.

![cell_velocity_example](github_media/velocity.gif)

&nbsp;

## Ommatidial lattice annotation

The ommatidial lattice emerges concurrently with the cell differentiation / local morphological changes triggered by passage of the morphogenetic furrow. A triangular lattice is dual to a hexoganl lattice, so what we are watching emerge here is the initial nucleation of the final ommatidial lattice, visible on the surface of the adult eye. Here, we demonstrate how to visualize the ommatidial lattice, defined using the centroids of the R8 class photoreceptor cells, using matlab's patch function.

![ommatidial_lattice_annotation_example](github_media/lattice.gif)

&nbsp;

## Ommatidial lattice annotation - column annotation

We define columns of ommatidia as being parallel to the morphogenetic furrow. In the publication, the morphogenetic furrow is oriented vertically with the anterior to the left and posterior to the right. Here, things are rotated 90 degrees and the morphogenetic furrow is horizontal with the anterior towards the bottom and posterior towards the top. Note how cells fated to belong to separate columns ommatidia posterior of the morphogenetic furrow are compressed into overlapping positions along the anterior-posterior axis inside the morphogenetic furrow. Here, we demonstrate how to color in R-cells according to their column identity within the ommatidial lattice.

![ommatidial_lattice_column_annotation_example](github_media/columns.gif)

&nbsp;

## Ommatidial lattice annotation - row annotation

Similar to column identity, we can also define rows of ommatidia as being perpendicular to the morphogenetic furrow. Note how the movement of ommatidia in adjacent rows inside the morphogenetic is antiphase: i.e. one row will flow with the morphogenetic furrow, while the other will become stationary and assimilate into the ommatidial lattice posterior of the morphogeneticu furrow. Here, we demonstrate how to color in R-cells according to their row identity within the ommatidial lattice.

![ommatidial_lattice_row_annotation_example](github_media/rows.gif)

&nbsp;
&nbsp;
&nbsp;

# PART 2: SEGMENTATION AND TRACKING (SEG_TRACKING_DRIVER.M)

The goal of this driver file is to segment and track images. The intention of this driver file is to document the pipeline we used to process data for our publication, but also to lay the framework for others to process their own data. Unless you find a way to achieve perfect pixel classification, this will unfortunately involve some manual correction. This driver file brings you from initial pixel classification using an external program (Ilastik or U-Net) through finding and tracking objects in your segmented images in matlab and using a GUI to discover and correct errors in your segmentation that lead to errors in cell tracking.

&nbsp;

## STEP 1: READ IN RAW IMAGES

While we technically will only be making measurements / doing analysis on a segmentation mask, it is helpful to view the segmentation mask as an overlay on top of the raw images. Therefore, we'll start by loading the raw images into MATLAB and storing them in a 3D tensor.

![raw_images_example](github_media/raw.gif)

&nbsp;

## STEP 2: PIXEL CLASSIFICATION & SEGMENTATION

There are two ways of performing pixel classification. Either option is completed outside of MATLAB and then loaded in prior to detection of cells.

&nbsp;

### OPTION 1: PIXEL CLASSIFICATION USING ILASTIK

The first option for pixel classification is using the pixel classification workflow in Ilastik (ilastik.org), which transforms the image from 8-bit space (or whatever bit depth you're in), where pixel value represents fluorescence intensity to a new 8-bit space where pixel value represents the probability of being either a cell edge or not (where 0s represent 100% probability that these pixels are background and 255s represent 100% probability that these pixels are cell edges).
https://www.ilastik.org/documentation/pixelclassification/pixelclassification

&nbsp;

Example of pixel classification process in Ilastik:
![ilastik_demo](github_media/ilastik_demo.gif)

&nbsp;

Example of pixel classification from Ilastik:
![ilastik_seg_example](github_media/ilastik_seg.gif)

&nbsp;

### OPTION 2: U-NET (or any other method of pixel classification that saves the result as a binary image)

The second option for pixel classification is using U-NET or any other pixel classification workflow that creates binary images where 0s are pixels classified as background or cell interiors and 1s are pixels classified as being a cell edge.

![u-net](github_media/U-Net.gif)

&nbsp;

## STEP 3: DETECT CELLS - watershed transform & bwlabel

After transforming our images into a space where 0s represent background pixels or cell interior and 1s represent cell edges, we next need to detect the location of cells. To do this, we are going to use a watershed transform ( https://en.wikipedia.org/wiki/Watershed_(image_processing) , https://www.mathworks.com/help/images/ref/watershed.html ) to define objects within the pixel classified images and clean up noise from the pixel classification, followed by a function called bwlabel that will assign identities to binary objects defined using a chosen 2D connectivity.

&nbsp;

## STEP 4: TRACKING - hungarian (munkres) algorithm

Bwlabel gives cells a unique identify for every time point they exist. To track cells across time, we must create a map that connects cells between adject time points. We will be using the munkres assignment algorithm (https://en.wikipedia.org/wiki/Hungarian_algorithm , https://www.mathworks.com/matlabcentral/fileexchange/20328-munkres-assignment-algorithm ). 

&nbsp;

## STEP 5: MANUAL CORRECTIONS - using the GUI

Try as we might, there is currently no methodology that can generate perfect segmentation. U-Net performed the best out of all methods we tested. However, it still had ~0.5% percent error in segmentation that, when tracked over 120 time points, compounded to over 10% error in tracking! Therefore, we developed a matlab GUI ('segmeter') that uses tracking errors to discover and correct the underlying segmentation errors. Tutorial video pending.

&nbsp;
&nbsp;
&nbsp;

# Installation
We recommend using MATLAB R2018a or newer. We have not tested our code with older versions of MATLAB. This code additionally requires MATLAB's Image Processing Toolbox.

&nbsp;

# License
The MIT License (MIT)

Copyright (c) 2021 Kevin D Gallagher

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

&nbsp;

# Questions
Reach out to me at kevin.d.gallagher2@gmail.com if you have any questions about this repository or code.
Contact Madhav Mani (madhav.mani@northwestern.edu) or Richard Carthew (r-carthew@northwestern.edu) with questions pertaining to the affiliated publication
