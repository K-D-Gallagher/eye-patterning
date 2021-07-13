# eye-morphogenesis [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Table of Contents
- [PART 1: SEGMENTATION AND TRACKING (SEG_TRACKING_DRIVER.M)](#part-1-segmentation-and-tracking-seg_tracking_driverm)
- [STEP 1: READ IN RAW IMAGES](#step-1-read-in-raw-images)
- [STEP 2: PIXEL CLASSIFICATION & SEGMENTATION](#step-2-pixel-classification--segmentation)
- [STEP 3: DETECT CELLS - watershed transform & bwlabel](#step-3-detect-cells---watershed-transform--bwlabel)
- [STEP 4: TRACKING - hungarian (munkres) algorithm](#step-4-tracking-hungarian-munkres-algorithm)
- [STEP 5: MANUAL CORRECTIONS - using the GUI](#step-5-manual-corrections-using-the-gui)
- [PART 2: ANALYSIS (ANALYSIS_DRIVER.M)](#part-2-analysis-analysis_driverm)
- [Installation](#installation)
- [License](#license)

--------------------------------------------------------------------------
--------------------------------------------------------------------------
--------------------------------------------------------------------------
## PART 1: SEGMENTATION AND TRACKING (SEG_TRACKING_DRIVER.M)
--------------------------------------------------------------------------
--------------------------------------------------------------------------
--------------------------------------------------------------------------

The goal of this driver file is to segment and track images. Unless you
find a way to achieve perfect pixel classification, this will
unfortunately involve some manual correction. This driver file brings you
from initial pixel classification using an external program (some options
are given) through finding and tracking objects in your segmented images
and using a GUI to discover and correct errors in your segmentation that
lead to errors in cell tracking.

--------------------------------------------------------------------------
STEP 1: READ IN RAW IMAGES
--------------------------------------------------------------------------

While we technically will only be making measurements / doing analysis on
a segmentation mask, it is helpful to view the segmentation mask as an
overlay on top of the raw images. Therefore, we'll start by loading the
raw images into MATLAB and storing them in a 3D tensor.

--------------------------------------------------------------------------
STEP 2: PIXEL CLASSIFICATION & SEGMENTATION
--------------------------------------------------------------------------

There are two ways of performing pixel classification. Either option is
completed outside of MATLAB and then loaded in prior to detection of
cells.

OPTION 1: PIXEL CLASSIFICATION USING ILASTIK

The first option for pixel classification is using the pixel
classification workflow in Ilastik (ilastik.org), which transforms the
image from 8-bit space (or whatever bit depth you're in), where pixel
value represents fluorescence intensity to a new 8-bit space where pixel
value represents the probability of being either a cell edge or not
(where 0s represent 100% probability that these pixels are background and
255s represent 100% probability that these pixels are cell edges).
https://www.ilastik.org/documentation/pixelclassification/pixelclassification


OPTION 2: U-NET (or any other method of pixel classification that saves
the result as a binary image)

The second option for pixel classification is using U-NET or any other
pixel classification workflow that creates binary images where 0s
are pixels classified as background or cell interiors and 1s are pixels
classified as being a cell edge.

--------------------------------------------------------------------------
STEP 3: DETECT CELLS - watershed transform & bwlabel
--------------------------------------------------------------------------

After transforming our images into a space where 0s represent background
pixels or cell interior and 1s represent cell edges, we next need to
detect the location of cells. To do this, we are going to use a watershed
transform ( https://en.wikipedia.org/wiki/Watershed_(image_processing) , 
https://www.mathworks.com/help/images/ref/watershed.html ) to define
cells and clean up noise from the pixel classification, followed by a
function called bwlabel that will assign identities to binary objects
defined using a defined 2D connectivity.

--------------------------------------------------------------------------
STEP 4: TRACKING - hungarian (munkres) algorithm
--------------------------------------------------------------------------

Bwlabel gives cells a unique identify for every time point they exist. To
track cells across time, we must create a map that connects cells between
adject time points. We will be using the munkres assignment algorithm
(https://en.wikipedia.org/wiki/Hungarian_algorithm , 
https://www.mathworks.com/matlabcentral/fileexchange/20328-munkres-assignment-algorithm ). 

--------------------------------------------------------------------------
STEP 5: MANUAL CORRECTIONS - using the GUI
--------------------------------------------------------------------------

Try as we might, there is currently no methodology that can generate
perfect segmentation. U-Net performed the best out of all methods we
tested. However, it still had ~0.5% percent error in segmentation that,
when tracked over 120 time points, compounded to over 10% error in
tracking! Therefore, we developed a matlab GUI ('segmeter') that uses
tracking errors to discover and correct the underlying segmentation
errors. Tutorial video pending.

--------------------------------------------------------------------------
--------------------------------------------------------------------------
--------------------------------------------------------------------------
## PART 2: ANALYSIS (ANALYSIS_DRIVER.M)
--------------------------------------------------------------------------
--------------------------------------------------------------------------
--------------------------------------------------------------------------

This driver file (ANALYSIS_DRIVER.m) contains code demonstrating how to 
access our data and the various piece of data annotation we relied on
for our analysis. We hope this provides the foundation for future analysis
and exploration of the data.

--------------------------------------------------------------------------
Cell tracking
--------------------------------------------------------------------------

Cells are

![tracking_example](github_media/tracking.gif)

--------------------------------------------------------------------------
Annotation of ommatidial preclusters and photoreceptor classes
--------------------------------------------------------------------------

![precluster_ID_example](github_media/preclusters.gif)

--------------------------------------------------------------------------
Morphogenetic furrow location
--------------------------------------------------------------------------

![morphogenetic_furrow_location_example](github_media/MF.gif)

--------------------------------------------------------------------------
Cell velocity
--------------------------------------------------------------------------

![cell_velocity_example](github_media/velocity.gif)

--------------------------------------------------------------------------
Ommatidial lattice column annotation
--------------------------------------------------------------------------

![ommatidial_lattice_column_annotation_example](github_media/columns.gif)

--------------------------------------------------------------------------
Ommatidial lattice row annotation
--------------------------------------------------------------------------

![ommatidial_lattice_row_annotation_example](github_media/rows.gif)

--------------------------------------------------------------------------
Ommatidial lattice annotation
--------------------------------------------------------------------------

![ommatidial_lattice_annotation_example](github_media/lattice.gif)


## Installation
We recommend using MATLAB R2018a or newer. We have not tested our code with older versions of MATLAB. This code additionally requires MATLAB's Image Processing Toolbox.

## License
The MIT License (MIT)

Copyright (c) 2021 Kevin D Gallagher

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

## Questions
Reach out to me at kevin.d.gallagher2@gmail.com if you have any questions about this repository or code.
Contact Madhav Mani (madhav.mani@northwestern.edu) or Richard Carthew (r-carthew@northwestern.edu) with questions pertaining to the affiliated publication
