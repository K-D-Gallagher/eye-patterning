function [nonborder_cells] = findNonBorderCells(L)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function finds all non-border cells after segmentation (generation
% of L).  While L does not display border cells if 'clean up borders' is
% selected when generating Struct, the border cells still exist in Struct
% and will be included in tracking.  To avoid this, I created this function
% that reads the unique indices of L with the border cells hidden, thereby
% creating a matrix containing non-border cells for each time point.  This
% matrix can then be used for tracking non-border cells
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% find number of time points
t_resolution = size(L,3);

% initialize cell array to store cell indices from each time point
nonborder_cells = cell(t_resolution,1);

% iterate through time
for ii = 1:t_resolution
    
    % find unique pixel values for current time-frame of L, which
    % represents the indices of all objects
    cells = unique(L(:,:,ii));
    
    % input everything in this list into the cell array storing cell
    % indices except for the first entry, which is always 0 and represents
    % the boundaries, and the second entry, which is always 1 and
    % represents the boundary
    nonborder_cells{ii} = cells(2:end);
    
end

end
