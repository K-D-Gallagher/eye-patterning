function [tracks, display_tracks] = generateTracksGUI(L, cell_area, cellPairs)

% Note, this is the same as 'generateTracksMinimal', except that it
% requires that cellPairs be passed in as an argument

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2: find non-border cells
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nonborder_cells = findNonBorderCells(L);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3: Generate cell array to store tracks. Before creating tracks, we need to
% populate it with the indices of cells in the final time point
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get dimensions of the image stack
x_resolution = size(L,1);
y_resolution = size(L,2);
t_resolution = size(L,3);

% initialize cell array to store tracks
tracking_array = cell(zeros,t_resolution);

% populate the last cell of tracks with the indices of cells from the last
% time point
for ii = 1:length(nonborder_cells{t_resolution})
    tracking_array{t_resolution}(ii) = nonborder_cells{t_resolution}(ii);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4: Create tracks from cellPairs.  Tracks will be represented as indices
% within the tracking_array.  If cells are successfully mapped between
% adjacent time points (by cellPairs), then their unique indices (from
% bwlabel) will be stored in the same track-index.  However, if they are
% not mapped by cellPairs, then a new track will be generated.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% iterate through time
for t = 1:t_resolution-1
    
    % store current time as an easy-to-read variable
    current_time = t_resolution-t;
    
    % create a counter to keep track of how many cells we'll be indexing -
    % this is specifically for adding cells to the tracking_array when
    % there is no match in cellPairs - we'll start by taking the number of
    % unique tracks from the previous time point (t+1)
    current_length = length(tracking_array{current_time+1});
    
    previous_tracks = tracking_array{current_time+1};
    
    % iterate through cells of current time point
    for ii = 1:length(nonborder_cells{current_time})
        
        % store index of current cell
        current_index = nonborder_cells{current_time}(ii);
        
        % store indices mapped from t+1 to current time point by cellPairs
        mapped_indices = cellPairs{current_time+1};
        
        % look for match between current index and indices mapped from t+1
        % by cellPairs
        test_match = ismember(current_index, mapped_indices);
        
        % if current index was successfully mapped by cellPairs
        if test_match ~= 0
            
            % find the index of the current cell's mapped partner from
            % the t+1 timepoint in cellPairs. this is the index of the
            % current cell in the t+1 time point.
            partner_index = find(mapped_indices == current_index);
            
            % find the location of the partner index (representing the t+1
            % index corresponding to the current_index) in tracks from the
            % previous time point, which we vectorized in the
            % inner-time-loop
            tracking_partner_index = find(previous_tracks == partner_index);
            
            % insert the current_index into the appropriate track.  note
            % that tracks are represented by indices within tracking_array
            tracking_array{current_time}(tracking_partner_index) = current_index;

        end
        
        % if current index IS NOT mapped in cellPairs
        if test_match == 0
            
            % update the counter to signify creation of a new track
            current_length = current_length+1;
            
            % create new track in tracking_array (counter 'current_length'
            % should be a greater than the largest entry in tracking_array
            % for the current time point)
            tracking_array{current_time}(current_length) = current_index;

        end
    end
end






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 5: Pull tracks out of cell array and store in matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%
% First, let's fine the largest number of cells we need to account for 
%%%%%

% initialize counter that records the largest number of tracks from
% tracking_array
max_cells = 0;

% iterate through time
for ii = 1:t_resolution
    
    % store the number of tracks from the current time point of
    % tracking_array
    current_length = length(tracking_array{ii});
    
    % if the current number of tracks is greater than the previous largest
    % number of tracks..
    if current_length > max_cells
        
        % then update max_cells to the current number of tracks
        max_cells = current_length;

    end
end

% initialize matrix to store tracks. it will be the greatest number of
% tracks x time
tracks = zeros(max_cells, t_resolution);    % (tracks x time)

%%%%%
% now let's populate the matrix with tracks from the remaining time points
%%%%%

% loop through time
for ii = 1:t_resolution
    
    % store list of current cell indices
    current_cell_list = tracking_array{ii};
    
    % loop through list of current indices
    for iii = 1:length(current_cell_list)
        
        % insert current index from tracking_array into the current
        % position within tracks
        tracks(iii,ii) = current_cell_list(iii);

    end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 6: Create a matrix to visualize the tracks.  This will be achieved by
% uniquely coloring each track with random numbers.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize matrix to display tracks
display_tracks = zeros(x_resolution, y_resolution, 3, t_resolution);

% create an RGB random color vector to uniquely color each track
random_colors = rand(length(tracks(:,1)),3);

% iterate through time
for t = 1:t_resolution
    
    % iterate through tracks
    for ii = 1:length(tracks(:,t))
        
        % store current index (unique index in L for cell at current time
        % in the current rack
        current_cell_index = tracks(ii,t);
        
        % check to make sure a cell actually exists at this time point
        if current_cell_index ~= 0 
            
            % if it does, then find the pixels indices representing its
            % area
            x_pixels = cell_area{t}{current_cell_index}(:,1);
            y_pixels = cell_area{t}{current_cell_index}(:,2);
            
            % iterate through this pixel-area indices
            for iii = 1:length(x_pixels)
                
                % color in the cell area based on the random-color vector
                display_tracks(x_pixels(iii),y_pixels(iii),1,t) = random_colors(ii,1);
                display_tracks(x_pixels(iii),y_pixels(iii),2,t) = random_colors(ii,2);
                display_tracks(x_pixels(iii),y_pixels(iii),3,t) = random_colors(ii,3);
                
            end
            
        end
        
    end
end



end
