function [Velocity,calculation_time_range] = velocityField(L, cell_area, cell_centroid, tracks ,...
    calculation_type, displacement_time, interp_subsample_factor, ...
    velocity_scaling_factor,interp_time_range)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% centroid based velocity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% calculate velocities and store in cell arrays - Note that it's important
% that we move the centroid coordinates WITH the velocity readings, as we
% will be plotting the velocities from the centroid positions.  Like with
% the rest of my code, the tracking index is equivalent to the tracking ID.
% ultimately, by creating new structures that contain velocities and their
% corresponding centroid positions, we are removing this information from
% the tracking infrastructure.  new code will need to be developed that
% stores velocity information alongside tracking information.  the purpose
% of this code is simply to allow us to visual the
% centroid-based-velocity-field

time = size(L,3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find range for interpolating velocity field
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
temp_x_min = zeros(interp_time_range,1);
temp_x_max = zeros(interp_time_range,1);
temp_y_min = zeros(interp_time_range,1);
temp_y_max = zeros(interp_time_range,1);

% find bounding box for each time point
for t = 1:interp_time_range
    fill_mask = regionprops(imcomplement(L(:,:,t)),'BoundingBox');
    temp_y_min(t) = round(fill_mask.BoundingBox(2))+30;
    temp_y_max(t) = temp_y_min(t) + fill_mask.BoundingBox(4)-50;  % adjusted to make range of 380
    temp_x_min(t) = round(fill_mask.BoundingBox(1))+30;
    temp_x_max(t) = temp_x_min(t) + round(fill_mask.BoundingBox(3))-50;
end

% find absolute min and max
x_min = min(temp_x_min) + 25;
x_max = max(temp_x_max) - 25;
y_min = min(temp_y_min) + 25;
y_max = min(temp_y_max) - 25;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate displacements of centroids
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% [x-component velocity, time] & [y-component velocity, time]
CVX = zeros(size(tracks,1),time);
CVY = zeros(size(tracks,1),time);

% [x-component centroid, time] & [y-component centroid, time]
CX = zeros(size(tracks,1),time);
CY = zeros(size(tracks,1),time);

% Note: because the dimensions of these matrices match those of the
% tracking matrix, tracking information is preserved

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALCULATION TYPE 1
% calculate velocity by central averaging across pairs of adjacent time
% points, then averaging all of these averages
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if calculation_type == 1
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % CALCULATION TYPE 1
    % DISPLACEMENT TIME 3
    % calculate velocity using +/- 1 time point
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if displacement_time == 3
        
        for t = 2:time-2
            
            future_centroids1 = cell_centroid{t+1};
            current_centroids = cell_centroid{t};
            previous_centroids = cell_centroid{t-1};
            
            % search through tracks at this time point
            for ii = 1:size(tracks,1)
                
                % store current cell
                future_cell1 = tracks(ii,t+1);
                current_cell = tracks(ii,t);
                previous_cell = tracks(ii,t-1);
                
                % the cell exists for both the current and previous time points,
                % then we'll compute it's velocity
                if future_cell1 > 0 && current_cell > 0 && previous_cell > 0
                    
                    deltaX1 = future_centroids1{future_cell1}(1) - current_centroids{current_cell}(1);
                    deltaY1 = future_centroids1{future_cell1}(2) - current_centroids{current_cell}(2);
                    
                    deltaX2 = current_centroids{current_cell}(1) - previous_centroids{previous_cell}(1);
                    deltaY2 = current_centroids{current_cell}(2) - previous_centroids{previous_cell}(2);
                    
                    % average the averages
                    deltaX = (deltaX1 + deltaX2)/2;
                    deltaY = (deltaY1 + deltaY2)/2;
                    
                    % store velocity components
                    CVX(ii,t) = deltaX;
                    CVY(ii,t) = deltaY;
                    
                    % create matrix storing centroids
                    CX(ii,t) = current_centroids{current_cell}(1);
                    CY(ii,t) = current_centroids{current_cell}(2);
                    
                end
            end
        end
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % CALCULATION TYPE 1
    % DISPLACEMENT TIME 5
    % calculate velocity using +/- 2 time points
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    elseif displacement_time == 5
        
        for t = 3:time-3
            
            future_centroids2 = cell_centroid{t+2};
            future_centroids1 = cell_centroid{t+1};
            current_centroids = cell_centroid{t};
            previous_centroids = cell_centroid{t-1};
            previous_centroids2 = cell_centroid{t-2};
            
            % search through tracks at this time point
            for ii = 1:size(tracks,1)
                
                % store current cell
                future_cell2 = tracks(ii,t+2);
                future_cell1 = tracks(ii,t+1);
                current_cell = tracks(ii,t);
                previous_cell = tracks(ii,t-1);
                previous_cell2 = tracks(ii,t-2);
                
                % the cell exists for both the current and previous time points,
                % then we'll compute it's velocity
                if future_cell2 > 0 && future_cell1 > 0 && current_cell > 0 && ...
                        previous_cell > 0 && previous_cell2 > 0
                    
                    deltaX1 = future_centroids2{future_cell2}(1) - future_centroids1{future_cell1}(1);
                    deltaY1 = future_centroids2{future_cell2}(2) - future_centroids1{future_cell1}(2);
                    
                    deltaX2 = future_centroids1{future_cell1}(1) - current_centroids{current_cell}(1);
                    deltaY2 = future_centroids1{future_cell1}(2) - current_centroids{current_cell}(2);
                    
                    deltaX3 = current_centroids{current_cell}(1) - previous_centroids{previous_cell}(1);
                    deltaY3 = current_centroids{current_cell}(2) - previous_centroids{previous_cell}(2);
                    
                    deltaX4 = previous_centroids{previous_cell}(1) - previous_centroids2{previous_cell2}(1);
                    deltaY4 = previous_centroids{previous_cell}(2) - previous_centroids2{previous_cell2}(2);
                    
                    % average the averages
                    deltaX = (deltaX1 + deltaX2 + deltaX3 + deltaX4)/4;
                    deltaY = (deltaY1 + deltaY2 + deltaY3 + deltaY4)/4;
                    
                    % store velocity components
                    CVX(ii,t) = deltaX;
                    CVY(ii,t) = deltaY;
                    
                    % create matrix storing centroids
                    CX(ii,t) = current_centroids{current_cell}(1);
                    CY(ii,t) = current_centroids{current_cell}(2);
                    
                end
            end
        end
        
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % CALCULATION TYPE 1
    % DISPLACEMENT TIME 7
    % calculate velocity using +/- 3 time points
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    elseif displacement_time == 7
        
        for t = 4:time-4
            
            future_centroids3 = cell_centroid{t+3};
            future_centroids2 = cell_centroid{t+2};
            future_centroids1 = cell_centroid{t+1};
            current_centroids = cell_centroid{t};
            previous_centroids = cell_centroid{t-1};
            previous_centroids2 = cell_centroid{t-2};
            previous_centroids3 = cell_centroid{t-3};
            
            % search through tracks at this time point
            for ii = 1:size(tracks,1)
                
                % store current cell
                future_cell3 = tracks(ii,t+3);
                future_cell2 = tracks(ii,t+2);
                future_cell1 = tracks(ii,t+1);
                current_cell = tracks(ii,t);
                previous_cell = tracks(ii,t-1);
                previous_cell2 = tracks(ii,t-2);
                previous_cell3 = tracks(ii,t-3);
                
                % the cell exists for both the current and previous time points,
                % then we'll compute it's velocity
                if future_cell3 > 0 && future_cell2 > 0 && future_cell1 > 0 && ...
                        current_cell > 0 && previous_cell > 0 && previous_cell2 > 0 && previous_cell3 > 0
                    
                    deltaX1 = future_centroids3{future_cell3}(1) - future_centroids2{future_cell2}(1);
                    deltaY1 = future_centroids3{future_cell3}(2) - future_centroids2{future_cell2}(2);
                    
                    deltaX2 = future_centroids2{future_cell2}(1) - future_centroids1{future_cell1}(1);
                    deltaY2 = future_centroids2{future_cell2}(2) - future_centroids1{future_cell1}(2);
                    
                    deltaX3 = future_centroids1{future_cell1}(1) - current_centroids{current_cell}(1);
                    deltaY3 = future_centroids1{future_cell1}(2) - current_centroids{current_cell}(2);
                    
                    deltaX4 = current_centroids{current_cell}(1) - previous_centroids{previous_cell}(1);
                    deltaY4 = current_centroids{current_cell}(2) - previous_centroids{previous_cell}(2);
                    
                    deltaX5 = previous_centroids{previous_cell}(1) - previous_centroids2{previous_cell2}(1);
                    deltaY5 = previous_centroids{previous_cell}(2) - previous_centroids2{previous_cell2}(2);
                    
                    deltaX6 = previous_centroids2{previous_cell2}(1) - previous_centroids3{previous_cell3}(1);
                    deltaY6 = previous_centroids2{previous_cell2}(2) - previous_centroids3{previous_cell3}(2);
                    
                    % average the averages
                    deltaX = (deltaX1 + deltaX2 + deltaX3 + deltaX4 + deltaX5 + deltaX6)/6;
                    deltaY = (deltaY1 + deltaY2 + deltaY3 + deltaY4 + deltaY5 + deltaY6)/6;
                    
                    % store velocity components
                    CVX(ii,t) = deltaX;
                    CVY(ii,t) = deltaY;
                    
                    % create matrix storing centroids
                    CX(ii,t) = current_centroids{current_cell}(1);
                    CY(ii,t) = current_centroids{current_cell}(2);
                    
                end
            end
        end
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % CALCULATION TYPE 1
    % DISPLACEMENT TIME 9
    % calculate velocity using +/- 4 time points
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    elseif displacement_time == 9
        
        for t = 5:time-5
            
            future_centroids4 = cell_centroid{t+4};
            future_centroids3 = cell_centroid{t+3};
            future_centroids2 = cell_centroid{t+2};
            future_centroids1 = cell_centroid{t+1};
            current_centroids = cell_centroid{t};
            previous_centroids = cell_centroid{t-1};
            previous_centroids2 = cell_centroid{t-2};
            previous_centroids3 = cell_centroid{t-3};
            previous_centroids4 = cell_centroid{t-4};
            
            % search through tracks at this time point
            for ii = 1:size(tracks,1)
                
                % store current cell
                future_cell4 = tracks(ii,t+4);
                future_cell3 = tracks(ii,t+3);
                future_cell2 = tracks(ii,t+2);
                future_cell1 = tracks(ii,t+1);
                current_cell = tracks(ii,t);
                previous_cell = tracks(ii,t-1);
                previous_cell2 = tracks(ii,t-2);
                previous_cell3 = tracks(ii,t-3);
                previous_cell4 = tracks(ii,t-4);
                
                % the cell exists for both the current and previous time points,
                % then we'll compute it's velocity
                if future_cell4 > 0 && future_cell3 > 0 && future_cell2 > 0 && ...
                        future_cell1 > 0 && current_cell > 0 && previous_cell > 0 && ...
                        previous_cell2 > 0 && previous_cell3 > 0 && previous_cell4 > 0
                    
                    deltaX1 = future_centroids4{future_cell4}(1) - future_centroids3{future_cell3}(1);
                    deltaY1 = future_centroids4{future_cell4}(2) - future_centroids3{future_cell3}(2);
                    
                    deltaX2 = future_centroids3{future_cell3}(1) - future_centroids2{future_cell2}(1);
                    deltaY2 = future_centroids3{future_cell3}(2) - future_centroids2{future_cell2}(2);
                    
                    deltaX3 = future_centroids2{future_cell2}(1) - future_centroids1{future_cell1}(1);
                    deltaY3 = future_centroids2{future_cell2}(2) - future_centroids1{future_cell1}(2);
                    
                    deltaX4 = future_centroids1{future_cell1}(1) - current_centroids{current_cell}(1);
                    deltaY4 = future_centroids1{future_cell1}(2) - current_centroids{current_cell}(2);
                    
                    deltaX5 = current_centroids{current_cell}(1) - previous_centroids{previous_cell}(1);
                    deltaY5 = current_centroids{current_cell}(2) - previous_centroids{previous_cell}(2);
                    
                    deltaX6 = previous_centroids{previous_cell}(1) - previous_centroids2{previous_cell2}(1);
                    deltaY6 = previous_centroids{previous_cell}(2) - previous_centroids2{previous_cell2}(2);
                    
                    deltaX7 = previous_centroids2{previous_cell2}(1) - previous_centroids3{previous_cell3}(1);
                    deltaY7 = previous_centroids2{previous_cell2}(2) - previous_centroids3{previous_cell3}(2);
                    
                    deltaX8 = previous_centroids3{previous_cell3}(1) - previous_centroids4{previous_cell4}(1);
                    deltaY8 = previous_centroids3{previous_cell3}(2) - previous_centroids4{previous_cell4}(2);
                    
                    % average the averages
                    deltaX = (deltaX1 + deltaX2 + deltaX3 + deltaX4 + deltaX5 + deltaX6 + deltaX7 + deltaX8)/8;
                    deltaY = (deltaY1 + deltaY2 + deltaY3 + deltaY4 + deltaY5 + deltaY6 + deltaY7 + deltaY8)/8;
                    
                    % store velocity components
                    CVX(ii,t) = deltaX;
                    CVY(ii,t) = deltaY;
                    
                    % create matrix storing centroids
                    CX(ii,t) = current_centroids{current_cell}(1);
                    CY(ii,t) = current_centroids{current_cell}(2);
                    
                end
            end
        end
    end

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALCULATION TYPE 2
% calculate centroid by looking at displacement across specified distance
% in time - i.e. if we say 9, then we'll simply look at the displacement
% from t-4 to t+4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif calculation_type == 2
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % CALCULATION TYPE 2
    % DISPLACEMENT TIME 3
    % calculate velocity using +/- 1 time point
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if displacement_time == 3
        
        for t = 2:time-2
            
            future_centroids1 = cell_centroid{t+1};
            current_centroids = cell_centroid{t};
            previous_centroids = cell_centroid{t-1};
            
            % search through tracks at this time point
            for ii = 1:size(tracks,1)
                
                % store current cell
                future_cell1 = tracks(ii,t+1);
                current_cell = tracks(ii,t);
                previous_cell = tracks(ii,t-1);
                
                % the cell exists for both the current and previous time points,
                % then we'll compute it's velocity
                if future_cell1 > 0 && current_cell > 0 && previous_cell > 0
                    
                    deltaX = future_centroids1{future_cell1}(1) - previous_centroids{previous_cell}(1);
                    deltaY = future_centroids1{future_cell1}(2) - previous_centroids{previous_cell}(2);
                    
                    % store velocity components
                    CVX(ii,t) = deltaX;
                    CVY(ii,t) = deltaY;
                    
                    % create matrix storing centroids
                    CX(ii,t) = current_centroids{current_cell}(1);
                    CY(ii,t) = current_centroids{current_cell}(2);
                    
                end
            end
        end
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % CALCULATION TYPE 2
    % DISPLACEMENT TIME 5
    % calculate velocity using +/- 2 time points
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    elseif displacement_time == 5
        
        for t = 3:time-3
            
            future_centroids1 = cell_centroid{t+2};
            current_centroids = cell_centroid{t};
            previous_centroids = cell_centroid{t-2};
            
            % search through tracks at this time point
            for ii = 1:size(tracks,1)
                
                % store current cell
                future_cell1 = tracks(ii,t+2);
                current_cell = tracks(ii,t);
                previous_cell = tracks(ii,t-2);
                
                % the cell exists for both the current and previous time points,
                % then we'll compute it's velocity
                if future_cell1 > 0 && current_cell > 0 && previous_cell > 0
                    
                    deltaX = future_centroids1{future_cell1}(1) - previous_centroids{previous_cell}(1);
                    deltaY = future_centroids1{future_cell1}(2) - previous_centroids{previous_cell}(2);
                    
                    % store velocity components
                    CVX(ii,t) = deltaX;
                    CVY(ii,t) = deltaY;
                    
                    % create matrix storing centroids
                    CX(ii,t) = current_centroids{current_cell}(1);
                    CY(ii,t) = current_centroids{current_cell}(2);
                    
                    
                end
            end
        end
        
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % CALCULATION TYPE 2
    % DISPLACEMENT TIME 7
    % calculate velocity using +/- 3 time points
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    elseif displacement_time == 7
        
        for t = 4:time-4
            
            future_centroids1 = cell_centroid{t+3};
            current_centroids = cell_centroid{t};
            previous_centroids = cell_centroid{t-3};
            
            % search through tracks at this time point
            for ii = 1:size(tracks,1)
                
                % store current cell
                future_cell1 = tracks(ii,t+3);
                current_cell = tracks(ii,t);
                previous_cell = tracks(ii,t-3);
                
                % the cell exists for both the current and previous time points,
                % then we'll compute it's velocity
                if future_cell1 > 0 && current_cell > 0 && previous_cell > 0
                    
                    deltaX = future_centroids1{future_cell1}(1) - previous_centroids{previous_cell}(1);
                    deltaY = future_centroids1{future_cell1}(2) - previous_centroids{previous_cell}(2);
                    
                    % store velocity components
                    CVX(ii,t) = deltaX;
                    CVY(ii,t) = deltaY;
                    
                    % create matrix storing centroids
                    CX(ii,t) = current_centroids{current_cell}(1);
                    CY(ii,t) = current_centroids{current_cell}(2);
                    
                end
            end
        end
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % CALCULATION TYPE 2
    % DISPLACEMENT TIME 9
    % calculate velocity using +/- 4 time points
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    elseif displacement_time == 9
        
        for t = 5:time-5
            
            future_centroids1 = cell_centroid{t+4};
            current_centroids = cell_centroid{t};
            previous_centroids = cell_centroid{t-4};
            
            % search through tracks at this time point
            for ii = 1:size(tracks,1)
                
                % store current cell
                future_cell1 = tracks(ii,t+4);
                current_cell = tracks(ii,t);
                previous_cell = tracks(ii,t-4);
                
                % the cell exists for both the current and previous time points,
                % then we'll compute it's velocity
                if future_cell1 > 0 && current_cell > 0 && previous_cell > 0
                    
                    deltaX = future_centroids1{future_cell1}(1) - previous_centroids{previous_cell}(1);
                    deltaY = future_centroids1{future_cell1}(2) - previous_centroids{previous_cell}(2);
                    
                    % store velocity components
                    CVX(ii,t) = deltaX;
                    CVY(ii,t) = deltaY;
                    
                    % create matrix storing centroids
                    CX(ii,t) = current_centroids{current_cell}(1);
                    CY(ii,t) = current_centroids{current_cell}(2);
                    
                end
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % CALCULATION TYPE 2
        % DISPLACEMENT TIME 11
        % calculate velocity using +/- 4 time points
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        elseif displacement_time == 11
        
        for t = 6:time-6
            
            future_centroids1 = cell_centroid{t+5};
            current_centroids = cell_centroid{t};
            previous_centroids = cell_centroid{t-5};
            
            % search through tracks at this time point
            for ii = 1:size(tracks,1)
                
                % store current cell
                future_cell1 = tracks(ii,t+5);
                current_cell = tracks(ii,t);
                previous_cell = tracks(ii,t-5);
                
                % the cell exists for both the current and previous time points,
                % then we'll compute it's velocity
                if future_cell1 > 0 && current_cell > 0 && previous_cell > 0
                    
                    deltaX = future_centroids1{future_cell1}(1) - previous_centroids{previous_cell}(1);
                    deltaY = future_centroids1{future_cell1}(2) - previous_centroids{previous_cell}(2);
                    
                    % store velocity components
                    CVX(ii,t) = deltaX;
                    CVY(ii,t) = deltaY;
                    
                    % create matrix storing centroids
                    CX(ii,t) = current_centroids{current_cell}(1);
                    CY(ii,t) = current_centroids{current_cell}(2);
                    
                end
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % CALCULATION TYPE 2
        % DISPLACEMENT TIME 13
        % calculate velocity using +/- 4 time points
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        elseif displacement_time == 13
        
        for t = 7:time-7
            
            future_centroids1 = cell_centroid{t+6};
            current_centroids = cell_centroid{t};
            previous_centroids = cell_centroid{t-6};
            
            % search through tracks at this time point
            for ii = 1:size(tracks,1)
                
                % store current cell
                future_cell1 = tracks(ii,t+6);
                current_cell = tracks(ii,t);
                previous_cell = tracks(ii,t-6);
                
                % the cell exists for both the current and previous time points,
                % then we'll compute it's velocity
                if future_cell1 > 0 && current_cell > 0 && previous_cell > 0
                    
                    deltaX = future_centroids1{future_cell1}(1) - previous_centroids{previous_cell}(1);
                    deltaY = future_centroids1{future_cell1}(2) - previous_centroids{previous_cell}(2);
                    
                    % store velocity components
                    CVX(ii,t) = deltaX;
                    CVY(ii,t) = deltaY;
                    
                    % create matrix storing centroids
                    CX(ii,t) = current_centroids{current_cell}(1);
                    CY(ii,t) = current_centroids{current_cell}(2);
                    
                end
            end
        end
    end
end

%%%%%
% set time ranges for calculations based on displacement_time
%%%%

if displacement_time == 3
    
    calculation_time_range = 2:time-2;
    
elseif displacement_time == 5
    
    calculation_time_range = 3:time-3;
    
elseif displacement_time == 7
    
    calculation_time_range = 4:time-4;
    
elseif displacement_time == 9
    
    calculation_time_range = 5:time-5;
    
elseif displacement_time == 11
    
    calculation_time_range = 6:time-6;
    
elseif displacement_time == 13
    
    calculation_time_range = 7:time-7;
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% store nonzero velocity readings, and the centroids corresponding to them,
% in cell arrays

% in order to quiver these plots, we need to get rid of all the zero
% entries

% additionally, we'll be propegating tracking information so that the
% centroids do not become divorced from the tracking infrastructure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize cell arrays [time, 2]
CVX2 = cell(size(tracks,2),1);
CVY2 = cell(size(tracks,2),1);
CX2 = cell(size(tracks,2),1);
CY2 = cell(size(tracks,2),1);
CA = cell(size(tracks,2),1);    % area


% iterate through time
for t = calculation_time_range
    
    count = 0;
    
    % iterate through tracks
    for ii = 1:size(tracks,1)
        
        % check to see if we have velocity readings for current track
        if CX(ii,t) > 0 %&& temp_CVY(ii,t) > 0
            
            % if so, update tracks
            count = count + 1;
            
            % store in cell array
            
            % velocity
            CVX2{t}{count}{1} = CVX(ii,t);
            CVX2{t}{count}{2} = ii;
            CVY2{t}{count}{1} = CVY(ii,t);
            CVY2{t}{count}{2} = ii;
            
            % centroid
            CX2{t}{count}{1} = CX(ii,t);
            CX2{t}{count}{2} = ii;
            CY2{t}{count}{1} = CY(ii,t);
            CY2{t}{count}{2} = ii;
            
            % area
            CA{t}{count}{1} = size(cell_area{t}{tracks(ii,t)},1);
            CA{t}{count}{2} = ii;
            
%             % store in cell array
%             CVX2{t}{count} = CVX(ii,t);
%             CVY2{t}{count} = CVY(ii,t);
%             CX2{t}{count} = CX(ii,t);
%             CY2{t}{count} = CY(ii,t);
            
        end
        
    end
    
end


% visual velocity field

% initialize structures to hold velocity information
x_velocities = cell(size(L,3),1);
y_velocities = cell(size(L,3),1);
x_coordinates = cell(size(L,3),1);
y_coordinates = cell(size(L,3),1);
vqX = cell(size(L,3),1);
vqY = cell(size(L,3),1);
x_query = cell(size(L,3),1);
y_query = cell(size(L,3),1);
interpArea = cell(size(L,3),1);
area = cell(size(L,3),1);

for t = calculation_time_range
    
    t
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % create vectors of the x/y components of the centroid positions and
    % velocity for the current time point (need to be pulled out of the
    % cell arrays)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % store number of velocities for this time point as a variable
    num_curr_velocities = length(CX2{t});
    
    % initialize vectors to store centroid coordinates and velocities
    temp_x_coordinates = zeros(num_curr_velocities,1);
    temp_y_coordinates = zeros(num_curr_velocities,1);
    temp_x_velocities = zeros(num_curr_velocities,1);
    temp_y_velocities = zeros(num_curr_velocities,1);
    temp_area = zeros(num_curr_velocities,1);
    
    % pull centroid coordinates and velocities for current time out of cell
    % arrays and put in vectors
    for ii = 1:length(temp_x_coordinates)
       temp_x_coordinates(ii) = CX2{t}{ii}{1};
       temp_y_coordinates(ii) = CY2{t}{ii}{1};
       temp_x_velocities(ii) = CVX2{t}{ii}{1};
       temp_y_velocities(ii) = CVY2{t}{ii}{1};
       temp_area(ii) = CA{t}{ii}{1};
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % interpolate centroid-based velocity vectors onto evenly spaced grid
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % calculate interpolation function
    FX = scatteredInterpolant(temp_x_coordinates,temp_y_coordinates,temp_x_velocities);
    FY = scatteredInterpolant(temp_x_coordinates,temp_y_coordinates,temp_y_velocities);
    FA = scatteredInterpolant(temp_x_coordinates,temp_y_coordinates,temp_area);
    
    % create vectors to define spacing at which we'll query the
    % interpolation function (x and y need to be flipped for some reason)
    x_query1 = linspace(x_min,x_max,length(x_min:x_max)/interp_subsample_factor);
    y_query1 = linspace(y_min,y_max,length(y_min:y_max)/interp_subsample_factor);
    
    % turn these vectors for interpolation spacing into meshgrids
    [temp_x_query,temp_y_query] = meshgrid(x_query1,y_query1);
    
    % pull out interpolated points along the mesh grid - these are
    % separated by x and y component, because so are the interpolation
    % functions
    temp_vqX = FX(temp_x_query,temp_y_query) .* velocity_scaling_factor;
    temp_vqY = FY(temp_x_query,temp_y_query) .* velocity_scaling_factor;
    temp_interpArea = FA(temp_x_query,temp_y_query);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % store data in cell array
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    x_velocities{t} = temp_x_velocities;
    y_velocities{t} = temp_y_velocities;
    x_coordinates{t} = temp_x_coordinates;
    y_coordinates{t} = temp_y_coordinates;
    vqX{t} = temp_vqX;
    vqY{t} = temp_vqY;
    x_query{t} = temp_x_query;
    y_query{t} = temp_y_query;
    area{t} = temp_area;
    interpArea{t} = temp_interpArea;
    
    
end

field1 = 'x_velocities'; value1 = x_velocities;
field2 = 'y_velocities'; value2 = y_velocities;
field3 = 'x_coordinates'; value3 = x_coordinates;
field4 = 'y_coordinates'; value4 = y_coordinates;
field5 = 'interp_x_velocities'; value5 = vqX;
field6 = 'interp_y_velocities'; value6 = vqY;
field7 = 'interp_x_coordinates'; value7 = x_query;
field8 = 'interp_y_coordinates'; value8 = y_query;
field9 = 'tracking_x_centroid'; value9 = CX2;
field10 = 'tracking_y_centroid'; value10 = CY2;
field11 = 'tracking_x_velocity'; value11 = CVX2;
field12 = 'tracking_y_velocity'; value12 = CVY2;
field13 = 'area'; value13 = area;
field14 = 'interpArea'; value14 = interpArea;

Velocity =  struct(field1,value1,field2,value2,field3,value3,field4,value4, ...
    field5,value5,field6,value6,field7,value7,field8,value8,field9,value9, ...
    field10,value10,field11,value11,field12,value12,field13,value13,...
    field14,value14);

end
