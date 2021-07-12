function [birth,death,total_born,total_dead,boundary_cells,birth_death_display] = birthDeath(L,tracks,cell_area)

% Function to detect cell births and deaths based on tracking information.
% In principle, if there are no errors in the tracking, then births and
% deaths should be detectable by the appearance and disappearance of cells.
% This will, of course, be tempered by boundary conditions.  We can't have
% cells being born and dying everytime they go on/off screen.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find boundary cells
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

boundary_cells = {};

for t = 1:size(L,3)
    
    % print time
    sprintf(strcat('Boundary Cells Time_', num2str(t)))
    
    %%%%%%%%%%%%%%%%%%%%%%%
    % find perimeter of FOV
    %%%%%%%%%%%%%%%%%%%%%%%
    tempT = imfill(imcomplement(L(:,:,t)));  % fill in seg mask
    se = strel('disk',2);
    tempT2 = imerode(tempT,se);              % erode it by 2 pixels
    tempT3 = tempT-tempT2;                   % subtract original mask from
    tempT4 = tempT3-tempT;                    % eroded mask to get perimeter
    
    % store perimeter as list of pixel indices
    temp_FOV_boundary = regionprops(tempT3,'PixelIdxList');
    [x_FOV_boundary,y_FOV_boundary] = ind2sub([size(L,1) size(L,2)],temp_FOV_boundary.PixelIdxList);
    
    %%%%%%%%%%%
    % record cell IDs that overlap with our perimeter pixel list
    %%%%%%%%%%%
    
    % initialize vector for storing the IDs of cells touching the perimeter
    matches = zeros(length(x_FOV_boundary),1);
    
    % search through perimeter pixel list
    for k = 2:length(x_FOV_boundary)
        % record the IDs of cells that overlap with the perimeter pixel list
        matches(k) = L(x_FOV_boundary(k),y_FOV_boundary(k),t);
    end
    
    % get rid of duplicate cell IDs
    matches = unique(matches);
    
    % trim off zero at beginning (from membrane overlap)
    matches = matches(2:end);
    
    % store vectors listing boundary cells in cell array
    boundary_cells{t} = matches;
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find cell births and deaths
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize cell-arrays to store birth/death info
temp_births = cell(size(L,3),1);
temp_deaths = cell(size(L,3),1);

% iterate through tracks over time
for t = 1:size(tracks,2)-1  % end one time point early because we can't
    % calculate birth/death for all time points
    
    % print time
    sprintf(strcat('Cell Birth/Death Time_', num2str(t)))
    
    % pull out this info, just because it's easier to index
    future_tracks = tracks(:,t+1);
    current_tracks = tracks(:,t);
    count_1 = 0;
    count_2 = 0;
    
    % go through all tracks per time point
    for k = 1:size(tracks,1)
        
        % make sure this isn't a boundary cell (don't wanna be counting
        % going on/off screen as births and deaths)
        if not(ismember(future_tracks(k),boundary_cells{t+1})) && ...
                not(ismember(current_tracks(k),boundary_cells{t}))
            
            % if a cell exists at the future time point, but not the
            % current, then a cell is born
            if future_tracks(k) ~= 0 && current_tracks(k) == 0
                count_1 = count_1 + 1;
                temp_births{t+1}{count_1} = future_tracks(k);
                
                % if a cell does not exist in the future time point, but does in
                % the current, then it has died
            elseif future_tracks(k) == 0 && current_tracks(k) ~= 0
                count_2 = count_2 + 1;
                temp_deaths{t}{count_2} = current_tracks(k);
                
            end
            
        end
        
    end
    
end

% store as vectors inside cell array
birth = cell(size(L,3),1);
death = cell(size(L,3),1);

for t = 1:size(L,3)
    
    temp = zeros(length(temp_births{t}),1);
    for k = 1:length(temp_births{t})
        temp(k) = temp_births{t}{k};
    end
    birth{t} = temp;
    
    temp = zeros(length(temp_deaths{t}),1);
    for k = 1:length(temp_deaths{t})
        temp(k) = temp_deaths{t}{k};
    end
    death{t} = temp;
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% collect all tracks of cells born and cells that die
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% find tracks of all cells that were born or that died
total_born = [];
total_dead = [];
count1 = 0;
count2 = 0;
for t = 1:size(L,3)
    % find tracks of cells born at this time point
    if not(isempty(birth{t}))
        for j = 1:length(birth{t})
            count1 = count1 + 1;
            [~, LOCB] = ismember(birth{t}(j),tracks(:,t));
            total_born(count1) = LOCB;
        end
    end
    % find tracks of cells that die at this time point
    if not(isempty(death{t}))
        for j = 1:length(death{t})
            count2 = count2 + 1;
            [~, LOCB] = ismember(death{t}(j),tracks(:,t));
            total_dead(count2) = LOCB;
        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% labels a cell at the time before it dies or the time after it is born
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


birth_death_display = uint8(zeros(size(L,1),size(L,2),3,size(L,3)));
for t = 1:size(L,3)
    
    % print time
    sprintf(strcat('Display Generation Time_', num2str(t)))
    
    % look through all cells at current time point
    for k = 2:length(cell_area{t})-1
        
        % check if cell is born
        if ismember(k,birth{t})
            
            % track?
            [~,track_ID] = ismember(k,tracks(:,t));
            
            
            %%%%
            % t
            %%%%
            curr_ID = tracks(track_ID,t);
            
            % find if exists
            if curr_ID > 0
                
                temp = cell_area{t}{curr_ID};
                
                for kk = 1:size(temp,1)
                    birth_death_display(temp(kk,1),temp(kk,2),2,t) = 255;
                    birth_death_display(temp(kk,1),temp(kk,2),3,t) = 255;
                end
                
            end
            
            

            
        elseif ismember(k,death{t})
            
            % track?
            [~,track_ID] = ismember(k,tracks(:,t));
            
            
            %%%
            % t
            %%%
            curr_ID = tracks(track_ID,t);
            
            % find if exists
            if curr_ID > 0
                
                temp = cell_area{t}{curr_ID};
                
                for kk = 1:size(temp,1)
                    birth_death_display(temp(kk,1),temp(kk,2),1,t) = 255;
                    birth_death_display(temp(kk,1),temp(kk,2),3,t) = 255;
                end
                
            end
            
            
        end
        
        
    end
    
    
end






end









