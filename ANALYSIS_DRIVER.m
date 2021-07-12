
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% analysis tools
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% There is example code below showing how to use the following variables
%
% 1) tracking
% 2) photorecepter IDs
% 3) column identity
% 4) row identity
% 5) R8 lattice



%% visualize tracking

% create colors for tracks
tracking_colors = [rand(size(tracks,1),1),rand(size(tracks,1),1),rand(size(tracks,1),1)];
display_tracks = uint8(zeros(size(L,1),size(L,2),3,size(L,3)));

% iterate through time
for t = 1:size(L,3)
    
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
                display_tracks(x_pixels(iii),y_pixels(iii),1,t) = tracking_colors(ii,1) * 255;
                display_tracks(x_pixels(iii),y_pixels(iii),2,t) = tracking_colors(ii,2) * 255;
                display_tracks(x_pixels(iii),y_pixels(iii),3,t) = tracking_colors(ii,3) * 255;
                
            end
            
        end
        
    end
    
    imshow(display_tracks(:,:,:,t))
    drawnow
    pause(0.1)
end

%% photorecepter IDs (R-cell IDs)

for t = 1:size(L,3)
    
    % binarize, dilate, 8bit seg-mask for display purposes.
    se = strel('disk',1);
    im = imdilate(not(L(:,:,t)),se);
    im = uint8(cat(3,im,im,im)) *255;
    
    % color in photoreceptors precluster cells
    for j = 1:length(preclusters)
        % iterate through precluster cells (always 5)
        for jj = 1:5
            % find unique ID of cell at current time point
            ID = tracks(preclusters{j}(jj),t);
            % make sure it exists at current time point (not off screen)
            if ID > 0
                % find area pixels
                pix = cell_area{t}{ID};
                % paint in cells according to subtype
                for k = 1:size(pix,1)
                    if jj == 1 %if R8
                        im(pix(k,1),pix(k,2),:) = [255 0 255];
                    elseif jj == 2  %if R2
                        im(pix(k,1),pix(k,2),:) = [200 0 200];
                    elseif jj == 3  %if R5
                        im(pix(k,1),pix(k,2),:) = [150 0 150];
                    elseif jj == 4  %if R3
                        im(pix(k,1),pix(k,2),:) = [100 0 100];
                    elseif jj == 5  %if R4
                        im(pix(k,1),pix(k,2),:) = [50 0 50];

                    end
                end
            end
        end
    end
    
    % print
    imshow(im)
    drawnow
    pause(0.1)
    
end

%% morphogenetic furrow location (anterior edge)

for t = 1:size(L,3)
    
    % binarize, dilate, 8bit seg-mask
    se = strel('disk',1);
    im = imdilate(not(L(:,:,t)),se);
    im = uint8(cat(3,im,im,im)) *255;
    
    % color in photoreceptors precluster cells
    for j = 1:length(preclusters)
        % iterate through precluster cells (always 5)
        for jj = 1:5
            % find unique ID of cell at current time point
            ID = tracks(preclusters{j}(jj),t);
            % make sure it exists at current time point (not off screen)
            if ID > 0
                % find area pixels
                pix = cell_area{t}{ID};
                % paint in cells according to subtype
                for k = 1:size(pix,1)
                    if jj == 1 %if R8
                        im(pix(k,1),pix(k,2),:) = [255 0 255];
                    elseif jj == 2  %if R2
                        im(pix(k,1),pix(k,2),:) = [200 0 200];
                    elseif jj == 3  %if R5
                        im(pix(k,1),pix(k,2),:) = [150 0 150];
                    elseif jj == 4  %if R3
                        im(pix(k,1),pix(k,2),:) = [100 0 100];
                    elseif jj == 5  %if R4
                        im(pix(k,1),pix(k,2),:) = [50 0 50];

                    end
                end
            end
        end
    end
    
    % print
    imshow(im)
    hold on
    
    % plot position of MF
    x = linspace(1,size(L,2),size(L,2));
    y = linspace(furrow_coordinate(t),furrow_coordinate(t),size(L,2));
    plot(x,y,'c-','LineWidth',3)
    
    hold off
    
    drawnow
    pause(0.1)
    
end

%% plot velocity (either centroid-based or interpolated) and MF position

for t = 1:size(L,3)
    
    % binarize, dilate, 8bit seg-mask for display purposes
    se = strel('disk',1);
    im = imdilate(not(L(:,:,t)),se);
    im = uint8(cat(3,im,im,im)) *100;
    
    % print
    imshow(im)
    hold on
    
    % plot position of MF
    x = linspace(1,size(L,2),size(L,2));
    y = linspace(furrow_coordinate(t),furrow_coordinate(t),size(L,2));
    plot(x,y,'r-','LineWidth',3)
    
    % plot centroid velocity
   quiver(Velocity(t).x_coordinates,Velocity(t).y_coordinates, ...
       Velocity(t).x_velocities,Velocity(t).y_velocities ,2,'c')
    
    % plot interpolated velocity
%     quiver(Velocity(t).interp_x_coordinates,Velocity(t).interp_y_coordinates, ...
%         Velocity(t).interp_x_velocities,Velocity(t).interp_y_velocities ,2,'c')

    drawnow
    
    pause(0.1)
    hold off
    
end


%% column identity

% ONLY WILDTYPE REPLICATE 1

for t = 1:size(L,3)
    
    L_columnDisplay = cat(3,L(:,:,t),L(:,:,t),L(:,:,t));
    colorz = lines(length(column_labels));
    for j = 1:length(column_labels)
        for jj = 1:length(column_labels{j})
            % check that R8 currently exists and then color in
            if tracks(preclusters{column_labels{j}(jj)}(1),t) > 0
                curr_area_R8 = cell_area{t}{tracks(preclusters{column_labels{j}(jj)}(1),t)};
                for jjj = 1:size(curr_area_R8,1)
                    L_columnDisplay(curr_area_R8(jjj,1),curr_area_R8(jjj,2),:) = colorz(j,:);
                end
            end
            % check that R2 currently exists and then color in
            if tracks(preclusters{column_labels{j}(jj)}(2),t) > 0
                curr_area_R2 = cell_area{t}{tracks(preclusters{column_labels{j}(jj)}(2),t)};
                for jjj = 1:size(curr_area_R2,1)
                    L_columnDisplay(curr_area_R2(jjj,1),curr_area_R2(jjj,2),:) = colorz(j,:);
                end
            end
            % check that R5 currently exists and then color in
            if tracks(preclusters{column_labels{j}(jj)}(3),t) > 0
                curr_area_R5 = cell_area{t}{tracks(preclusters{column_labels{j}(jj)}(3),t)};
                for jjj = 1:size(curr_area_R5,1)
                    L_columnDisplay(curr_area_R5(jjj,1),curr_area_R5(jjj,2),:) = colorz(j,:);
                end
            end
            % check that R3 currently exists and then color in
            if tracks(preclusters{column_labels{j}(jj)}(4),t) > 0
                curr_area_R3 = cell_area{t}{tracks(preclusters{column_labels{j}(jj)}(4),t)};
                for jjj = 1:size(curr_area_R3,1)
                    L_columnDisplay(curr_area_R3(jjj,1),curr_area_R3(jjj,2),:) = colorz(j,:);
                end
            end
            % check that R4 currently exists and then color in
            if tracks(preclusters{column_labels{j}(jj)}(5),t) > 0
                curr_area_R4 = cell_area{t}{tracks(preclusters{column_labels{j}(jj)}(5),t)};
                for jjj = 1:size(curr_area_R4,1)
                    L_columnDisplay(curr_area_R4(jjj,1),curr_area_R4(jjj,2),:) = colorz(j,:);
                end
            end
        end
    end
    imshow(L_columnDisplay)
    
end


%% row identity

% ONLY WILDTYPE REPLICATE 1

for t = 1:size(L,3)
    
    L_rowDisplay = cat(3,L(:,:,t),L(:,:,t),L(:,:,t));
    colorz = lines(length(row_labels));
    for j = 1:length(row_labels)
        for jj = 1:length(row_labels{j})
            % check that R8 currently exists and then color in
            if tracks(preclusters{row_labels{j}(jj)}(1),t) > 0
                curr_area_R8 = cell_area{t}{tracks(preclusters{row_labels{j}(jj)}(1),t)};
                for jjj = 1:size(curr_area_R8,1)
                    L_rowDisplay(curr_area_R8(jjj,1),curr_area_R8(jjj,2),:) = colorz(j,:);
                end
            end
            % check that R2 currently exists
            if tracks(preclusters{row_labels{j}(jj)}(2),t) > 0
                curr_area_R2 = cell_area{t}{tracks(preclusters{row_labels{j}(jj)}(2),t)};
                for jjj = 1:size(curr_area_R2,1)
                    L_rowDisplay(curr_area_R2(jjj,1),curr_area_R2(jjj,2),:) = colorz(j,:);
                end
            end
            % check that R5 currently exists
            if tracks(preclusters{row_labels{j}(jj)}(3),t) > 0
                curr_area_R5 = cell_area{t}{tracks(preclusters{row_labels{j}(jj)}(3),t)};
                for jjj = 1:size(curr_area_R5,1)
                    L_rowDisplay(curr_area_R5(jjj,1),curr_area_R5(jjj,2),:) = colorz(j,:);
                end
            end
            % check that R3 currently exists
            if tracks(preclusters{row_labels{j}(jj)}(4),t) > 0
                curr_area_R3 = cell_area{t}{tracks(preclusters{row_labels{j}(jj)}(4),t)};
                for jjj = 1:size(curr_area_R3,1)
                    L_rowDisplay(curr_area_R3(jjj,1),curr_area_R3(jjj,2),:) = colorz(j,:);
                end
            end
            % check that R4 currently exists
            if tracks(preclusters{row_labels{j}(jj)}(5),t) > 0
                curr_area_R4 = cell_area{t}{tracks(preclusters{row_labels{j}(jj)}(5),t)};
                for jjj = 1:size(curr_area_R4,1)
                    L_rowDisplay(curr_area_R4(jjj,1),curr_area_R4(jjj,2),:) = colorz(j,:);
                end
            end
        end
    end
    imshow(L_rowDisplay)
    
end


%% plot velocity, lattice, and MF position

% ONLY WILDTYPE REPLICATE 1

for t = 1:size(L,3)
    
    % binarize, dilate, 8bit seg-mask for display purposes
    se = strel('disk',1);
    im = imdilate(not(L(:,:,t)),se);
    im = uint8(cat(3,im,im,im)) *100;
    
    % print
    imshow(im)
    hold on
    
    % plot position of MF
    x = linspace(1,size(L,2),size(L,2));
    y = linspace(furrow_coordinate(t),furrow_coordinate(t),size(L,2));
    plot(x,y,'r-','LineWidth',3)
    
    % plot lattice
    patch('Faces',master_R8_lattice,'Vertices',lattice_points_overTime(:,:,t),...
        'FaceColor','none',...
        'EdgeColor','yellow',...
        'LineWidth',2)
    
    % plot centroid velocity
   quiver(Velocity(t).x_coordinates,Velocity(t).y_coordinates, ...
       Velocity(t).x_velocities,Velocity(t).y_velocities ,2,'c')
    
    % plot interpolated velocity
%     quiver(Velocity(t).interp_x_coordinates,Velocity(t).interp_y_coordinates, ...
%         Velocity(t).interp_x_velocities,Velocity(t).interp_y_velocities ,2,'c')

    drawnow
    
    pause(0.1)
    hold off
    
end


