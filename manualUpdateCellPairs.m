function cellPairs = manualUpdateCellPairs(cellPairs, track_links)



for j = 2:size(track_links,1)
    
    % check which time point is larger
    if track_links(j,2) > track_links(j,6)
        test1 = 0;
    elseif track_links(j,2) < track_links(j,6)
        test1 = 1;
    end
    
    % pull out centroid components to ID cells
    if test1
        cell1_x = track_links(j,3);
        cell1_y = track_links(j,4);
        cell2_x = track_links(j,7);
        cell2_y = track_links(j,8);
        cell1_time = track_links(j,2);
        cell2_time = track_links(j,6);
    else
        cell1_x = track_links(j,7);
        cell1_y = track_links(j,8);
        cell2_x = track_links(j,3);
        cell2_y = track_links(j,4);
        cell1_time = track_links(j,6);
        cell2_time = track_links(j,2);
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % find cell IDs based on centroid
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % cell1
    for jj = 1:length(cell_centroid{cell1_time})-1
        track_links_centroid = [cell1_x cell1_y];
        centroid_query = cell_centroid{cell1_time}{jj};
        if centroid_query == track_links_centroid
            cell1_ID = jj;
        end
    end
    
    % cell2
    for jj = 1:length(cell_centroid{cell2_time})-1
        track_links_centroid = [cell2_x cell2_y];
        centroid_query = cell_centroid{cell2_time}{jj};
        if centroid_query == track_links_centroid
            cell2_ID = jj;
        end
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%
    % update track entries
    %%%%%%%%%%%%%%%%%%%%%%
    
    cellPairs{cell2_time}(cell2_ID) = cell1_ID;
    
    % clear all other cells that are linked to cell1_ID
    tempCellPairs = zeros(length(cellPairs{cell2_time}));
    for e = 1:length(cellPairs{cell2_time})
        if e ~= cell2_ID
            if cellPairs{cell2_time}(e) == cell1_ID
                cellPairs{cell2_time}(e) = 0;
            end
        end
    end
    
end

end