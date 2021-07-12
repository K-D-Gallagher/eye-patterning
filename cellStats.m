function [cell_area, cell_centroid] = cellStats(L)

%time points
timePoints = size(L,3);

% Prealocate cell arrays for storing centroid and PixelIdxList
cell_area = cell(timePoints,1);
cell_centroid = cell(timePoints,1);
for ii = 1:timePoints
    cell_area{ii,1} = cell(length(unique(L(:,:,ii))),1);
    cell_centroid{ii,1} = cell(length(unique(L(:,:,ii))),1);
end


% computer cell area and centroid and store in cell arrays
for ii = 1:size(L,3)
    
    % regionprops
    stats = regionprops(L(:,:,ii),'centroid','PixelIdxList');
    
    % record in cell arrays
    for iii = 1:size(stats,1)
        cell_centroid{ii,1}{iii,1} = stats(iii).Centroid;
        [x_coord, y_coord] = ind2sub(size(L,1),stats(iii).PixelIdxList);
        cell_area{ii,1}{iii,1} = [x_coord,y_coord];
    end
    
end


end



