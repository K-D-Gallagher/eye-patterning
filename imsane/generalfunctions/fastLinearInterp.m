function result = fastLinearInterp(inputGrids,image,ss)


%[X,Y] = meshgrid(0:100);

%Z = X;
%ss = 2;

%ssX = X(1:ss:end,1:ss:end);
%ssY = Y(1:ss:end,1:ss:end);

%ssZ = Z(1:ss:end,1:ss:end);

%
%inputGrids = {X,Y};
%ssGrids    = {ssX,ssY};
%image      = uint16(ssZ);

% 

imageInterp = zeros(size(image,1)*ss,size(image,2)*ss,'single');

divSum = imageInterp+1;
imageInterp(1:2:end,1:2:end)   = image;

% 
% interpolate colums;
imageInterp(1:2:end,2:2:end)   = imageInterp(1:2:end,2:2:end)+image;
imageInterp(1:2:end,2:2:end-1) = imageInterp(1:2:end,2:2:end-1)+image(:,2:end);
divSum(1:2:end,2:2:end-1)      = divSum(1:2:end,2:2:end-1) +1;

% interpolate rows; 
imageInterp(2:2:end,1:2:end)   = imageInterp(2:2:end,1:2:end)+image;
imageInterp(2:2:end-1,1:2:end) = imageInterp(2:2:end-1,1:2:end)+image(2:end,:);
divSum(2:2:end-1,1:2:end)      = divSum(2:2:end-1,1:2:end) +1;

% final step: interpolate row & colum; 
imageInterp(2:2:end,2:2:end)     = imageInterp(2:2:end,2:2:end)+image;
%
imageInterp(2:2:end-1,2:2:end-1) = imageInterp(2:2:end-1,2:2:end-1)+image(2:end,2:end);
divSum(2:2:end-1,2:2:end-1)      = divSum(2:2:end-1,2:2:end-1) +1;


% 
result = imageInterp./divSum;

result = result(1:size(inputGrids{1},1),1:size(inputGrids{1},2));

%surf(double(result)-Z), shading interp
end