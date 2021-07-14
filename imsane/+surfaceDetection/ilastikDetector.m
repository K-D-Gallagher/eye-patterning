classdef ilastikDetector < surfaceDetection.surfaceDetector
    % detect surface by maximal radial derivative
    % options:  channel : channel of the stack to use
    %           sigma : width of the Gaussian 
    %           ssfactor : sub-sampling factor
    %           rmRadialOutliers : remove radial outlier >
    %               rmRadialOutliers*sigma (set zero for none)
    %           rmIntensityOutliers : remove edge intensity outliers
    %           thresh : threshold of intensity to be foreground
    %           amin : minimal connectec foreground pixel size
    %           bgdisc : radius of disc for background estimation
    %           dildisc : radius of dilation disc
    %           sp: stack dimension permutation     
    
    %---------------------------------------------------------------------
    % license
    %---------------------------------------------------------------------

    % Copyright 2014 Idse Heemskerk and Sebastian Streichan
    %
    % This file is part of ImSAnE.
    % 
    % ImSAnE is free software: you can redistribute it and/or modify
    % it under the terms of the GNU General Public License as published by
    % the Free Software Foundation, either version 3 of the License, or
    % (at your option) any later version.
    % 
    % ImSAnE is distributed in the hope that it will be useful,
    % but WITHOUT ANY WARRANTY; without even the implied warranty of
    % MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    % GNU General Public License for more details.
    % 
    % You should have received a copy of the GNU General Public License
    % along with ImSAnE.  If not, see <http://www.gnu.org/licenses/>.
    %
    % We make an explicit exception giving permission to link this code
    % with GPL-incompatible matlab facilities.
    
    %---------------------------------------------------------------------
    % properties
    %---------------------------------------------------------------------
    
    properties (Constant)
        % default detector options
        defaultOptions = struct('channel', 1, 'sigma', 2, 'ssfactor', 4,...
            'rmRadialOutliers', 2, 'rmIntensityOutliers', 2,...
            'thresh',25,'amin',20,'bgdisc',0,'dildisc',20,'sp',[1 2 3],'fileName',[]);      
    end
    
    %---------------------------------------------------------------------
    % public methods
    %---------------------------------------------------------------------    

    methods
         
        % ------------------------------------------------------
        % constructor
        % ------------------------------------------------------
        
        function this = ilastikDetector()
            % Constructor
            %
            % radialEdgeDetector()
            
            this = this@surfaceDetection.surfaceDetector();
        end 
        
        % ------------------------------------------------------
        % surface detection
        % ------------------------------------------------------
        
        function detectSurface(this, stack)
            % Detect surface in the stack with preset options.
            %
            % detectSurface(stack)
            
            % Detect the surface using cylindricalRadialD function from
            % generalFunctions/filters with specified options.
            % First down-sample the stack, then corect for background, 
            % detect edges, slice by slice compare intensity to thresh and
            % apply filters specified in the options. 
            % With resulting point cloud correct for radial outliers
            % and intensity outliers out. 
            
            opts = this.options;
            
            debugMsg(1, ['radialEdge.detectSurface() : channel='...
                num2str(opts.channel) ', sigma=' num2str(opts.sigma)...
                ', ssfactor=' num2str(opts.ssfactor)...
                ', rmRadialOutliers=' num2str(opts.rmRadialOutliers)...
                ', rmIntensityOutliers=' num2str(opts.rmIntensityOutliers)...
                ', threshold =' num2str(opts.thresh)...
                ', min area =' num2str(opts.amin)...
                ', background = ' num2str(opts.bgdisc) ...
                ', dilation disc =' num2str(opts.dildisc) '\n']);
            
            if ~isempty(stack)
                data = stack.image.apply{opts.channel};
            else
                error('stack is empty');
            end
            
            if sum( (opts.sp-[1 2 3]).^2 )~=0
                % permute the stack; 
                data = permute(data,opts.sp([2 1 3]));
            end
            
            if opts.ssfactor > 1
                data = data(1:opts.ssfactor:end, 1:opts.ssfactor:end, 1:opts.ssfactor:end);
            else
%                 if numel(data) > 10^8
%                     crazy = questdlg('Are you sure you want to call the detector with no sub-sampling?',...
%                         'Big data, no sub-sampling, potentially exceeding available memory, slow progress!',...
%                         'Yes','No','No');
%                     if strcmp(crazy, 'No')
%                         return;
%                     end
%                 end
            end
            
            
            %---------------------------------
            % edge detection 
            %---------------------------------

            % is based on ilastik results; 
            
            fileName = opts.fileName;
            
            file = h5read(fileName,'/exported_data');
            
            %data = permute(file,[4,3,2,1]);
            
            
            %file = h5read(name,'/volume/prediction/');
    pred = permute(file,[3,4,2,1]);
    pred = pred(:,:,:,2);
    pred = uint8(255*pred);
    
    
    zmin = 1;
            zmax = size(pred, 3)*4;
            
            ySize = size(pred,1)*4;
                xSize = size(pred,2)*4;
                zSize = size(pred,3)*4;

    %file = h5read(name2,'/volume/prediction/');
    %pred2 = permute(file,[4,3,2,1]);
    %pred2 = pred2(:,:,:,2);
    %pred2 = uint8(255*pred2); 
    %
    points = struct('x',[],'y',[],'z',[]);
    finMask = 0*pred;
    for z = 1 : size(pred,2)

        im = squeeze(pred(:,z,:));


        mask = bwareaopen(im > 255*4/8,200);
        mask = imfill(mask,'holes');


        %mask = bwconvhull(mask);
        mask([1,end],:) = 0;
        mask(:,[1,end]) = 0;
        mask = imdilate(mask,strel('disk',10));
        %if sum(mask(:))>0
        %
        lambda     = 10;   % Curvature parameter
        u_in       = 100;   % Foreground intensity 
        u_out      = 50;     % Background intensity 
        iterations = 500;     % Number of iterations
        seg = mask;

        seg = imerode(seg,strel('disk',10));

        finMask(:,z,:) = seg;
        per = bwperim(seg);
        if z > 1
        ind = find(per==1);
        else
            ind = find(seg == 1);
        end
        [I,J] = ind2sub(size(seg),ind);

        points(z).x = I;
        points(z).y = J;
        points(z).z = ones(size(I))*z;
            
       end
%close all;
    % 
    x = 4*cat(1,points.x);
    y = 4*cat(1,points.y);
    z = 4*cat(1,points.z);     
            
            points = [z,x,y];
            
           pointCloud = points;
            %---------------------------------
            % remove radial outliers
            %---------------------------------
            redPointCloud = [];
            if opts.rmRadialOutliers > 0
                
                for z = zmin:zmax
                    
                    pcSlice = pointCloud(pointCloud(:,3) == z,:);
                    CM = mean(pcSlice, 1);
                    Rsq = zeros([size(pcSlice,1) 1]);
                    for i = 1:3
                        Rsq = Rsq + (double(pcSlice(:,i)) - CM(i)).^2;
                    end
                    R = sqrt(Rsq);
                    good =  abs(R - mean(R)) < opts.rmRadialOutliers*std(R);
                    pcSlice = pcSlice(good, :);
                    redPointCloud = [redPointCloud;pcSlice];
                end
            end
            pointCloud = redPointCloud;
            
            %---------------------------------
            % scale point cloud to full size
            %---------------------------------
            
            largePC = zeros(size(pointCloud));
            for i = 1:3
                largePC(:,i) = (pointCloud(:,i)-1)*opts.ssfactor + 1;
            end
            
            largePC = sortrows(largePC, 3);
            if sum( (opts.sp-[1 2 3]).^2 )~=0
                % permute back; 
                inverseorder(opts.sp) = 1:numel(opts.sp);
                largePC = largePC(:,inverseorder([2 1 3]));
            end
            
            ROI = surfaceDetection.RegionOfInterest(eye(4),eye(4));
            % ranges
            xpRange = [1, xSize*opts.ssfactor];
            ypRange = [1, ySize*opts.ssfactor];
            zpRange = [1, zSize*opts.ssfactor];
            ROI.setRanges(xpRange,ypRange,zpRange);
            this.pointCloud = surfaceDetection.PointCloud(largePC,ROI);
            this.pointCloud.determineROI(15);
            %this.pointCloud = surfaceDetection.PointCloud(largePC);

        end
        
        function inspectQuality(this, inspectOpts, stack)
            %   inspect quality of fit in single slice in dimension specified 
            %   by options and display image.
            %
            %   inspectQuality(inspectOpts, stack)
            %
            %   override inspectQuality to deal with subsampled point cloud
            
            ssfac = this.options.ssfactor;
            
            if ssfac ~= 1 && rem(inspectOpts.value-1, ssfac) ~= 0
                inspectOpts.value = round(inspectOpts.value/ssfac)*ssfac + 1;
                disp([]);
                debugMsg(2, ['WARNING: sub-sampled point cloud taken from different'...
                    ' plane, may look a little off\n']);
            end

            inspectQuality@surfaceDetection.surfaceDetector(this, inspectOpts, stack);
        end
        
    end
end
