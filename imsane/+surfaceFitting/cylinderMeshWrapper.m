classdef cylinderMeshWrapper < surfaceFitting.meshWrapper
    % Generate surface from a point cloud using thin plate spline fit
    %
    % The thin plate spline method fits a surface that behaves like a plate
    % with some bending stiffness. This surface remains parametrized by
    % coordinates xy on a plane, which we are also calling planar, so this
    % method should be used to fit point clouds that correspond to planar
    % surfaces.
       
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
   
    
    %---------------------------------------------------------------------
    % public methods
    %---------------------------------------------------------------------
    
    methods
        
        % ------------------------------------------------------
        % constructor
        % ------------------------------------------------------
        
        function this = cylinderMeshWrapper()
            % MESHWRAPPER Create a SOI from a generic mesh
            %
            % meshWrapper()
            %

            
            % call superclass constructor
            
            this = this@surfaceFitting.meshWrapper();
            
            % overload charts, since second level child classes run into
            % issues with instantiating the abstract property charts; We
            % may want to not declare an abstract property in the
            % superclass surfaceFitter and therfore instead overload in all
            % children - This is a do to. 
            
            % charts - charts this fitter can produce and their properties
            %
            % the rows have the structure: 
            % {name, description, stepsize, intersections, fundamental, desired}
            %
            % stepsize:         2d vector of stepsize
            % intersections:    lists indices of charts with overlapping domain
            % fundamental:      boolean, does chart define a set in the topology
            % desired:          boolean indicating whether to produce it
            
            this.charts = struct(...
            'name', {'cylinder1', 'cylinder2',...
                        'cylinder1_proper', 'cylinder2_proper', ...
                        'anteriorEquidistant',      'posteriorEquidistant'},...
            'description', {'Cylinder coordinates z, phi',...
                            'Cylinder shifted by pi',...
                            'Proper cylinder coordinates',...
                            'Shifted proper cylinder coordinates',...
                            'Proper coordinates in lower z planes',...
                            'Proper coordinates in upper z planes',},...
            'stepSize', {[1 .1], [1 .1], [1 1], [1 1], [1 1], [1 1]},...
            'intersections', {[],[],[],[],[],[]},...%{[2 3 4], [1 3 4], [], [], [], []},...
            'fundamental', {1, 1, 0, 0, 1 ,1},...
            'desired', {1, 1, 1, 1, 0, 0});
            
            % initialize fitOptions
            this.fitOptions.chartSeeds      = [];     
            this.fitOptions.transitionWidth = 100;
            this.fitOptions.fixAxis         = 1;
            this.fitOptions.rotation        = eye(4);
            this.fitOptions.translation     = eye(4);
            this.fitOptions.fixResolution   = 0;
            this.fitOptions.resolution      = [];
            
            % initialize fittedParam
            this.fittedParam = struct('mesh', [], 'submeshes', [],...
                                      'submVidx', [], 'maxDist',[],...
                                      'intersectIdx',{cell(0)},'intersects',[]);
        end
        
        
        % ------------------------------------------------------
        %  generate embedding
        % ------------------------------------------------------
               
        function [embeddings, charts] = generateEmbedding(this, chartName)
            % GENERATEEMBEDDING Create the embedding in some chart
            %
            % [embedding, chart] = generateEmbedding(chartName)
            %
            % The embedding is generated using the result of a fit, stored
            % in fittedParam, so fitSurface should be called first.
            % The available charts are listed in the charts property of the
            % SurfaceFitter object.
            % 
            % See also surfaceFitting.tpsFitter.fitSurface
            
            % number of seeds
            %seeds = this.fitOptions.chartSeeds;
            %nSeeds = numel(seeds);
            
            desCharts = struct();
            for i = 1:length(this.charts)
                desCharts.(this.charts(i).name) = this.charts(i).desired;
            end
            
            embeddings = {};
            charts = {};
            
            % --------------------------------------------
            % Define Cylinder 1 & 2 Embedding
            % --------------------------------------------
            
            this.fittedParam.submeshes{1}.u = cell(1,sum(cat(1,this.charts.desired)));
            
            name  = 'cylinder1'; 
            if strcmp(chartName, name) 
            
                % here, we parametrize the mesh using a cylindrical
                % coordinate system. For this we orient the data such that
                % the long axis coincides with z; 
                
                % extract the point cloud to determine its orientation. 
                
                points = this.fittedParam.mesh.v;
                
                
                if all([sum(sum((this.fitOptions.rotation-eye(4)).^2)),...
                        sum(sum((this.fitOptions.translation-eye(4)).^2))]==0) || ...
                        this.fitOptions.fixAxis == 0
                    
                    debugMsg(2, 'Updating orientation\n');
                    
                    pc = surfaceDetection.PointCloud(points);
                    pc.determineROI(5);
                    rotation    = pc.ROI.rotation;
                    translation = pc.ROI.translation;

                    this.fitOptions.rotation    = rotation;
                    this.fitOptions.translation = translation;
                end
                
                % from setRoi in the PointCloud class ;
                R = this.fitOptions.rotation;
                T = this.fitOptions.translation;

                
                pointsHomogenous = [points ones([size(points,1) 1])];
                alignedPoints = (R*T*pointsHomogenous')';
                alignedPoints = alignedPoints(:,1:3);
                
    
                z = alignedPoints(:,3); %  u
                phi = atan2(alignedPoints(:,2),alignedPoints(:,1))+pi; % v
                r   = sqrt( (alignedPoints(:,1)).^2+(alignedPoints(:,2)).^2);
                
                
                this.fittedParam.submeshes{1}.u{1} = [z(:),phi(:)];
                
                % determine step size in phi;
                r_max = max(r);
                %dz   = 1; % this will be the user fed input. 
                %dphi = dz/r_max; 
                
                c1step = this.getChartResolution('cylinder1');
                
                if isempty(this.fitOptions.resolution)
                    dz   = c1step(1);
                    dphi = dz/r_max;
                    this.fitOptions.resolution = [dz dphi];
                else
                    dz   = this.fitOptions.resolution(1);
                    dphi = this.fitOptions.resolution(2);
                end
                % define the stepsize in the charts; 
                stepSize = [dz dphi];
                this.setChartResolution('cylinder1', stepSize);
                this.setChartResolution('cylinder2', stepSize);
                
                % define the boundary of the future chart;
%                boundary = {[min(round(z)), max(round(z))], [0, 2*pi]};
%                 boundary = {[min(round(z)), max(round(z))],...
%                             [min(phi) max(phi)]};
                boundary = {[min(round(z)), max(round(z))],[pi/2 3*pi/2]};
            
                        
                [zG,phiG] = meshgrid(boundary{1}(1):dz:boundary{1}(2),...
                    boundary{2}(1):dphi:boundary{2}(2));
                % phiG is the phi coordinate in cylinder1; 
                x = alignedPoints(:,1);
                y = alignedPoints(:,2);
                   
                
                % finally we are in the position to generate the embedding
                % grids for x and y; 
                % we multiply everywhere with r_max to get phi into real
                % distances, rather than angles; Otherweise the matlab
                % internal routine to compute the triangulation produces
                % errors.
                if exist('scatteredInterpolant') > 0
                    F = scatteredInterpolant(z,r_max*phi,x,'linear','none');
                else 
                    F = TriScatteredInterp(z,r_max*phi,x);
                end
                xG = F(zG,r_max*phiG);

                if exist('scatteredInterpolant') > 0
                    F = scatteredInterpolant(z,r_max*phi,y,'linear','none');
                else 
                    F = TriScatteredInterp(z,r_max*phi,y);
                end
                yG = F(zG,r_max*phiG);
                
                % don't forget to rotate [xG,yG,zG] back into the embedding space!
            
                grid = [xG(:),yG(:),zG(:)];
                
                gridHomogenous = [grid ones([size(grid,1) 1])];
                alignedGrid = (inv(T)*inv(R)*gridHomogenous')';
                alignedGrid = alignedGrid(:,1:3);
                
                % this is the embeddding grid
                grids{1} = reshape(alignedGrid(:,1),size(xG,1),size(xG,2));
                grids{2} = reshape(alignedGrid(:,2),size(xG,1),size(xG,2));
                grids{3} = reshape(alignedGrid(:,3),size(xG,1),size(xG,2));
                
                
                % and now we deal with the ImSAnE infrastructure, i.e. make 
                % charts, embedding etc. ...
                
                % define the image and domain of the chart
                % boundary was defined above; so was stepsize
                chartim  = diffgeometry.FiniteSet(chartName,...
                                                boundary, stepSize);
                domain   = chartim.makeIndexSet();

                % chart: name_index -> name
                % the domain is the image of the chart, the image is the domain of
                % the fit, i.e. the embedding space
                chart = diffgeometry.CoordinateMap(domain, chartim,...
                                                chartim.makeHandles());

                % for charts made with FiniteSet.makeHandles we can also set an
                % analytic inverse (for speed and accuracy)
                chartInv = diffgeometry.CoordinateMap(chartim, domain,...
                                        chartim.makeInverseHandles());
                chart.setInverse(chartInv);

                % embedding: cylinder1 -> targetSpace
                embedding = diffgeometry.CoordinateMap(chartim,...
                                         this.fitDomain, grids);

                % embedding \circ chart: cylinder1_index -> targetSpace
                embedding = embedding.compose(chart); 

                %
                embeddings{1} = embedding;
                charts{1} = chart;
                
                
                
                
                % --------------------------------------------
                % Define Cylinder2 embedding
                % --------------------------------------------

                if desCharts.('cylinder2') == 1    
                
                name  = 'cylinder2'; 
                % shift phi and repeat the above steps for cylinder2
                % charts; 
                
                % u remained z, just v = phi2 changes!
                phi2 = mod(phi + (phi<pi)*pi + (phi>pi)*pi,2*pi); % v,

                this.fittedParam.submeshes{1}.u{2} = [z(:),phi2(:)];
                
                boundary = {[min(round(z)), max(round(z))], [0, 2*pi]};
                
                [zG,phiG] = meshgrid(boundary{1}(1):dz:boundary{1}(2),...
                    boundary{2}(1):dphi:boundary{2}(2));
                % phiG is the phi coordinate in cylinder1; 
                x = alignedPoints(:,1);
                y = alignedPoints(:,2);
                   
                
                % finally we are in the position to generate the embedding
                % grids for x and y; 
                % we multiply everywhere with r_max to get phi into real
                % distances, rather than angles; Otherweise the matlab
                % internal routine to compute the triangulation produces
                % errors.
                if exist('scatteredInterpolant') > 0
                    F = scatteredInterpolant(z,r_max*phi2,x);
                else 
                    F = TriScatteredInterp(z,r_max*phi2,x);
                end
                xG = F(zG,r_max*phiG);

                if exist('scatteredInterpolant') > 0
                    F = scatteredInterpolant(z,r_max*phi2,y);
                else 
                    F = TriScatteredInterp(z,r_max*phi2,y);
                end
                yG = F(zG,r_max*phiG);
                
                % don't forget to rotate [xG,yG,zG] back into the embedding space!
            
                grid = [xG(:),yG(:),zG(:)];
                
                gridHomogenous = [grid ones([size(grid,1) 1])];
                alignedGrid = (inv(T)*inv(R)*gridHomogenous')';
                alignedGrid = alignedGrid(:,1:3);
                
                % this is the embeddding grid
                grids{1} = reshape(alignedGrid(:,1),size(xG,1),size(xG,2));
                grids{2} = reshape(alignedGrid(:,2),size(xG,1),size(xG,2));
                grids{3} = reshape(alignedGrid(:,3),size(xG,1),size(xG,2));
                
                
                % and now we deal with the ImSAnE infrastructure, i.e. make 
                % charts, embedding etc. ...
                
                % define the image and domain of the chart
                % boundary was defined above; so was stepsize
                chartim  = diffgeometry.FiniteSet(name,...
                                                boundary, stepSize);
                domain   = chartim.makeIndexSet();

                % chart: name_index -> name
                % the domain is the image of the chart, the image is the domain of
                % the fit, i.e. the embedding space
                chart = diffgeometry.CoordinateMap(domain, chartim,...
                                                chartim.makeHandles());

                % for charts made with FiniteSet.makeHandles we can also set an
                % analytic inverse (for speed and accuracy)
                chartInv = diffgeometry.CoordinateMap(chartim, domain,...
                                        chartim.makeInverseHandles());
                chart.setInverse(chartInv);

                % embedding: cylinder1 -> targetSpace
                embedding = diffgeometry.CoordinateMap(chartim,...
                                         this.fitDomain, grids);

                % embedding \circ chart: cylinder1_index -> targetSpace
                embedding = embedding.compose(chart); 

                %
                embeddings{2} = embedding;
                charts{2} = chart;
                
                
                disp('Currently no transition map between the cylinders implemented!')
                
                end
                
                %[tmap12, tmap21] = this.makeTMaps(2, 1, charts{2}, charts{1});
                
                %SOI.atlas(gti).addTransitionMap(tmap12);
                %SOI.atlas(gti).addTransitionMap(tmap21);
                
            end
            
            
        end
        
        
        
        % ------------------------------------------------------
        %  populate SOI
        % ------------------------------------------------------
        
        function  populateSOI(this, SOI, varargin)
            % POPULATESOI Add fit result to SurfaceOfInterest object
            %
            % populateSOI(SOI)
            % populateSOI(SOI, t)
            %
            % SOI:  SurfaceOfInterest object
            % t:    time, needs to be provided if SOI.dynamic = true
            %
            % Generates chart domain, chart and embedding using the result 
            % of a fit, adds these to the SOI.
            % Fit results are stored in fittedParam, so fitSurface should 
            % be called first.
            % 
            % See also surfaceFitting.tpsFitter.fitSurface
            
            % gti : geometric time index (one for static geometry, time index for
            % dynamic geometry)
            if SOI.dynamic == false
                gti = 1;
                if length(varargin) == 1
                    debugMsg(1, 'static SOI: ignoring time provided\n');
                end
            else
                if length(varargin) == 1
                    gti = SOI.tIdx(varargin{1});
                else
                    error('for dynamic SOI, time argument needs to be provided');
                end
            end
            
            desCharts = struct();
            for i = 1:length(this.charts)
                desCharts.(this.charts(i).name) = this.charts(i).desired;
            end

            %----------------------------------------------------------
            % Define charts
            %----------------------------------------------------------

            % define the image of the chart
%             for j = 1:numel(this.charts)
%                 
%                 name = this.charts(j).name;
%                 
%                 if desCharts.(name) == 1
                    
                    name = 'cylinder1';
                    debugMsg(2, ['generating ' name ' charts \n']);
                    
                    [embeddings, charts] = this.generateEmbedding(name);

                    % TODO: topology = overlap, transition maps
                    for i = 1:numel(embeddings)

                        intersects = {};
                        SOI.topologicalSpace(gti).addSet(charts{i}.domain, intersects);
                        SOI.atlas(gti).addChart(charts{i});
                        SOI.embedding(gti).addPatch(embeddings{i});
                    end
%                 end
%             end
        
        
        

            %-----------------------------------------------------------------
            % Calculate induced metric
            %-----------------------------------------------------------------

            % after having added the embeddings in different charts we can now
            % calculate the induced metric
            
            % here we need to 
            SOI.NCalcInducedMetric(SOI.timePoints(gti));

            %-----------------------------------------------------------------
            % Proper charts cylinder1
            %-----------------------------------------------------------------

            origChartNames = {};

            if desCharts.('cylinder1_proper') == 1 
                origChartNames = [origChartNames, 'cylinder1'];
            end

            if desCharts.('cylinder2_proper') == 1    
                origChartNames = [origChartNames, 'cylinder2'];
            end
            
            
            
            
            
            if ~isempty(origChartNames)
        
                for i = 1:length(origChartNames)

                    debugMsg(2, ['generating ' origChartNames{i} ' chart\n']);

                    origChart = SOI.atlas(gti).getChart(origChartNames{i});
                    domain = origChart.domain;

                    dz = origChart.image.stepSize(1);
                    df = origChart.image.stepSize(2);
                    
                    gzz = SOI.g(gti).getPatch(domain.name).cmp({1,1});
                    gff = SOI.g(gti).getPatch(domain.name).cmp({2,2});

                    % remove nans; 
                    gzz(isnan(gzz)) = 0;
                    gff(isnan(gff)) = 0;
                    
                    % there is some numerical issues on the edges, where
                    % charts are anyways not designed for. -> remove
                    gzz(1:9,:) = 0;
                    gzz(end-9:end,:) = 0;
                    gzz(:,1:9) = 0;
                    gzz(:,end-9:end) = 0;
                    gff(1:9,:) = 0;
                    gff(end-9:end,:) = 0;
                    gff(:,1:9) = 0;
                    gff(:,end-9:end) = 0;
                    
                    % zp = int_0^z dz sqrt(gzz) 
                    % z is first coordinate so second index
                    zp = real(cumsum(sqrt(gzz)*dz, 2)); 

                    % non-monotonous behaviour of z' as a function of z for phip ==
                    % const lead to a re-definiton of z': z' = int
                    % sqrt(gzz(z,phi_with_max_distance)) dz;
                    % Identify the phip with maximum distance zp;
                    [~,maxind] =  max(zp(:));
                    [II,~] = ind2sub(size(zp),maxind);
                    % generate a grid that replicates this zp; 
                    [zp,~] = meshgrid(zp(II,:),1:size(zp,1));

                    % fp = int_0^f df sqrt(gff) + c(z) (see egggeometry)
                    intdphi = real(cumsum(sqrt(gff)*df, 1));
                    mid = round(size(gff,1)/2); % mid ~ pi
                    cz = pi - real(sum(sqrt(gff(1:mid,:))*df, 1));
                    cz = repmat(cz, size(intdphi,1), 1);
                    fp = intdphi + cz; 

                    this.fittedParam.submeshes{1}.u{2+i} = [zp(:),fp(:)];
                    
                    
                    
                    
                    % use the same stepsize for the proper chart image
                    bdry = {[min(zp(:)), max(zp(:))], [min(fp(:)), max(fp(:))]};
                    
                    image = diffgeometry.FiniteSet([origChartNames{i} '_proper'], bdry, [dz dz]);
                    pumpkin = diffgeometry.CoordinateMap(domain, image, {zp, fp});

                    SOI.atlas(gti).addChart(pumpkin);

                    %-----------------------------------------------------------------
                    % Proper charts at poles
                    %-----------------------------------------------------------------
                    % These charts are based on the proper distances we get through
                    % the integration of the metric during generation of 
                    % cylinder1_proper charts. To avoid recomputation, we put the
                    % generation of anteriorEquidistant and posteriorEquidistant
                    % charts in this position. 
                    if i == 1 

                        if desCharts.('anteriorEquidistant') == 1

                            
                             disp('Then implement anterior Equidistant charts!')
                            
%                             %%%%%
%                             % anteriorEquidistant charts;         
%                             %%%%%                        
% 
%                             debugMsg(2, 'generating anteriorEquidistant chart\n');
% 
% 
%                             % define a new chart, that contains the z' as z and phi as phi; 
%                             phi      = origChart.apply{2};
%                             z        = origChart.apply{1};
% 
%                             % perform a tri scattered interpolation to obtain z(zp,phi)
%                             % and so on.
% 
%                             z_zp_phi = TriScatteredInterp(double(zp(:)),double(phi(:)),double(z(:)));
%                             eGrids   = SOI.embedding.getPatch(origChart.domain.name).apply();
%                             x_zp_phi = TriScatteredInterp(double(zp(:)),double(phi(:)),double(eGrids{1}(:)));
%                             y_zp_phi = TriScatteredInterp(double(zp(:)),double(phi(:)),double(eGrids{2}(:)));
% 
% 
%                             % now define a meshgrid xG,yG, compute rG and phiG (raduis and angle)
%                             % in this meshgrid and plug thisinto z_zp_phi, thus relating
%                             % zp with r_G
% 
%                             fac = .6;% defines reach of charts in fraction of max(zp). .6 is a good value 
%                             % for overlap with other charts.
% 
%                             [xG,yG] = meshgrid(-max(zp(:))*fac:dz:max(zp(:))*fac);
%                             % get r and phi
%                             [phiG,rG] = cart2pol(xG,yG);
%                             phiG = phiG+pi;
%                             phiG(phiG>max(phi(:))) = max(phi(:));
%                             zG = z_zp_phi(double(rG),double(phiG));
% 
%                             xGE = x_zp_phi(double(rG),double(phiG));
%                             yGE = y_zp_phi(double(rG),double(phiG));
% 
%                             clear z_zp_phi;
%                             clear x_zp_phi;
%                             clear y_yp_phi;
%                             % this is the definition of the embedding.
%                             zG  = griddata(double(xG(~isnan(zG))),double(yG(~isnan(zG))),double(zG(~isnan(zG))),double(xG),double(yG));
%                             xGE = griddata(double(xG(~isnan(xGE))),double(yG(~isnan(xGE))),double(xGE(~isnan(xGE))),double(xG),double(yG));
%                             yGE = griddata(double(xG(~isnan(yGE))),double(yG(~isnan(yGE))),double(yGE(~isnan(yGE))),double(xG),double(yG));
% 
%                             zG(isnan(zG)) = 0;
% 
%                             bdry = {[1, size(zG,1)], [1, size(zG,2)]};
%                             domain = diffgeometry.FiniteSet('anteriorEquidistant_index', bdry, [1 1]);
% 
%                             bdry = {[-max(zp(:))*fac, max(zp(:))*fac], [-max(zp(:))*fac, max(zp(:))*fac]};
%                             image = diffgeometry.FiniteSet('anteriorEquidistant', bdry, [1 1]);
%                             chartPLZ = diffgeometry.CoordinateMap(domain, image, {xG, yG});
% 
% 
%                             %bdry = {[min(z(:)), max(z(:))], [0, 2*pi]};
%                             %image = diffgeometry.FiniteSet(['polarLowerZ2_proper'], bdry, [dz df]);
%                             %properPole = diffgeometry.CoordinateMap(domain, image, {zG, phiG});
% 
% 
%                             bdry      = {[1, size(zG,1)], [1, size(zG,2)]};
%                             domain    = diffgeometry.FiniteSet('anteriorEquidistant_index', bdry, [1 1]);
% 
%                             bdry      = ({[min(xG(:)),max(xG(:))],[min(yG(:)),max(yG(:))]});
%                             image     = diffgeometry.FiniteSet('anteriorEquidistant',bdry,[1 1]);
%                             chart     = diffgeometry.CoordinateMap(domain, image, {xG, yG});
% 
% 
%                             embeddingPLZ = diffgeometry.CoordinateMap(chart.domain, this.fitDomain,...
%                                          {xGE, yGE, zG});
% 
%                             intersects = {};
%                             SOI.topologicalSpace(gti).addSet(chartPLZ.domain, intersects);
%                             SOI.atlas(gti).addChart(chartPLZ);
%                             SOI.embedding(gti).addPatch(embeddingPLZ);         
% 
% 
%                             % generate transition map between polarLowerZ2 and cylinder1    
%                             % this takes x,y to z, phi
%                             % z(x,y) and, phi(x,y) still have to be
%                             % computed
% 
%                             % the embedding is defined as
%                             % X = R cos(phi) + X0, Y = R sin(phi)/sqrt(1-e^2) + Y0
%                             cosphi = xGE - X0(zG);
%                             sinphi = sqrt( 1- e(zG).^2 ).*( yGE - Y0(zG) );
%                             % this phi is the angle, the next phi is a transition map!
%                             phiPLZ = atan2(sinphi, cosphi);
% 
%                             % there seems to be an offset between this phi and the phi
%                             % of the cylinder1 chart. We solve this by determination of
%                             % a phi in cylinder1, and finding the best match in the
%                             % two embeddings. The difference to phiPLZ will be
%                             % considered an offset. 
%                             bla = chart1.apply();
%                             ff = bla{2}(round(end/2),100);
% 
%                             embedding1Grids = embedding1.apply();
% 
%                             diff = (xGE-embedding1Grids{1}(round(end/2),100)).^2+(yGE-embedding1Grids{2}(round(end/2),100)).^2+...
%                                     (zG-embedding1Grids{3}(round(end/2),100)).^2;
%                             [~,ii] = min(diff(:));
% 
%                             phiPLZ = mod(phiPLZ,2*pi);
%                             phioff = ff - phiPLZ(ii);
% 
%                             phiPLZ = mod(phiPLZ+phioff,2*pi);
% 
%                             % tmapLZ1 : polarLowerZ (x,y) -> cylinder1 (z, phi)
%                             phi_PLZ1 = { zG, phiPLZ};           
%                             tmapPLZ1   = diffgeometry.CoordinateMap(chartPLZ.image, chart1.image, phi_PLZ1); 
% 
%                             % now generate the transition map from cylinder1 to polarLowerZ2,
%                             % By composing pumkin with inverse cylinder, we get a map from z into zp, 
%                             % then the relations of the x,y plane are simply zp*cos(phi). 
% 
%                             bla = origChart.apply();
%                             tmap11P = pumpkin.composeInverse(chart1);
%                             bla2 = tmap11P.apply();
%                             x = (bla2{1}).*cos(bla{2});
%                             y = (bla2{1}).*sin(bla{2});
% 
% 
%                             % tmapLZ1 : cylinder1 -> polarLowerZ
%                             phi_1PLZ = {x,y};
%                             tmap1PPLZ   = diffgeometry.CoordinateMap(origChart.image, chartPLZ.image, phi_1PLZ); 
% 
%                             % now obtain the transition map between cylinder and
%                             % pumpkin. Then compose the result with the transition map
%                             % from pumpkin to PLZ to get a transition map from cylinder
%                             % to PLZ.
% 
%                             % curChart is the pumpkin1 Chart, chart1 the cylinder1
%                             %tmap11P = pumpkin.composeInverse(chart1);
%                             %tmap1PLZ = tmap1PPLZ.compose(tmap11P);
% 
% 
%                             % add both transition maps to the atlas
%                             SOI.atlas(gti).addTransitionMap(tmapPLZ1);
%                             SOI.atlas(gti).addTransitionMap(tmap1PPLZ);
% 
%                             % mow define transition map between cylinder2 and polarLowerZ2 
%                             % this is much shorter, since we can use the tmap from cylinder1
%                             % to cylinder2. 
% 
%                             % tmapLZ2: polarLowerZ -> cylinder2
%                             tmapPLZ2   = tmap12.compose(tmapPLZ1);
%                             % tmap2LZ: cylinder2 -> polarLowerZ
%                             tmap2PLZ   = tmap1PPLZ.compose(tmap21);
% 
%                             % add to the atlas
%                             SOI.atlas(gti).addTransitionMap(tmapPLZ2);
%                             SOI.atlas(gti).addTransitionMap(tmap2PLZ);

                        end


                        %%%%%
                        % posteriorEquidistant charts;         
                        %%%%% 


                        if desCharts.('posteriorEquidistant') == 1      

                            disp('Then implement posterior Equidistant charts!')
%                             debugMsg(2, 'generating posteriorEquidistant chart\n');
% 
%                             zp = zp(:,end:-1:1);
%                             % we now formally define a chart with an inverted proper z.
%                             % this won't be added to the atlas. 
%                             bdry = {[min(zp(:)), max(zp(:))], [min(fp(:)), max(fp(:))]};
% 
%                             image = diffgeometry.FiniteSet([origChartNames{i} '_proper_inv'], bdry, [dz dz]);
%                             pumpkin_semiInv = diffgeometry.CoordinateMap(origChart.domain, image, {zp, fp});
% 
% 
%                             z_zp_phi = TriScatteredInterp(double(zp(:)),double(phi(:)),double(z(:)));
%                             x_zp_phi = TriScatteredInterp(double(zp(:)),double(phi(:)),double(eGrids{1}(:)));
%                             y_zp_phi = TriScatteredInterp(double(zp(:)),double(phi(:)),double(eGrids{2}(:)));
% 
% 
%                             % now define a meshgrid xG,yG, compute rG and phiG in this meshgrid and plug this 
%                             % into z_zp_phi
%                             [xG,yG] = meshgrid(-max(zp(:))*fac:dz:max(zp(:))*fac);
%                             % get r and phi
%                             [phiG,rG] = cart2pol(xG,yG);
%                             phiG = phiG+pi;
%                             phiG(phiG>max(phi(:))) = max(phi(:));
%                             zG = z_zp_phi(double(rG),double(phiG));
% 
%                             xGE = x_zp_phi(double(rG),double(phiG));
%                             yGE = y_zp_phi(double(rG),double(phiG));
% 
%                             clear z_zp_phi;
%                             clear x_zp_phi;
%                             clear y_yp_phi;
%                             zG  = griddata(double(xG(~isnan(zG))),double(yG(~isnan(zG))),double(zG(~isnan(zG))),double(xG),double(yG));
%                             xGE = griddata(double(xG(~isnan(xGE))),double(yG(~isnan(xGE))),double(xGE(~isnan(xGE))),double(xG),double(yG));
%                             yGE = griddata(double(xG(~isnan(yGE))),double(yG(~isnan(yGE))),double(yGE(~isnan(yGE))),double(xG),double(yG));
% 
%                             zG(isnan(zG)) = 0;
% 
%                             bdry = {[1, size(zG,1)], [1, size(zG,2)]};
%                             domain = diffgeometry.FiniteSet('posteriorEquidistant_index', bdry, [1 1]);
% 
%                             bdry = {[-max(zp(:))*fac, max(zp(:))*fac], [-max(zp(:))*fac, max(zp(:))*fac]};
%                             image = diffgeometry.FiniteSet('posteriorEquidistant', bdry, [dz dz]);
%                             chartPUZ = diffgeometry.CoordinateMap(domain, image, {xG, yG});
% 
%                             bdry      = {[1, size(zG,1)], [1, size(zG,2)]};
%                             domain    = diffgeometry.FiniteSet('posteriorEquidistant_index', bdry, [1 1]);
% 
%                             bdry      = ({[min(xG(:)),max(xG(:))],[min(yG(:)),max(yG(:))]});
%                             image     = diffgeometry.FiniteSet('posteriorEquidistant',bdry,[1 1]);
%                             chart     = diffgeometry.CoordinateMap(domain, image, {xG, yG});
% 
% 
%                             embeddingPUZ = diffgeometry.CoordinateMap(chart.domain, this.fitDomain,...
%                                          {xGE, yGE, zG});   
% 
% 
%                             intersects = {};
%                             SOI.topologicalSpace(gti).addSet(chartPUZ.domain, intersects);
%                             SOI.atlas(gti).addChart(chartPUZ);
%                             SOI.embedding(gti).addPatch(embeddingPUZ);         
% 
% 
% 
%                             % Transition maps
%                             cosphi = xGE - X0(zG);
%                             sinphi = sqrt( 1- e(zG).^2 ).*( yGE - Y0(zG) );
%                             % this phi is the angle, the next phi is a transition map!
%                             phiPUZ = atan2(sinphi, cosphi);
%                             phiPUZ = mod(phiPUZ,2*pi);
%                             phi_PUZ1 = { zG, phiPUZ};           
%                             tmapPUZ1   = diffgeometry.CoordinateMap(chartPUZ.image, chart1.image, phi_PUZ1); 
% 
%                              % now generate the transition map from cylinder1 to polarLowerZ2,
%                             % By composing pumkin with inverse cylinder, we get a map from z into zp, 
%                             % then the relations of the x,y plane are simply zp*cos(phi). 
% 
%                             bla = origChart.apply();
%                             tmap11P = pumpkin_semiInv.composeInverse(chart1);
%                             bla2 = tmap11P.apply();
%                             x = (bla2{1}).*cos(bla{2});
%                             y = (bla2{1}).*sin(bla{2});
% 
% 
%                             % tmapUZ1 : cylinder1 -> polarUpperZ
%                             phi_1PUZ = {x,y};
%                             tmap1PPUZ   = diffgeometry.CoordinateMap(origChart.image, chartPUZ.image, phi_1PUZ); 
% 
%                             SOI.atlas(gti).addTransitionMap(tmapPUZ1);
%                             SOI.atlas(gti).addTransitionMap(tmap1PPUZ);
% 
%                             % mow define transition map between cylinder2 and polarUpperZ2 
%                             % this is much shorter, since we can use the tmap from cylinder1
%                             % to cylinder2. 
% 
%                             % tmapUZ2: polarUpperZ -> cylinder2
%                             tmapPUZ2   = tmap12.compose(tmapPUZ1);
%                             % tmap2UZ: cylinder2 -> polarLowerZ
%                             tmap2PUZ   = tmap1PPUZ.compose(tmap21);
% 
%                             % add to the atlas
%                             SOI.atlas(gti).addTransitionMap(tmapPUZ2);
%                             SOI.atlas(gti).addTransitionMap(tmap2PUZ);
                        end
                    end
                end
            end
            
            
        
        end
        
    end
end