classdef SurfaceOfInterest < diffgeometry.Manifold2D
    % SurfaceOfInterest extends Manifold, adding methods that go beyond the
    % basic mathematical structure
   
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
    % dependent properties
    %---------------------------------------------------------------------
    
    properties (SetAccess = protected, Dependent = true)
        
        data               % data pullback
    end
    
    %---------------------------------------------------------------------
    % properties
    %---------------------------------------------------------------------
    
    properties (SetAccess = protected)
        
        nLayers            % number of layers
    end
    
    %---------------------------------------------------------------------
    % public methods
    %---------------------------------------------------------------------
    
    methods
        
        %------------------------------------------------------
        % constructor
        %------------------------------------------------------
        
        function this = SurfaceOfInterest(varargin)
            % SURFACEOFINTEREST Create a new SOI
            % 
            % SurfaceOfInterest(dataSpace, embeddingSpace, timePoints, dynamic)
            % SurfaceOfInterest(datadir)
            %
            % embeddingSpace:   FiniteSet
            % dataSpace:        space of image data
            % timePoints:       vector of integer timepoints for which fields
            %                   should be defined
            % dynamic:          if false, the geometry is time independent
            % datadir:          directory holding a saved manifold
            
            if nargin == 4
                dataSpace = varargin{1};
                varargin = varargin(2:4);
            end
            
            this = this@diffgeometry.Manifold2D(varargin{:});
            
            this.nLayers = 0;
            
            % there should always be a data field
            if nargin == 1 && all(this.getField('data')==0)
                
                dataDir = varargin{1};
                DOMnode = xmlread(fullfile(dataDir, 'SOI.xml'));
                node = xpathNode(DOMnode, 'SOI/fields/Field[name=''data'']/targetSpace/Set');
                dataSpace = diffgeometry.FiniteSet(node);
                
                this.createField('data', 'diffgeometry.TensorPatch', dataSpace, true);
                
            elseif nargin == 4 
                
                this.createField('data', 'diffgeometry.TensorPatch', dataSpace, true);
            end
        end
        
        %------------------------------------------------------
        % pullbackStack
        %------------------------------------------------------
        
        function pullbackStack(this, stack, ROI, time, onionOpts)
            % PULLBACKSTACK Pull the stack data back to the surface
            %
            % pullbackStack(stack, ROI, time)
            %
            % stack:    Stack object containing the data
            % ROI:      RegionOfInterest object containing affine
            %           transformation from embedding coordinates to stack
            %           ROI can be left empty: ROI = []
            % time:     the time of the pullback, for storing in the SOI

            % We have:
            % - a chart phi: M -> R^2
            % - an embedding X: M -> R^3
            % - a transformation T : R^3 -> R^3'
            % - an image Im: R^3' -> I 
            % then produce  Im \circ T \circ X \circ \phi^{-1}: R^2 -> I
            % T is the inverse of the transformation we used to
            % transform the raw stack to the stack we used for surface
            % detection
            
            % time index for data
            ti = this.tIdx(time);
            
            % initialize onion options
            if nargin == 4
                onionOpts = struct('nLayers', 1, 'layerDistance', 1,...
                                            'sigma', 1, 'makeMIP', false);
            else
                assert(isfield(onionOpts, 'nLayers') && rem(onionOpts.nLayers, 2),...
                'onionOptions should have field nLayers with odd integer value')
                if ~isfield(onionOpts, 'layerDistance')
                    onionOpts.layerDistance = 1;
                end
                if ~isfield(onionOpts, 'sigma'); onionOpts.sigma = 1; end
                if ~isfield(onionOpts, 'makeMIP'); onionOpts.makeMIP = 0; end
            end

            % store number of layers
            this.nLayers = onionOpts.nLayers;
            halfLayers = (this.nLayers - 1)/2;

            % create fields for onion layers
            for li = 1:halfLayers;

                name = 'data_layer_';
                PT = this.data.patchClass;
                TS = this.data.targetSpace;
                    
                if this.getField([name 'p' num2str(li)]) == 0 
                    this.createField([name 'p' num2str(li)], PT, TS, true);
                    this.createField([name 'm' num2str(li)], PT, TS, true);
                end
            end
                    
            % time index for geometry: can be different because for a
            % static geometry the index should be one for every time but
            % the data will still be time dependent
            if this.dynamic == false,   gti = 1;
            else                        gti = this.tIdx(time); end
            
            %-------------------------------------------------------
            % first transform the embedding if necessary: T \circ X
            %-------------------------------------------------------
            
            if ~isempty(ROI)
                % TODO: more precise input testing
                                
                % define map from the aligned space to the stack space
                imageVolume = stack.image.domain;
                alignedImageVolume = this.embedding.targetSpace;
                T = ROI.translation;
                R = ROI.rotation;
                %TJitter = alignment.TJitter;
                stackTransform = diffgeometry.AffineMap(alignedImageVolume,...
                                                        imageVolume, inv(T)*inv(R));
                
                % transformed embedding, has imageVolume as target
                embedding = diffgeometry.Field(this.embedding(gti).name,...
                    this.embedding(gti).patchClass, imageVolume,...
                    this.embedding(gti).topology);
                
                for i = 1:length(this.embedding(gti).patches)
                    newPatch = stackTransform.compose(this.embedding(gti).patches{i});
                    embedding.addPatch(newPatch);
                end
            else
                embedding = this.embedding(gti);
            end
            
            %-------------------------------------------------------
            % then pull back the image in each chart
            %-------------------------------------------------------
            
            for i=1:length(this.atlas(gti).charts)
                                
                curChart = this.atlas(gti).charts{i};
                curDom = curChart.domain;

                debugMsg(2, ['Pulling back stack to chart ' curChart.image.name '\n']);
                
                % (T \circ X) \circ \phi^{-1}
                EmbPatch = embedding.getPatch(curDom.name);
                chartInv = curChart.getInverse();
                chartEmb = EmbPatch.compose(chartInv);
                                
                %---- mask interpolation outside data ---------------------

                % boundary curve: image of index set boundary
                ub = [curChart.apply{1}(1,:) curChart.apply{1}(:,end)'...
                    fliplr(curChart.apply{1}(end,:)) fliplr(curChart.apply{1}(:,1)')];
                vb = [curChart.apply{2}(1,:) curChart.apply{2}(:,end)'...
                    fliplr(curChart.apply{2}(end,:)) fliplr(curChart.apply{2}(:,1)')];
                
                umin = curChart.image.boundary{1}(1);
                vmin = curChart.image.boundary{2}(1);

                ub = (ub-umin)./curChart.image.stepSize(1);
                vb = (vb-vmin)./curChart.image.stepSize(2);

                usize = curChart.image.gridSize(1);
                vsize = curChart.image.gridSize(2);

                mask = poly2mask(double(ub), double(vb), vsize, usize);

                %---- create onion layers ---------------------
                
                projection = {};

                X = chartEmb.apply;
                dN = onionOpts.layerDistance;

                % compute normal
                [Nx,Ny,Nz] = GaussNormal(X{1},X{2},X{3},onionOpts.sigma);

                % normal displacement unit
                dX = {dN*Nx,dN*Ny,dN*Nz}; 
                for ni = 1:numel(dX); dX{ni}(isnan(dX{ni})) = 0; end

                % now loop through the layers
                for li = 1:this.nLayers

                    idx = li - halfLayers - 1;

                    % normally evolved embedding
                    def = {X{1} + idx*dX{1}, X{2} + idx*dX{2}, X{3} + idx*dX{3}};
                    embLayer = diffgeometry.CoordinateMap(chartEmb.domain,...
                                            chartEmb.image, def);

                    % Im \circ (T \circ X \circ \phi^{-1})
                    projection = stack.image.compose(embLayer);
                    grids = projection.apply;
                    
                    % go through channels to mask
                    for ci = 1:numel(grids)
                        grids{ci}(~mask) = NaN;
                    end

                    % to deal with multi-channel pullbacks
                    if numel(grids) > 1
                        type = [0]; % non-tensorial index: multiple scalars
                    else
                        type = []; % single scalar
                    end
                    
                    % convert index to field name
                    idx = li - halfLayers - 1;
                    if idx < 0
                        fieldName = ['data_layer_m' num2str(-idx)];
                    elseif idx > 0
                        fieldName = ['data_layer_p' num2str(idx)];
                    else
                        fieldName = 'data';
                    end

                    debugMsg(2, ['pulling back ' fieldName '\n']);
                    
                    dataField = this.getField(fieldName);
                    dataPatch = dataField(ti).getPatch(curChart.domain.name);

                    if dataPatch ~= 0

                        debugMsg(2, ['field already defined on ' curChart.domain.name...
                            ' adding it as transformation\n']);

                        dataPatch.setTransform(curChart.image.name, grids, type);

                    else 
                        dataPatch = diffgeometry.TensorPatch(curChart.domain,...
                            projection.image, grids, type,...
                            curChart.image.name, this.atlas(gti));

                        dataField(ti).addPatch(dataPatch);
                    end  
                end
            end
            
            % make MIP if wanted
            if onionOpts.makeMIP
                this.makeMIP(time);
            end
        end
        
        %------------------------------------------------------
        % make MIP
        %------------------------------------------------------
        
        function makeMIP(this, time)
            % make maximal intensity projection of multiplayer pullback
            %
            % makeMIP(time)
             
            assert(this.nLayers > 1, 'need multiple layers to make MIP');
            halfLayers = (this.nLayers - 1)/2;
            
            ti = this.tIdx(time);
            if this.dynamic == false,   gti = 1;
            else                        gti = ti; end
            
            % make field if necessary
            fieldName = 'data_MIP';
            if this.getField(fieldName) == 0 
                PT = this.data.patchClass;
                TS = this.data.targetSpace;
                this.createField(fieldName, PT, TS, true);
            end
            
            MIP = this.getField(fieldName);
            
            % now loop through charts and make mips
            for i=1:length(this.atlas(gti).charts)
                                
                curChart = this.atlas(gti).charts{i};
                domName = curChart.domain.name;
                
                MIPgrids = this.data(ti).getPatch(domName).apply;
                
                % loop through the layers
                for li = 1:this.nLayers

                    idx = li - halfLayers - 1;

                    % convert index to field name
                    if idx < 0
                        fieldName = ['data_layer_m' num2str(-idx)];
                    elseif idx > 0
                        fieldName = ['data_layer_p' num2str(idx)];
                    else
                        fieldName = 'data';
                    end

                    debugMsg(2, ['pulling back ' fieldName '\n']);

                    dataField = this.getField(fieldName);
                    dataPatch = dataField(ti).getPatch(domName);
                    grids = dataPatch.apply;
                    
                    for ci = 1:numel(MIPgrids)
                        MIPgrids{ci} = max(MIPgrids{ci}, grids{ci});
                    end
                end
                
                MIPpatch = MIP(ti).getPatch(domName);
                
                if MIPpatch ~= 0

                    debugMsg(2, ['field already defined on ' curChart.domain.name...
                        ' adding it as transformation\n']);

                    MIPpatch.setTransform(curChart.image.name, MIPgrids, dataPatch.type);

                else 
                    MIPpatch = diffgeometry.TensorPatch(dataPatch.domain,...
                        dataPatch.image, MIPgrids, dataPatch.type,...
                        curChart.image.name, this.atlas(gti));

                    MIP(ti).addPatch(MIPpatch);
                end  
            end
        end
        
        %------------------------------------------------------
        % measurements
        %------------------------------------------------------
        
        function lproper = properLength(this, time, curve, chartName)
            % calculate proper length of curve
            %
            % lproper = properLength(time, curve, chartName)
            %
            % time:         integer time point
            % curve:        Nx2 surface coordinates along curve
            % chartName:    chart in which curve is specified
            
            tidx = this.tIdx(time);
            
            u = curve;
            ushift = circshift(u, [1 0]);
            du = ushift - u;
            
            if isempty(this.g(tidx).patches)
                this.NCalcInducedMetric(time);
            end
            
            domName = this.atlas(tidx).getChart(chartName).domain.name;
            g = this.g(tidx).getPatch(domName).getTransform(chartName);
            gcurve = g.apply({u(:,1),u(:,2)});

            dlsq = du(:,1).*(gcurve{1,1}.*du(:,1) + gcurve{1,2}.*du(:,2))...
                    +du(:,2).*(gcurve{2,1}.*du(:,1) + gcurve{2,2}.*du(:,2));
            lproper = sum(sqrt(dlsq));
        end
        
        function Aproper = properArea(this, time, mask, chartName)
            % calculate proper area of region
            % 
            % Aproper = properArea(time, mask, chartName)
            %
            % time:         integer time point
            % curve:        binary mask with size matching chart
            % chartName:    chart in which mask is specified
            
            tidx = this.tIdx(time);
            
            domName = this.atlas(tidx).getChart(chartName).domain.name;
            g = this.g(tidx).getPatch(domName);
            
            detg = g.determinant().apply{1};
            Aproper = sum(detg(mask));
        end
        
        %------------------------------------------------------
        % load segmentation
        %------------------------------------------------------
        
        function loadSegmentation(this, dataDir, segDir)
            % LOADSEGMENTATION load segmented data into SOI
            %
            % loadSegmentation(dataDir, segDir)
            %
            % dataDir:  directory containing saved version of this SOI
            % segDir:   directory containing the segmentation, this should
            %           be structured as a the data subdirectory in fields,
            %           i.e. have subdirectories for patches etc with the
            %           same file names as the raw counterparts
            
            % basically reuses code from Manifold2D.load but loads the
            % segmented files with the data field metadata
            
            DOMnode = xmlread(fullfile(dataDir, 'SOI.xml'));
            fieldsNode = xpathNode(DOMnode, 'SOI/fields');
            fieldNode = fieldsNode.getFirstChild;

            while ~isempty(fieldNode)

                if  strcmp(fieldNode.getNodeName, 'Field') &&...
                    strcmp(char(xpathText(fieldNode, 'name')), 'data')

                    t = str2num(fieldNode.getAttribute('time'));
                    patchClass = char(xpathText(fieldNode, 'patchClass'));
                    
                    node = xpathNode(fieldNode, 'targetSpace/Set');
                    tspaceClass= str2func(char(node.getAttribute('subclass')));
                    targetSpace = tspaceClass(node);

                    % create the segmentation field
                    fDynamic = 1; 
                    this.createField('segmentation', patchClass, targetSpace, fDynamic);
                    segField = this.getField('segmentation');

                    % the atlas that needs to be passed to patch
                    % constructor
                    if this.dynamic, gti = this.tIdx(t);
                    else gti = 1; end
                    patlas = this.atlas(gti);
                    
                    % patches: contain the actual field data
                    patchesNode = xpathNode(fieldNode, 'patches');
                    patchNode = patchesNode.getFirstChild;
                    
                    while ~isempty(patchNode)

                        if strcmp(patchNode.getNodeName, 'Map')

                            patchConstructor = str2func(patchClass);
                            domName = xpathText(patchNode, 'domain/Set/name');
                            pdir = fullfile(segDir, domName);
                            patch = patchConstructor(patchNode, pdir, patlas);
                            
                            % if the data cannot be found, phi is left empty
                            if ~isempty(patch.phi)
                                segField(this.tIdx(t)).addPatch(patch);
                            end
                        end
                        
                        patchNode = patchNode.getNextSibling;
                    end
                end
                fieldNode = fieldNode.getNextSibling;
            end
        end
        
        %------------------------------------------------------
        % getters for dependent properties
        %------------------------------------------------------
        
        function data = get.data(this)
            % return the metric tensor
            data = this.getField('data');
        end
    end
 
end