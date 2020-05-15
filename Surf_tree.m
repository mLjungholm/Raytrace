classdef Surf_tree < handle
%   Space partition object for acellerated ray tracing. Object takes in
%   points and triangle faces from sufaces and subpartition space fot them.
    
%   Created by Mikael Ljungholm - 2016-04
%   Point partition algorithm is based on original code by Sven Holcombe
%   1.2     - 2016-04 Modified to encompas surface triangles.
%   1.3     - 2016-05 VoxTree traversal algorithm added.
%   1.4     - 2018-06 Added code for multiple surfaces/volumes


    properties
        octree_properties;
        surf_nums;
        
        points;
        point_surf_index;
        point_bins;
        num_pts;
        
        faces;
        face_surf_index;
        face_n;
        face_count;
        
        surf_refract_index;
        surface_blocking;
        surf_p_count;
        surf_f_count;
        
        bin_count;
        bin_boundaries;
        bin_depths;
        bin_parents;
        bin_childs;
        bin_capacity = 3;
        
        bin_triangles;
        temp_tri_face_index;
        
        refract_order;
        surf_absorbing;
%         surface_absorbing
    end
    
    methods
        function this = Surf_tree(varargin)
            % This is the OcTree header line  
            this.surf_nums = nargin;
            tmp_p = cell(nargin,1);            %   points
            tmp_psi = cell(nargin,1);          %   point_surf_index
            tmp_spc = zeros(nargin,1,'uint16');%   surface_point_count
            tmp_f = cell(nargin,1);            %   faces
            tmp_fsi = cell(nargin,1);          %   face_surf_index
            tmp_sfc = zeros(nargin,1,'uint16');%   surface_face_count
            tmp_fn = cell(nargin,1);           %   surface normals
                       
            for i = 1:nargin
                tmp_p{i} = single(varargin{i}.v);               % store points form all surfaces in cell
                tmp_spc(i) = varargin{i}.v_count;               % Number of points for every surface
                tmp_psi{i} =  uint16(ones(tmp_spc(i),1)*i);     % Surface index for each point.
                
                tmp_f{i} = uint16(varargin{i}.f);               % Sorce faces in arrays
                tmp_sfc(i) = varargin{i}.f_count;               % Number of faces for each surface
                tmp_fsi{i} = uint16(ones(tmp_sfc(i),1)*i);      % Surface index for each face
                
                tmp_fn{i} = single(varargin{i}.n);
            end
                        
            % Modify the face indicies to the new point list
            p_count = tmp_spc(1);
            for i = 2:nargin
                tmp_f{i} = tmp_f{i} + p_count;
                p_count = p_count + tmp_spc(i);
            end
            
            % Store values
            this.points = cell2mat(tmp_p);
            this.faces = cell2mat(tmp_f);
            this.point_surf_index = cell2mat(tmp_psi);
            this.face_surf_index = cell2mat(tmp_fsi);
            this.surf_p_count = tmp_spc;
            this.surf_f_count = tmp_sfc;
            this.num_pts = sum(tmp_spc);
            this.face_n = cell2mat(tmp_fn);
            this.face_count = sum(tmp_sfc);

            % Initialise a single bin surrounding all given points
            this.bin_boundaries = [min(this.points,[],1) max(this.points,[],1)];
            this.point_bins = ones(this.num_pts,1,'uint16'); % All points start beloning to bin #1
            this.bin_depths = uint16(0); % Bin depth starts at 0
            this.bin_parents(1) = uint16(0); % Bin #1 hass no parent
            this.bin_count = uint16(1);
                        
            % Start dividing!
            this.preallocate_space;
            this.divide(1);
            this.deallocate_space;
            this.mk_childs;
            this.bin_parents = uint16(this.bin_parents);
            this.mk_tri_bin;
%             this.big_bin;
            
        end
        
        % MATLAB performs better if arrays that grow are initialised,
        % rather than grown during a loop. These two functions do just that
        % before and after the identification of new bins.
        function preallocate_space(this)
            num_bins = ceil(2*this.num_pts/this.bin_capacity);
            this.bin_depths(num_bins) = 0; % All bins start att 0
            this.bin_parents(num_bins) = 0;   % all binns belong to #0
            this.bin_boundaries(num_bins,1) = 0;  % all bins have boundary 0
            this.bin_childs(num_bins) = 0;  % No lower bins
        end
        function deallocate_space(this)
            this.bin_depths(this.bin_count+1:end) = [];
            this.bin_parents(this.bin_count+1:end) = [];
            this.bin_boundaries(this.bin_count+1:end,:) = [];
        end
        function divide(this, starting_bins)
            % Loop over each bin we will consider for division
            for i = 1:length(starting_bins)
                bin_no = starting_bins(i);
                % Prevent dividing beyond the maximum depth
%                 if this.BinDepths(binNo)+1 >= this.Properties.maxDepth
%                     continue;
%                 end
                
                % Prevent dividing beyond a minimum size                
%                 thisBounds = this.bin_boundaries(binNo,:);
%                 binEdgeSize = diff(thisBounds([1:3;4:6]));
%                 minEdgeSize = min(binEdgeSize);
%                 maxEdgeSize = max(binEdgeSize);
%                 if minEdgeSize < this.Properties.minSize
%                     continue;
%                 end
                
                % There are two conditions under which we should divide
                % this bin. 1: It's bigger than maxSize. 2: It contains
                % more points than binCapacity.
                old_count = this.bin_count;
                if nnz(this.point_bins == bin_no) > this.bin_capacity
                    this.divide_bin(bin_no);
                    this.divide(old_count+1:this.bin_count);
                    continue;
                end
%                 if maxEdgeSize>this.Properties.maxSize
%                     this.divideBin(bin_no);
%                     this.divide(old_count+1:this.bin_count);
%                     continue;
%                 end
            end
        end
        
        function divide_bin(this,bin_no)
            % Gather the new points (a bit more efficient to copy once)
            binPtMask = this.point_bins==bin_no;
            thisBinsPoints = this.points(binPtMask,:);
            
            % Get the old corner points and the new division point
            oldMin = this.bin_boundaries(bin_no,1:3);
            oldMax = this.bin_boundaries(bin_no,4:6);
%             if strcmp('weighted',this.Properties.style) && any(binPtMask)
%                 newDiv = mean(thisBinsPoints,1);
% else
                newDiv = mean([oldMin; oldMax], 1);
%             end
            
            % Build the new boundaries of our 8 subdivisions
            minMidMax = [oldMin newDiv oldMax];
            newBounds = minMidMax([...
                1 2 3 4 5 6;
                1 2 6 4 5 9;
                1 5 3 4 8 6;
                1 5 6 4 8 9;
                4 2 3 7 5 6;
                4 2 6 7 5 9;
                4 5 3 7 8 6;
                4 5 6 7 8 9]);
            
            % Determine to which of these 8 bins each current point belongs
            binMap = cat(3,[0 0 0],[0 0 1],[0 1 0],[0 1 1],...
                [1 0 0],[1 0 1],[1 1 0],[1 1 1]);
            gtMask = bsxfun(@gt, thisBinsPoints, newDiv);
            [~,binAssignment] = max(all(bsxfun(@eq,gtMask,binMap),2),[],3);
            % [~, binAssignment] = ismember(gtMask,binMap,'rows'); % A little slower than above.
            
            % Make the new bins and reassign old points to them
            newBinInds = this.bin_count+1:this.bin_count+8;
            this.bin_boundaries(newBinInds,:) = newBounds;
            this.bin_depths(newBinInds) = this.bin_depths(bin_no)+1;
            this.bin_parents(newBinInds) = bin_no;
            this.point_bins(binPtMask) = newBinInds(binAssignment);
            this.bin_count = this.bin_count + 8;
        end
        
        
        function binNos = query(this, newPts, queryDepth)
            % Get the OcTree bins that new query points belong to.
            %
            % BINS = OT.query(NEWPTS) searches the OcTree object OT and
            % returns an N-by-1 vector of BINS giving the bin index in
            % which each of the N points in NEWPTS is contained. For any
            % query points outside all bins in OT, the index -1 is
            % returned.
            %
            % BINS = OT.query(NEWPTS,DEPTH) restricts the search to DEPTH
            % levels in the OT bin tree. Note that the first bin
            % (containing all other bins in OT) has DEPTH = 1.

            if nargin<3
                queryDepth = max(this.bin_depths);
            end
            
            numPts = size(newPts,1);
            newPts = permute(newPts,[3 2 1]);
            binNos = ones(numPts,1)*-1;
                        
            bin_children = arrayfun(@(i)find(this.bin_parents==i),1:this.bin_count,'Un',0)';
            binIsLeaf = cellfun(@isempty, bin_children);
            ptQuery_recurse(1:numPts, this.bin_parents==0, 0)
            
            function ptQuery_recurse(newIndsToCheck_, binsToCheck, depth)
                % Build a list of all points that fall within one of the
                % bins to be checked, and the list of which point falls in
                % which bin.
                boundsToCheck = this.bin_boundaries(binsToCheck,:);
                [ptInBounds, subBinNo] = max(all(...
                    bsxfun(@ge, newPts(:,:,newIndsToCheck_), boundsToCheck(:,1:3)) & ...
                    bsxfun(@le, newPts(:,:,newIndsToCheck_), boundsToCheck(:,4:6))...
                    ,2),[],1);
            
                if ~all(ptInBounds)
                    % Special case usually when depth=0, where a point may
                    % fall outside the bins entirely. This should only
                    % happen once so let's fix it once and let subsequent
                    % code rely on all points being in bounds
                    binNos(newIndsToCheck_(~ptInBounds)) = -1;
                    newIndsToCheck_(~ptInBounds) = [];
                    subBinNo(~ptInBounds) = [];
                end
                binNosToAssign = binsToCheck(subBinNo);
                newIndsToAssign = newIndsToCheck_;
                binNos(newIndsToAssign) = binNosToAssign;
                
                % Allow a free exit when we reach a certain depth
                if depth>=queryDepth
                    return;
                end
                
                % Otherwise, for all of the points we just placed into
                % bins, check which of the children of those bins those
                % same points fall into
                [unqBinNos, ~, unqGrpNos] = unique(binNosToAssign);
                for i = 1:length(unqBinNos)
                    thisPtMask = unqGrpNos==i;
                    if ~binIsLeaf(unqBinNos(i))
                        ptQuery_recurse(newIndsToCheck_(thisPtMask), bin_children{unqBinNos(i)}, depth+1)
                    end
                end
                
            end
        end
        function plot_bin(this,bin_index)
            hold on;
            binMinMax = this.bin_boundaries(bin_index,:);
            pts = cat(1, binMinMax([...
                1 2 3; 4 2 3; 4 5 3; 1 5 3; 1 2 3;...
                1 2 6; 4 2 6; 4 5 6; 1 5 6; 1 2 6; 1 2 3]),...
                nan(1,3), binMinMax([4 2 3; 4 2 6]),...
                nan(1,3), binMinMax([4 5 3; 4 5 6]),...
                nan(1,3), binMinMax([1 5 3; 1 5 6]));
            plot3(pts(:,1),pts(:,2),pts(:,3),'r');
        end
        
        function h = plot(this,varargin)
            % OcTree.plot plots bin bounding boxes of an OcTree object
            %
            % H = OT.plot('name',value,...) allows you to specify any
            % properties of the bounding box lines that you would normally
            % supply to a plot(...,'name',value) command, and returns plot
            % object handles (one per bin) to H.
            hold on;
            h = zeros(this.bin_count,1);
            for i = 1:this.bin_count
                binMinMax = this.bin_boundaries(i,:);
                pts = cat(1, binMinMax([...
                    1 2 3; 4 2 3; 4 5 3; 1 5 3; 1 2 3;...
                    1 2 6; 4 2 6; 4 5 6; 1 5 6; 1 2 6; 1 2 3]),...
                    nan(1,3), binMinMax([4 2 3; 4 2 6]),...
                    nan(1,3), binMinMax([4 5 3; 4 5 6]),...
                    nan(1,3), binMinMax([1 5 3; 1 5 6]));
                h(i) = plot3(pts(:,1),pts(:,2),pts(:,3),varargin{:});
            end
        end
        function plot_surface(this,surface)
            figure(1)
            hold on
            axis equal
            n = surface;
            if isequal(surface,'all')
                n = 1:this.surf_nums;
            end
            for i = 1:size(n,2)
%                 rind = find(this.face_surf_index == i);
                tf = this.faces(this.face_surf_index == n(i),:);
                trisurf(tf,this.points(:,1),this.points(:,2),this.points(:,3),'Facecolor','b','FaceAlpha',0.1,'EdgeAlpha',0.3)
            end
        end     
        function plot_tri(this,tri_list)
            tf = this.faces(tri_list,:);
            %             for i = 1:size(n,2)
            %                 rind = find(this.face_surf_index == i);
            %                 tf = this.faces(this.face_surf_index == n(i),:);
            trisurf(tf,this.points(:,1),this.points(:,2),this.points(:,3),'Facecolor','r','FaceAlpha',0.1,'EdgeAlpha',0.3)
            %             end
        end
        function h = plot3(this,varargin)
            % OcTree.plot plots bin bounding boxes of an OcTree
            %
            % See also OcTree.plot
            h = this.plot(varargin{:});
        end
        
        function mk_childs(this)
            % creates a matrix indicating what if any binns are contained
            % in the current bin.
            tempChild = zeros(8,this.bin_count,'uint16');
            for i = 2:this.bin_count
                binNo = this.bin_parents(i);
                check = true;
                n = 1;
                while check
                    if tempChild(n,binNo) == 0
                        tempChild(n,binNo) = i;
                        break;
                    else
                        n = n +1;
                        if n == 9
                            break;
                        end
                    end
                end
            end
            this.bin_childs = tempChild;
        end
        
        function mk_tri_bin(this)
%             tempvar = arrayfun(@this.find_tri_bins,this.faces(:,1),this.faces(:,2),this.faces(:,3),'UniformOutput',false);
            this.temp_tri_face_index = cell2mat(arrayfun(@this.find_tri_bins,this.faces(:,1),this.faces(:,2),this.faces(:,3),'UniformOutput',false));
%             ttemp = arrayfun(@this.match_bin_with_tri,1:this.bin_count,'UniformOutput',false);
            this.bin_triangles = cell2mat(arrayfun(@this.match_bin_with_tri,1:this.bin_count,'UniformOutput',false));
            
%             this.bin_triangles = cell2mat(ttemp);
            this.temp_tri_face_index = [];
        end
        
        
        function found_list = find_tri_bins(this, face_ind1,face_ind2,face_ind3)
            % Function loops though all bins and checks if current
            % face/triangle matches
            found_list = uint16(zeros(1,this.bin_count));   % List of matched bins
            exclude_list = found_list;                      % List of bins known to miss
            current_bin = 2;                                % Index of current bin to test

            
            tri_vert = [this.points(face_ind1,:); 
                        this.points(face_ind2,:);
                        this.points(face_ind3,:)];
%             tri_vert = [this.point(this.face(face_ind,1),:); 
%                         this.point(this.face(face_ind,2),:);
%                         this.point(this.face(face_ind,3),:)];
            
            while current_bin <= this.bin_count
                % jump over aldready excluded bins
                if exclude_list(current_bin)
                    current_bin = current_bin + 1;
                    continue;
                end
                
%                 box_edges = this.bin_boundaries(current_bin,:); Not using
%                 this variable
                % Check if the current bin is a bottom bin.
                if nnz(this.bin_childs(:,current_bin))
                    % If it is a parent bin. check only bounding box
                    overlap = 1;
                    sort_tri = sort(tri_vert);
                    if sort_tri(1,1) > this.bin_boundaries(current_bin,4) || sort_tri(3,1) < this.bin_boundaries(current_bin,1)
                        overlap = 0;
                    elseif sort_tri(1,2) > this.bin_boundaries(current_bin,5) || sort_tri(3,2) < this.bin_boundaries(current_bin,2)
                        overlap = 0;
                    elseif sort_tri(1,3) > this.bin_boundaries(current_bin,6) || sort_tri(3,3) < this.bin_boundaries(current_bin,3)
                        overlap = 0;
                    end
                    
                    % If there is no overlap then add child bins to
                    % exclusion.
                    if ~overlap % cloud also use found_list for this. Ushure which is faster.
                        add_lower_bins(current_bin);
                        current_bin = current_bin + 1;
                        continue;
                    end 
                end
                
                % if it is a bottom bin -> do a full box_triangle_test.
                if ~nnz(this.bin_childs(:,current_bin))
                    found_list(current_bin) = BoxTriOverlap(this.bin_boundaries(current_bin,:), tri_vert);
                end
                current_bin = current_bin + 1;
            end
            
            % Function that itterates over all lower bins.
            function add_lower_bins(bin_no)
                lower_bin_list = this.bin_childs(:,bin_no);
                for i = 1:8
                    if nnz(this.bin_childs(:,lower_bin_list(i)))
                        add_lower_bins(lower_bin_list(i));
                    else
                        exclude_list(lower_bin_list(i)) = 1;
                    end
                end
            end
%             pad_list = uint16(zeros(10,1));
            temp_list = find(found_list);
            found_list = uint16(zeros(1,30));
            found_list(1:size(temp_list,2)) = temp_list;
        end
        function tri_list = match_bin_with_tri(this,bin_index)
            [row,~] = find(this.temp_tri_face_index == bin_index);
            tri_list = uint16(zeros(30,1));
            if ~isempty(row)
                tri_list(1:size(row,1)) = row;
            end
        end
        
    end
end


%    OcTree is used to create a tree data structure of bins containing 3D
%    points. Each bin may be recursively decomposed into 8 child bins.
%
%    OT = OcTree(PTS) creates an OcTree from an N-by-3 matrix of point
%    coordinates.
%
%    OT = OcTree(...,'PropertyName',VALUE,...) takes any of the following
%    property values:
%
%     binCapacity - Maximum number of points a bin may contain. If more
%                   points exist, the bin will be recursively subdivided.
%                   Defaults to ceil(numPts/10).
%     maxDepth    - Maximum number of times a bin may be subdivided.
%                   Defaults to INF.
%     maxSize     - Maximum size of a bin edge. If any dimension of a bin 
%                   exceeds maxSize, it will be recursively subdivided.
%                   Defaults to INF.
%     minSize     - Minimum size of a bin edge. Subdivision will stop after 
%                   any dimension of a bin gets smaller than minSize.
%                   Defaults to 1000*eps.
%
%    Example 1: Decompose 200 random points into bins of 20 points or less,
%             then display each bin with its points in a separate colour.
%        pts = (rand(200,3)-0.5).^2;
%        OT = OcTree(pts,'binCapacity',20);        
%        figure
%        boxH = OT.plot;
%        cols = lines(OT.bin_count);
%        doplot3 = @(p,varargin)plot3(p(:,1),p(:,2),p(:,3),varargin{:});
%        for i = 1:OT.bin_count
%            set(boxH(i),'Color',cols(i,:),'LineWidth', 1+OT.BinDepths(i))
%            doplot3(pts(OT.PointBins==i,:),'.','Color',cols(i,:))
%        end
%        axis image, view(3)
%
%
%
% OcTree methods:
%     shrink            - Shrink each bin to tightly encompass its children
%     query             - Ask which bins a new set of points belong to.
%     plot, plot3       - Plots bin bounding boxes to the current axes.
%
% OcTree properties:
%     Points            - The coordinate of points in the decomposition.
%     PointBins         - Indices of the bin that each point belongs to.
%     bin_count          - Total number of bins created.
%     bin_boundaries     - bin_count-by-6 [MIN MAX] coordinates of bin edges.
%     BinDepths         - The # of subdivisions to reach each bin.
%     bin_parents        - Indices of the bin that each bin belongs to.
%     Properties        - Name/Val pairs used for creation (see help above)
%     BinContent        - List of vetricis contained in each bin.
%     BinTri            - List of trianlge faces contained in each bin.
%     BinChilds         - List of child bins for each parent bin.
%     BinPrim           - List of primitives for each parent bin. 


% Possible additional info "BinContent" - idices of the points in each bin
