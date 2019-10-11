classdef OcTree2 < handle
%   Space partition object for acellerated ray tracing. Object takes in
%   points and triangle faces from sufaces and subpartition space fot them.
    
%   Created by Mikael Ljungholm - 2016-04
%   Point partition algorithm is based on original code by Sven Holcombe
%   1.2     - 2016-04 Modified to encompas surface triangles.
%   1.3     - 2016-05 VoxTree traversal algorithm added.
%   1.4     - 2018-06 Added code for multiple surfaces/volumes


    properties
        octree_properties;
        
        points;
        point_surf_index;
        point_bins;
        
        faces;
        faces_surf_index;
        
        surf_refractive_index;
        surface_blocking;

        absorption_volumes;
        volume_abs_coefficient;
        volumes_total_absorption;
        
        bin_count;
        bin_boundaries;
        bin_depths;
        bin_parents;
        bin_childs;
        
        bin_triangles;
        bin_absorption_volumes;
    end
    
    methods
        function this = OcTree(type, varargin)
            % This is the OcTree header line
%             validateattributes(pts,{'numeric'},...
%                 {'real','finite','nonnan','ncols', 3},...
%                 mfilename,'PTS')
            switch type
                case 'surface'
                    
                case 'abs_volume'
                    
                otherwise
                    fprintf(strcat('ERROR:',string(type),' is not an accepted type. \n'))
                    return;
            end
            % Allow custom setting of Properties
            
            
            % Initialise a single bin surrounding all given points
            numPts = size(pts,1);
            this.BinBoundaries = [min(pts,[],1) max(pts,[],1)];
            this.Points = pts;
            this.PointBins = ones(numPts,1); %All points start beloning to bin #1
            this.BinDepths = 0; % Bin depth starts at 0
            this.BinParents(1) = 0; % Bin #1 hass no parent
            this.BinCount = 1;
            

            
            % Return on empty or trivial bins
            if numPts<2, return; end
            
            % Start dividing!
            this.preallocateSpace;
            this.divide(1);
            this.deallocateSpace;
            
            
        end
        
        % MATLAB performs better if arrays that grow are initialised,
        % rather than grown during a loop. These two functions do just that
        % before and after the identification of new beens.
        function preallocateSpace(this)
            numPts = size(this.Points,1);
            numBins = numPts;
            if isfinite(this.Properties.binCapacity)
                numBins = ceil(2*numPts/this.Properties.binCapacity);
            end
            this.BinDepths(numBins) = 0; % All bins start att 0
            this.BinParents(numBins) = 0;   % all binns belong to #0
            this.BinBoundaries(numBins,1) = 0;  % all bins have boundary 0
        end
        function deallocateSpace(this)
            this.BinDepths(this.BinCount+1:end) = [];
            this.BinParents(this.BinCount+1:end) = [];
            this.BinBoundaries(this.BinCount+1:end,:) = [];
        end
        
        function divide(this, startingBins)
            % Loop over each bin we will consider for division
            for i = 1:length(startingBins)
                binNo = startingBins(i);
                % Prevent dividing beyond the maximum depth
                if this.BinDepths(binNo)+1 >= this.Properties.maxDepth
                    continue;
                end
                
                % Prevent dividing beyond a minimum size                
                thisBounds = this.BinBoundaries(binNo,:);
                binEdgeSize = diff(thisBounds([1:3;4:6]));
                minEdgeSize = min(binEdgeSize);
                maxEdgeSize = max(binEdgeSize);
                if minEdgeSize < this.Properties.minSize
                    continue;
                end
                
                % There are two conditions under which we should divide
                % this bin. 1: It's bigger than maxSize. 2: It contains
                % more points than binCapacity.
                oldCount = this.BinCount;
                if nnz(this.PointBins==binNo) > this.Properties.binCapacity
                    this.divideBin(binNo);
                    this.divide(oldCount+1:this.BinCount);
                    continue;
                end
                if maxEdgeSize>this.Properties.maxSize
                    this.divideBin(binNo);
                    this.divide(oldCount+1:this.BinCount);
                    continue;
                end
            end
        end
        
        function divideBin(this,binNo)
            % Gather the new points (a bit more efficient to copy once)
            binPtMask = this.PointBins==binNo;
            thisBinsPoints = this.Points(binPtMask,:);
            
            % Get the old corner points and the new division point
            oldMin = this.BinBoundaries(binNo,1:3);
            oldMax = this.BinBoundaries(binNo,4:6);
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
            newBinInds = this.BinCount+1:this.BinCount+8;
            this.BinBoundaries(newBinInds,:) = newBounds;
            this.BinDepths(newBinInds) = this.BinDepths(binNo)+1;
            this.BinParents(newBinInds) = binNo;
            this.PointBins(binPtMask) = newBinInds(binAssignment);
            this.BinCount = this.BinCount + 8;
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
                queryDepth = max(this.BinDepths);
            end
            
            numPts = size(newPts,1);
            newPts = permute(newPts,[3 2 1]);
            binNos = ones(numPts,1)*-1;
                        
            binChildren = arrayfun(@(i)find(this.BinParents==i),1:this.BinCount,'Un',0)';
            binIsLeaf = cellfun(@isempty, binChildren);
            ptQuery_recurse(1:numPts, this.BinParents==0, 0)
            
            function ptQuery_recurse(newIndsToCheck_, binsToCheck, depth)
                % Build a list of all points that fall within one of the
                % bins to be checked, and the list of which point falls in
                % which bin.
                boundsToCheck = this.BinBoundaries(binsToCheck,:);
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
                        ptQuery_recurse(newIndsToCheck_(thisPtMask), binChildren{unqBinNos(i)}, depth+1)
                    end
                end
                
            end
        end
        function plot_bin(this,bin_index)
            hold on;
            binMinMax = this.BinBoundaries(bin_index,:);
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
            h = zeros(this.BinCount,1);
            for i = 1:this.BinCount
                binMinMax = this.BinBoundaries(i,:);
                pts = cat(1, binMinMax([...
                    1 2 3; 4 2 3; 4 5 3; 1 5 3; 1 2 3;...
                    1 2 6; 4 2 6; 4 5 6; 1 5 6; 1 2 6; 1 2 3]),...
                    nan(1,3), binMinMax([4 2 3; 4 2 6]),...
                    nan(1,3), binMinMax([4 5 3; 4 5 6]),...
                    nan(1,3), binMinMax([1 5 3; 1 5 6]));
                h(i) = plot3(pts(:,1),pts(:,2),pts(:,3),varargin{:});
            end
        end
        function h = plot3(this,varargin)
            % OcTree.plot plots bin bounding boxes of an OcTree
            %
            % See also OcTree.plot
            h = this.plot(varargin{:});
        end
        
        
        function mkContent(this)
            % Creates a matrix that displays the points contained in each
            % bin. Only the lowest level bins has a content. 
            tab = tabulate(this.PointBins);
            maxPoints = max(tab(:,2));
            tempCont = zeros(maxPoints,this.BinCount);
            for i = 1:size(this.Points,1)
                check = true;
                n = 1;
                while check                   
                    if tempCont(n,this.PointBins(i)) == 0
                        tempCont(n,this.PointBins(i)) = i;
                        check = false;
                    else
                        n = n + 1;
                    end
                end
            end
            this.BinContent = tempCont;      
        end
        
        function mkChilds(this)
            % creates a matrix indicating what if any binns are contained
            % in the current bin.
            tempChild = zeros(8,this.BinCount);
            for i = 2:this.BinCount
                binNo = this.BinParents(i);
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
            this.BinChilds = tempChild;
        end
        
        function mkBinTri(this)
            % Assigns the triangles to each Bin.
            tempTri = zeros(30,this.BinCount);
            step = 1/this.BinCount;
            process = 0;
            disp('Triangle partition process: ');
                

                % Trying to make the triangle binning paralell
%             for i = 1:this.BinCount
%                 triNums = zeros(size(this.faces,1),1);
%                 triNums = arrayfun(@BoxTriOverlap2, this.BinBoundaries,...
%                     this.Points(this.faces(:,1),:),...
%                     this.Points(this.faces(:,2),:),...
%                     this.Points(this.faces(:,3),:))
%             end


            for j = 1:this.BinCount
                n = 1;
                if this.BinChilds(1,j) == 0
                for i = 1:size(this.faces,1)
                    testTri = [this.Points(this.faces(i,1),:); this.Points(this.faces(i,2),:); this.Points(this.faces(i,3),:)];
                    if BoxTriOverlap(this.BinBoundaries(j,:), testTri)
                        tempTri(n,j) = i;
                        n = n + 1;
                    end
                end
                end
                
                process = process+step;
                disp(process);
            end
            this.BinTri = tempTri;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function ray = OcTreeTravel(this, ray)
% OcTree traversal algorithm. The function determines the traveled path 
% trhough the octree and when to intersect a surface.
%
%     Created by Mikael Ljungholm
%
%     v 1.2           - 2016-04-29. Next bin selection by logical indexes
%                     Current version does not care about "normal"
%                     direction of the plane (back and front). 
%

% IMPORTANT: Must fix the posibility of a ray intersecting a triangle in
% another bin than the same bin the triangle was found. Unless traversal is
% restarted from bin 1 after fefraction. 

% ERROR: Current version does not work if the ray starts within the
% volume.


%------------- Init Values -----------------------------------%
v0 = ray.v0;
v = ray.v;
n1 = ray.n;
n2 = this.refN;
binNo = 1;
% binDepth = 0;
rayEnd = 'Not Ended';
listOfLowerBins = zeros(8,1);
% octreeIndex [z,x,y]
ocIndex = [2 6; 1 5];
ocIndex(:,:,2) = [4 8; 3 7];
currentOcIndex = 0;
orientationMatrix = [0 0 0 5 3 2; 0 0 1 6 4 0; 0 1 0 7 0 4;...
    0 2 3 8 0 0; 1 0 0 0 7 6; 2 0 5 0 8 0; 3 5 0 0 0 8;...
    4 6 7 0 0 0];


%------------- Check to hit boundary volume ---------------------%
[tmin, flag] = rayBoxGPU(this.BinBoundaries(binNo,:),v0(end,:),v);
if flag == 0 % Return if missed volume
    ray.status = 'Ray miss';
    return;
end
bIP = v0 + v*tmin; % Bin intersection point
binDepth = 0;

% testvar = 1;

%------------- Main run loop ---------------------%
missed = 0;
jumpStep = 0;
run = 1;
while run
    if missed == 1 || jumpStep == 1
        n = 3;
        missed = 0;
        jumpStep = 0;
    else
        n = NextMove(); % Find what to do next
    end
    switch n
        case 1
            binNo = FindLowerBin();
            binDepth = binDepth  + 1;
%             testvar = 1;
        case 2
            [iP, intersection, vNew] = FindIntersect();
            if intersection == 1
                v = vNew;
                v0 = [v0 ; iP];
                bIP = iP;
                ray.status = 'Ray intersected';
                run = 0;
            else
                missed = 1;
            end
            % Do refraction.
%             testvar = 1;
        case 3
            [bIP, currentOcIndex] = TraverseBin();
            testvar = 1;
            if currentOcIndex == 0
                jumpStep = 1;
                binDepth = binDepth  - 1;
                binNo = this.BinParents(binNo);  
%                 testvar = 1;
                if binNo == 1
                    ray.status = 'Ray lost in space';
                    % Ray exit volume
                    run = 0;
                else
                    tempInds = this.BinChilds(:,this.BinParents(binNo));
                    currentOcIndex = find(tempInds == binNo);
%                     testvar = 1;
                end
            else
                binNo = this.BinChilds(currentOcIndex,this.BinParents(binNo));
            end
%             testvar = 1;
    end
end
ray.v0 = v0;
ray.v = v;

%------------ Sub functions -----------------------%
    function flag = NextMove()
        % Determining the next step. Fist checks if there are lower bins in
        % the current one. Next checks if the bin is empty. If both are
        % empty then move to edge of bin. 
%         testvar = 1;
        if this.BinChilds(1,binNo) ~= 0 
            flag = 1; 
        elseif this.BinTri(1,binNo) ~= 0
            flag = 2; 
        else
            flag = 3;
        end
%         testvar = 1;
    end


    function nextBin = FindLowerBin()
        % Uses current intersection point with the local logical indices to
        % determine the next lower bin.
        xm = 1;
        ym = 1;
        zm = 2;
        xh = (this.BinBoundaries(binNo,4)-this.BinBoundaries(binNo,1))/2;
        yh = (this.BinBoundaries(binNo,5)-this.BinBoundaries(binNo,2))/2;
        zh = (this.BinBoundaries(binNo,6)-this.BinBoundaries(binNo,3))/2;
%         testvar = 1;
        if (bIP(1)-xh) > this.BinBoundaries(binNo,1)
            xm = 2;
        end
        if (bIP(2)-yh) > this.BinBoundaries(binNo,2)
            ym = 2;
        end
        if (bIP(3)-zh) > this.BinBoundaries(binNo,3)
            zm = 1;
        end
        % octreeIndex [z,x,y]
        currentOcIndex = ocIndex(zm,xm,ym);
        nextBin = this.BinChilds(currentOcIndex,binNo);
%         testvar = 1;
    end

    function [iP, intersection, vNew] = FindIntersect()
        % Creates a list of triangles in the current bin and checks for any
        % intersections. Calculate refraction if intersection is found.
        triList = nonzeros(this.BinTri(:,binNo));
        closest = inf;
        closestIp = 0;
        triNum = 0;
        intersection = 0;
        vNew = v;
        for i = 1:length(triList)
            Tv0 = this.Points(this.faces(triList(i),1),:);
            Tv1 = this.Points(this.faces(triList(i),2),:);
            Tv2 = this.Points(this.faces(triList(i),3),:);
            u = Tv1-Tv0;
            w = Tv2-Tv0;
            N = cross(u,w);
            N = N./sqrt(N(1)^2 + N(2)^2 +N(3)^2);
            [iP, intersect] = RayTriIntersect(Tv0,Tv1,Tv2,N, v, v0(end,:));
            if intersect == 1
                dist = norm(iP-v0(end,:));
                if  dist < closest
%                     testvar = 1;
                    closest = dist;
                    closestIp = iP;
                    triNum = i;
                    intersection = 1;
                    triN = N;
                end
            end
        end
        iP = closestIp;
        if intersection == 1
            [vNew, reflected] = Snell(v, triN, n1, n2);
        end
    end

    function [nextP, nextInd] = TraverseBin()
        % Moves from current point in bin to bin edge. Then cheks if there 
        % are any more bins in the same level, otherwise moves up one 
        % level. if this is max level, then exits the volume.
        
%         parenBin = OT.BinParents(binNo);
        % currentOcIndex. use it to determine next bin or parent borders.
        
        tx = 0; ty = 0; tz = 0;
        d = zeros(1,3);
        if v(1) >= 0
            tx = (this.BinBoundaries(binNo,4)-bIP(1))/v(1);
            d(1) = 3;
        else
            tx = (this.BinBoundaries(binNo,1)-bIP(1))/v(1);
        end
        if v(2) >= 0
            ty = (this.BinBoundaries(binNo,5)-bIP(2))/v(2);
            d(2) = 3;
        else
            ty = (this.BinBoundaries(binNo,2)-bIP(2))/v(2);
        end
        if v(3) >= 0
            tz = (this.BinBoundaries(binNo,6)-bIP(3))/v(3);
            d(3) = 3;
        else
            tz = (this.BinBoundaries(binNo,3)-bIP(3))/v(3);
        end
        [t, I] = min([tx,ty,tz]);
        nextInd = orientationMatrix(currentOcIndex,I+d(I));
%         [nexOcInd, oob] = Orientation();
        nextP = bIP + v*t;
    end

%     function [nexOcInd, oob] = Orientation()
%         % Uses the orientation matrix to determine what the next OcTree
%         % index and if the ray exits the parent bin. (oob = out of bounds)
%     end
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
%        cols = lines(OT.BinCount);
%        doplot3 = @(p,varargin)plot3(p(:,1),p(:,2),p(:,3),varargin{:});
%        for i = 1:OT.BinCount
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
%     BinCount          - Total number of bins created.
%     BinBoundaries     - BinCount-by-6 [MIN MAX] coordinates of bin edges.
%     BinDepths         - The # of subdivisions to reach each bin.
%     BinParents        - Indices of the bin that each bin belongs to.
%     Properties        - Name/Val pairs used for creation (see help above)
%     BinContent        - List of vetricis contained in each bin.
%     BinTri            - List of trianlge faces contained in each bin.
%     BinChilds         - List of child bins for each parent bin.
%     BinPrim           - List of primitives for each parent bin. 


% Possible additional info "BinContent" - idices of the points in each bin
