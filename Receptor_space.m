% This is the class for holding information on the retina. Including
% receptors, shapes, positions, absorption coeficients etc. It is used with
% the receptor trace function to calculate absorption in the retina.

%                           !Note!
%           For now this only works with cone shaped receptors

% Input variables:
%       Alt1: (Origin_points/base_points, end_points, 
%             [base_diameter, secondary diameter])



classdef Receptor_space < handle
    
    properties
        base_pos    % The base positions of the recepotrs. This is what determines the location of the receptor.
        end_pos     % End position of the receprors.
        box_min     % Min positions of the receptor grid
        box_max     % Max size of the receptor grid
        step_size   % Size of the receptor voxels
        receptor_nums       % Number of receoptors
        volume_resolution   % The number of steps in the receptor grid
        volume_id           % What volume id is acosiated with the retina.
        receptor_grid       % The 3D cell array where each indicies contains a list of the receptors contained in that voxel.
        absorption_coeff     % Absorption coeficcent for the receptors
        absorbed_val        % Total absorbed value for each receptor
        x_grid              % X coordinates for the voxel centers
        y_grid              % Y coordinates for the voxel centers
        z_grid              % Z coordinates for the voxel centers
        base_main           % Base radius of the main cone
        base_sec            % Base radius of the seconday cone
        alpha               % Apperature angle for the main cone
        beta                % Apperature angle for the secondary cone
        unmatched_receptors
        receptor_points     % Cell array with all points cointained withen each receptor (temporary test variable)
    end    
    
    methods
        function this = Receptor_space(varargin)
            % Header line for Receptor space
            indim = ndims(varargin{1});
            switch indim
                case 2 % && size(varargin{1},2) == 3
                    % origin receptor points
                    if size(varargin{2}) ~= size(varargin{1})
                        error('Second input(end positions) must be the same size as the first input(start position)')
                    end
                    this.base_pos = varargin{1};
                    this.end_pos = varargin{2};
                    this.receptor_nums = size(this.base_pos,1);
                    this.absorbed_val = zeros(this.receptor_nums,1);
                    this.box_min = min([this.end_pos;this.base_pos],[],1);
                    this. box_max = max([this.end_pos;this.base_pos],[],1);
                    if ndims(varargin{3}) > 1
                        this.base_main = varargin{3}(1);
                        this.base_sec = varargin{3}(2);
                    else
                        this.base_main = varargin{3};
                        this.base_sec = 0;
                    end
                    this.alpha = single(zeros(this.receptor_nums,1));
                    this.beta = single(zeros(this.receptor_nums,1));
                otherwise
                    errmsg = strcat('Unsuported input type : ', num2str(indim));
                    error(errmsg)
            end
        end
            
        % Allocation function for creating the rreceptor space with a
        % defined grid fineness/resolution.
        function allocate_space(this, volume_resolution)
            this.volume_resolution = volume_resolution;
            this.step_size = (this.box_max-this.box_min)./(volume_resolution-1);            
            this.x_grid = this.box_min(1):this.step_size(1):this.box_max(1);
            this.y_grid = this.box_min(2):this.step_size(2):this.box_max(2);
            this.z_grid = this.box_min(3):this.step_size(3):this.box_max(3);            
            % The cell array is created with the y coordinate fist so that the matrix
            % indexing mathces the physical carteesian coordinate structure.
            this.receptor_grid = cell(volume_resolution(2),volume_resolution(1),volume_resolution(3)); 
        end        
        
        % Fitt all the receptors into the allocated receptor grid.
        function fit_receptors2bins(this)
            if ~isa(this.receptor_grid,'cell')
                error('Receptor grid has not been allocated')
            end
            this.receptor_points = cell(this.receptor_nums,1); % temporary variable
            cone_dir = single(zeros(this.receptor_nums,3));
            aP = single(zeros(this.receptor_nums,3));
            h = single(zeros(this.receptor_nums,1));
            % Calculate the length, direction and angles of the cones
            for j = 1:this.receptor_nums
                v = this.base_pos(j,:)-this.end_pos(j,:);
                h(j) = norm(v);
                cone_dir(j,:) = v./h(j);
                aP(j,:) = (this.end_pos(j,:)-v);
                this.alpha(j) = atand(this.base_main/(2*h(j)));
                this.beta(j) = atand(this.base_sec/(2*h(j)));
            end
            
            hWaitBar = waitbar(0, 'Fitting retina Volume', 'CreateCancelBtn', ...
                @(src, event) setappdata(gcbf(), 'Cancelled', true));
            setappdata(hWaitBar, 'Cancelled', false);
            numComplete = 0;
            % Initiate fitting
            pcount = 1:this.receptor_nums;
            pcount = pcount';
            % Fitting function for the points and cones
            this.unmatched_receptors = single(zeros(this.receptor_nums,1));
            
            % Parralell fitt each cone ot all points
            try
                arrayfun(@point_cone_fit, pcount);
            catch
                delete(hWaitBar);
                error('unexpected error in fitting')
            end
            delete(hWaitBar);
            this.unmatched_receptors = this.unmatched_receptors(any(this.unmatched_receptors,2));
            
            function point_cone_fit(p_ind)
                found_any = 0;
                for z_ind = 1:this.volume_resolution(3)
                    for y_ind = 1:this.volume_resolution(2)
                        for x_ind = 1:this.volume_resolution(1)
                            % calculate the coodinates for the index
                            % point
                            point = [this.x_grid(x_ind) this.y_grid(y_ind) this.z_grid(z_ind)];
                            l = sqrt((point(1)-aP(p_ind,1))^2 + (point(2)-aP(p_ind,2))^2 + (point(3)-aP(p_ind,3))^2);
                            if  l > h(p_ind) && l < h(p_ind)*2
                                % Vector from cone apex to grid point
                                pvec = (point-aP(p_ind,:))./norm(point-aP(p_ind,:));
                                ang = acosd(pvec(1)*cone_dir(p_ind,1) + pvec(2)*cone_dir(p_ind,2) + pvec(3)*cone_dir(p_ind,3));
                                % Check if the intervector angle is inside
                                % the cone apex angle.
                                %  if  ang <= beta(i) && ang >= alpha(i)
                                if  ang <= this.alpha(p_ind)
                                    % This changes the size for the vectors
                                    % in the cell array. It might be slow.
                                    this.receptor_points{p_ind} = [this.receptor_points{p_ind}; point]; 
                                    this.receptor_grid{y_ind,x_ind,z_ind} = [this.receptor_grid{y_ind,x_ind,z_ind} p_ind];
%                                     this.receptor_grid{end-y_ind+1,x_ind,z_ind} = [this.receptor_grid{end-y_ind+1,x_ind,z_ind} p_ind];
                                    found_any = 1;
                                end
                            end
                        end
                    end
                end
                if ~found_any
                    this.unmatched_receptors(p_ind) = p_ind;
                end
                % Uppdate the waitbar
                numComplete = numComplete + 1;
                fractionComplete = numComplete / this.receptor_nums;
                waitbar(fractionComplete, hWaitBar);
            end
        end
        
        % Creates a series of control variables
        function [all_points, not_empty, empty_points] = check_fit_errors(this)
            minP = min([this.end_pos;this.base_pos],[],1);
            maxP = max([this.end_pos;this.base_pos],[],1);
            [x,y,z] = meshgrid(minP(1):this.step_size(1):maxP(1),minP(2):this.step_size(2):maxP(2),minP(3):this.step_size(3):maxP(3));
            all_points = [reshape(x,[],1),reshape(y,[],1),reshape(z,[],1)];
            all_points = single(all_points);
            
            not_empty = single(zeros(size(all_points,1),3));
            i = 1;
            for z_ind = 1:this.volume_resolution(3)
                for y_ind = 1:this.volume_resolution(2)
                    for x_ind = 1:this.volume_resolution(1)
                        if size(this.receptor_grid{y_ind,x_ind,z_ind},1) ~= 0
%                         if size(this.receptor_grid{end-y_ind+1,x_ind,z_ind},1) ~= 0
                            not_empty(i,:) = [this.x_grid(x_ind),this.y_grid(y_ind),this.z_grid(z_ind)];
                            i = i + 1;
                        end
                    end
                end
            end
            mapof_empty = ~ismember(not_empty,[0 0 0],'rows');
            not_empty = not_empty(mapof_empty,:);
            empty_map = ~ismember(all_points,not_empty,'rows');
            empty_points = all_points(empty_map,:);
        end
        
        function plot_cone(this,conenr,figurenr,color)
            figure(figurenr)
            hold on
            if size(conenr,2) > 1
                for i = 1:size(conenr,2)
                    plot_cone(this.base_pos(conenr(i),:),this.end_pos(conenr(i),:),this.base_main, figurenr,color)
                end
            else
                plot_cone(this.base_pos(conenr,:),this.end_pos(conenr,:),this.base_main, figurenr,color)
            end
            axis equal
        end
        
        function plot_conePoints(this,conenr,figurenr)
            % Guess max points
            found_points = zeros(this.volume_resolution(1)*this.volume_resolution(2)*this.volume_resolution(3)*0.1,3);
            i = 1;
            for zi = 1:this.volume_resolution(3)
                for yi = 1:this.volume_resolution(2)
                    for xi = 1:this.volume_resolution(1)
                        if ~size(this.receptor_grid{yi,xi,zi},1) == 0
%                         if ~size(this.receptor_grid{end-yi+1,xi,zi},1) == 0
                            if any(this.receptor_grid{yi,xi,zi}(this.receptor_grid{yi,xi,zi} == conenr))
%                             if any(this.receptor_grid{end-yi+1,xi,zi}(this.receptor_grid{end-yi+1,xi,zi} == conenr))
                                found_points(i,:) = [this.x_grid(xi),this.y_grid(yi),this.z_grid(zi)];
                                i = i + 1;
                            end
                        end
                    end
                end
            end
            map = ~ismember(found_points,[0 0 0],'rows');
            found_points = found_points(map,:);
            figure(figurenr)
            scatter3(found_points(:,1),found_points(:,2),found_points(:,3),'.')
        end
        
        function plot_conePoints2(this,conenr,figurenr)
            figure(figurenr)
            scatter3(this.receptor_points{conenr}(:,1),this.receptor_points{conenr}(:,2),this.receptor_points{conenr}(:,3),'.')
        end
        
        function plot(this, varargin)
        % Inputs: plot(receptors, type, figureNr)
        %         1: List of receptor indices or 'all'
        %         2: type of plot. 'line', 'base' or 'end'
        %         3: FigureNr
            switch varargin{2}
                case 'base'
                    if isequal(varargin{1}, 'all')
                        figure(varargin{3})
                        axis equal
                        scatter3(this.base_pos(:,1),this.base_pos(:,2),this.base_pos(:,3),'r.')
                    elseif ismatrix(varargin{1})                        
                        figure(varargin{3})
                        axis equal
                        scatter3(this.base_pos(varargin{1}',1),this.base_pos(varargin{1}',2),this.base_pos(varargin{1}',3),'r.')
                    else
                        error('Receptors must be list of indicies or "all"')
                    end                    
                case 'end'
                    if isequal(varargin{1}, 'all')
                        figure(varargin{3})
                        axis equal
                        scatter3(this.end_pos(:,1),this.end_pos(:,2),this.end_pos(:,3),'r.')
                    elseif ismatrix(varargin{1})                        
                        figure(varargin{3})
                        axis equal
                        scatter3(this.end_pos(varargin{1}',1),this.end_pos(varargin{1}',2),this.end_pos(varargin{1}',3),'r.')
                    else
                        error('Receptors must be list of indicies or "all"')
                    end
                case 'line'
                    if isequal(varargin{1}, 'all')
                        figure(varargin{3})
                        axis equal
                        hold on
                        for i = 1:this.receptor_nums
                            l = [this.base_pos(i,:); this.end_pos(i,:)];
                            plot3(l(:,1),l(:,2),l(:,3),'r')
                        end
                    elseif ismatrix(varargin{1})                        
                        figure(varargin{3})
                        axis equal
                        hold on
                        for i = 1:size(varargin{1},2)
                            l = [this.base_pos(varargin{1}(i),:); this.end_pos(varargin{1}(i),:)];
                            plot3(l(:,1),l(:,2),l(:,3),'r')
                        end
                    else
                        error('Receptors must be list of indicies or "all"')
                    end
                otherwise
                    error('Unsupported plot type')
            end
        end
            
            
            
%             if isequal(varargin{1}, 'all')
%                 % plott all
%             elseif ismatrix(varargin{1})
%                 % plot list
%             else
%                 error('Unsupported inpput type')
%             end
%             
%             figure(figureNr)
%             hold on
%             axis equal
%             scatter3(this.base_pos(:,1),this.base_pos(:,2),this.base_pos(:,3),'r.')
        
    end
end