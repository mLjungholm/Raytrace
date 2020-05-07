% This is the class for holding information on the retina. Including
% receptors, shapes, positions, absorption coeficients etc. It is used with
% the receptor trace function to calculate absorption in the retina.

%                           !Note!
%           For now this only works with cone shaped receptors

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
        aborption_coeff     % Absorption coeficcent for the receptors
        absorbed_val        % Total absorbed value for each receptor
        x_grid              % X coordinates for the voxel centers
        y_grid              % Y coordinates for the voxel centers
        z_grid              % Z coordinates for the voxel centers
        base_main           % Base diameter of the main cone
        base_sec            % Base diameter of the seconday cone
        alpha               % Apperature angle for the main cone
        beta                % Apperature angle for the secondary cone
    end
    
    methods
        function this = Receptor_space(varargin)
            % Header line for Receptor space
            switch varargin{1}
                case ndim(varargin{1}) == 2 && size(varargin{1},2) == 3
                    % origin receptor points
                    if size(varargin{2}) ~= size(varargin{1})
                        error('Second input(end positions) must be the same size as the first input(start position)')
                    end
                    this.base_pos = varargin{1};
                    this.end_pos = varargin{2};
                    this.receptor_nums = size(this.base_pos,1);                    
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
                    error('Unsuported input type')
            end
        end
            
        % Allocation function for creating the rreceptor space with a
        % defined grid fineness/resolution.
        function allocate_space(this, volume_resolution)
            this.volume_resolution = volume_resolution;
            this.step_size = (this.box_max-this.box_min)./(volume_resolution-1);            
            this.x_grid = minP(1):stepS(1):maxP(1);
            this.y_grid = minP(2):stepS(2):maxP(2);
            this.z_grid = minP(3):stepS(3):maxP(3);            
            % The cell array is created with the y coordinate fist so that the matrix
            % indexing mathces the physical carteesian coordinate structure.
            this.receptor_grid = cell(steps(2),steps(1),steps(3)); 
        end        
        
        % Fitt all the receptors into the allocated receptor grid.
        function fit_receptors2bins(this)
            if isnan(this.receptor_grid)
                error('Receptor grid has not been allocated')
            end
            cone_dir = single(zeros(this.receptor_nums,3));
            aP = single(zeros(this.receptor_nums,3));
            h = single(zeros(this.receptor_nums,1));
            % Calculate the length, direction and angles of the cones
            for j = 1:size(this.receptor_nums,1)                
                v = this.base_pos(j,:)-this.end_pos(j,:);
                h(j) = norm(v);
                cone_dir(j,:) = v./h(j);
                aP(j,:) = (this.end_pos(j,:)-2*v);
                this.alpha(j) = atand(4/(3*h(j)));
                this.beta(j) = atand(0.7/(3*h(j)));
            end
            
            hWaitBar = waitbar(0, 'Fitting retina Volume', 'CreateCancelBtn', ...
                @(src, event) setappdata(gcbf(), 'Cancelled', true));
            setappdata(hWaitBar, 'Cancelled', false);
            numComplete = 0;
            % Initiate fitting
            pcount = 1:this.receptor_nums;
            pcount = pcount';
            
            % Parralell fitt each cone ot all points
            try
                arrayfun(@point_cone_fit, pcount);
            catch
                delete(hWaitBar);
                error('unexpected error in fitting')
            end
            delete(hWaitBar);
            
            x_nums = this.volume_resolution(1);
            y_nums = this.volume_resolution(2);
            z_nums = this.volume_resolution(3);
            function match = point_cone_fit(p_ind)
                %                     match = single(zeros(1,100)); % Guessing no more than hundred cones for a point
                n = 1;
                for z_ind = 1:z_nums
                    for y_ind = 1:y_nums
                        for x_ind = 1:x_nums
                            l = sqrt((this.x_grid(x_ind)-aP(p_ind,1))^2 + (this.y_grid(y_ind)-aP(p_ind,2))^2 + (this.z_grid(z_ind)-aP(p_ind,3))^2);
                            if  l > h(p_ind)*2 && l < h(p_ind)*3
                                
                                
                                
                                %%% Do this next %%%
                                pvec = (points(p_ind,:)-cone_ap(x_ind,:))./norm(points(p_ind,:)-cone_ap(x_ind,:));
                                ang = acosd(pvec(1)*cone_dir(x_ind,1) + pvec(2)*cone_dir(x_ind,2) + pvec(3)*cone_dir(x_ind,3));
                        %                 if  ang <= beta(i) && ang >= alpha(i)
                        if  ang <= alpha(x_ind)
                            match(n) = x_ind;
                            n = n + 1;
                        end
                    end
                end
                numComplete = numComplete + 1;
                fractionComplete = numComplete / this.receptor_nums;
                waitbar(fractionComplete, hWaitBar);
            end
        end
    end
end