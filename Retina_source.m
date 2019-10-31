% The source contains arrays for all ray-values

% The assumes max path =  10. If you know there will be more refractions
% than 10 or if any volume is inhomogenious then this value needs to be
% increased

classdef Retina_source < handle
    properties
        
        origin;     % Origin of source
        source_dir; % Direction of source
        num_rays;    % Number of rays
        wavelength; % Wavelength of the rays
        v;          % (nx3) direction of n rays [x,y,z].
        path_x;
        path_y;
        path_z;
        v0;         % (nx3xm) path of ray. m = number of stepps. sett to a hig number to avoid reconstructing matrix
        steps;      % (nxm) number of steps taken. starts att 1
        absorption  % (nx3xm) m = number of steps.
        r_index;    % current refractive index of each ray.
        status;
        ray_alive;
    end
    
    methods
        function this = Retina_source(grid,target)
            this.v = (target-grid)./norm(target-grid);
            this.v0 = grid;
            [this.path_x, this.path_y, this.path_z] = deal(zeros(this.num_rays,10));
            this.status = strings(this.num_rays,1);
            this.ray_alive = ones(this.num_rays,1);
            this.steps = ones(this.num_rays,1);
            this.r_index = ones(size(grid,1),1);
            this.path_x = this.v0(:,1);
            this.path_y = this.v0(:,2);
            this.path_z = this.v0(:,3);
            this.num_rays = size(grid,1);
            this.ray_alive = ones(size(grid,1),1);
            this.steps = ones(this.num_rays,1);
            
        end
        
        function plot_grid(this,figureNr)
            figure(figureNr)
            hold on
            axis equal
            for i = 1:this.num_rays
                plot3(this.v0(i,1),this.v0(i,2),this.v0(i,3),'.r')
            end
        end
        
        function plot(this, figureNr)
            figure(figureNr)
            hold on
            axis equal
            for i = 1:this.num_rays
                temp_line = [this.path_x(i,1:this.steps(i))',this.path_y(i,1:this.steps(i))',this.path_z(i,1:this.steps(i))'];
                plot3(temp_line(:,1),temp_line(:,2),temp_line(:,3),'r')
            end
        end
        function plot_stray(this,figureNr, stray_length)
            figure(figureNr)
            hold on
            axis equal
            for i = 1:this.num_rays
                temp_line = zeros(2,3);
                temp_line(1,:) = [this.path_x(i,this.steps(i))',this.path_y(i,this.steps(i))',this.path_z(i,this.steps(i))'];
                temp_line(2,:) = temp_line(1,:) + this.v(i,:).*stray_length;
                plot3(temp_line(:,1),temp_line(:,2),temp_line(:,3),'r--')
            end
        end
        
    end
end