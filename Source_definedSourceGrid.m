% The source contains arrays for all ray-values

% The assumes max path =  10. If you know there will be more refractions
% than 10 or if any volume is inhomogenious then this value needs to be
% increased

% Found an issue where an extra surface like the cornia can block the rays
% between two ther surfaces (like retina and lens)

classdef Source_definedSourceGrid < handle
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
        refract_index;    % current refractive index of each ray.
        status;
        ray_alive;
    end
    
    methods
        function this = Source_definedSourceGrid(sourceGrid, refract_index, dir_Point)
            this.num_rays = size(sourceGrid,1);
            this.v0 = sourceGrid;
            [this.path_x, this.path_y, this.path_z] = deal(zeros(this.num_rays,10));
            this.status = strings(this.num_rays,1);
            this.ray_alive = ones(this.num_rays,1);
            this.steps = ones(this.num_rays,1);
            for i = 1:this.num_rays
                this.v(i,:) = (dir_Point - this.v0(i,:))./norm(dir_Point - this.v0(i,:));
                this.refract_index(i) = refract_index;
                this.path_x(i) = this.v0(i,1);
                this.path_y(i) = this.v0(i,2);
                this.path_z(i) = this.v0(i,3);
            end
        end
        
        function plot_grid(this,figureNr)
            figure(figureNr)
            hold on
            axis equal
            for i = 1:this.num_rays
                plot3(this.v0(i,1,1),this.v0(i,2,1),this.v0(i,3,1),'.r')
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
        
        function end_points = export_endpoints(this)
            end_points = zeros(this.num_rays,3);
            for i = 1:this.num_rays
                if ~this.steps(i) == 0
                    end_points(i,:) = [this.path_x(i,this.steps(i)),this.path_y(i,this.steps(i)),this.path_z(i,this.steps(i))];
                end
            end
        end
        
    end
end