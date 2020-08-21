% The source contains arrays for all ray-values

% The assumes max path =  10. If you know there will be more refractions
% than 10 or if any volume is inhomogenious then this value needs to be
% increased

classdef Source < handle
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
%         abs_vol     % Indicies of the absorption volumes (if any).
        step_absVol    % indicies for the absorbation volume at each step index. ex [0 0 1 1] for volume 1 at step 3 and  4
        refract_index;    % current refractive index of each ray.
        status;
        ray_alive;
        source_id;
    end
    
    methods
        function this = Source(dir, origin, hsize, raysPerAxis, refract_index)
            for i = 1:3
                if dir(i) == 0
                    dir(i) = 0;
                end
            end
            this.source_dir = dir./norm(dir);
            this.origin = origin;
            if raysPerAxis == 1
                this.num_rays = 1;
                this.v = this.source_dir;
                this.v0 = origin;
                this.refract_index = refract_index;
            else
                dir = dir./norm(dir);
                if dir(2) == 0 && dir(3) == 0
                    u = [-dir(3) 0 dir(1)];
                else
                    u = [0 -dir(3) dir(2)];
                end
                w = cross(u,dir);
                u1 = u/norm(u);
                w1 = w/norm(w);
                du = u1.*hsize.*2./(raysPerAxis-1);
                dw = w1.*hsize.*2./(raysPerAxis-1);
                tP = zeros(raysPerAxis * raysPerAxis,3);
                n = 1;
                for i = 1:raysPerAxis
                    for j = 1:raysPerAxis
                        tP(n,:) = origin -hsize.*u1 + du*(i-1) -hsize.*w1 + dw*(j-1);
                        if norm(tP(n,:)-origin) > hsize
                            tP(n,:) = nan;
                        end
                        n = n + 1;
                    end
                end
                tP = tP(:);
                tP = tP(~isnan(tP));
                tP = reshape(tP,[],3);
                this.num_rays = size(tP,1);
                [this.v, this.v0] = deal(zeros(this.num_rays,3));
                [this.path_x, this.path_y, this.path_z] = deal(zeros(this.num_rays,10));
%                 this.refract_index = ones(this.num_rays,1);
                this.status = strings(this.num_rays,1);
                this.ray_alive = ones(this.num_rays,1);
                this.steps = ones(this.num_rays,1);
                for i = 1:this.num_rays
                    this.v(i,:) = dir;
                    this.v0(i,:) = tP(i,:);
                    this.refract_index(i) = refract_index;
                    this.path_x(i) = tP(i,1);
                    this.path_y(i) = tP(i,2);
                    this.path_z(i) = tP(i,3);
                    this.step_absVol = zeros(this.num_rays,5);
                end
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
        function plot_single(this,raynr,stray_length)
            figure(1)
            hold on
            temp_line(1,:) = [this.path_x(raynr,this.steps(raynr))',this.path_y(raynr,this.steps(raynr))',this.path_z(raynr,this.steps(raynr))'];
            temp_line(2,:) = temp_line(1,:) + this.v(raynr,:).*stray_length;
            plot3(temp_line(:,1),temp_line(:,2),temp_line(:,3),'r--')
        end
        function unblock(this)
            for ray_ind = 1:this.num_rays
                if isequal(this.status(ray_ind),'Ray ended')
                    this.status(ray_ind) = 'refract';
                    this.ray_alive(ray_ind) = 1;
                end
            end
        end
        function plot_spec(this, figureNr, color, linewidth)
            figure(figureNr)
            hold on
            axis equal
            for i = 1:this.num_rays
                temp_line = [this.path_x(i,1:this.steps(i))',this.path_y(i,1:this.steps(i))',this.path_z(i,1:this.steps(i))'];
                plot3(temp_line(:,1),temp_line(:,2),temp_line(:,3),color,'linewidth',linewidth)
            end
        end
        
    end
end