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
        r_index;    % current refractive index of each ray.
        status;
        ray_alive;
    end
    
    methods
        function this = Source(dir, origin, hsize, raysPerAxis, r_index)
            this.source_dir = dir./norm(dir);
            this.origin = origin;
            if raysPerAxis == 1
                this.num_rays = 1;
                this.v = this.source_dir;
                this.v0 = origin;
                this.r_index = r_index;
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
%                 this.r_index = ones(this.num_rays,1);
                this.status = strings(this.num_rays,1);
                this.ray_alive = ones(this.num_rays,1);
                this.steps = ones(this.num_rays,1);
                for i = 1:this.num_rays
                    this.v(i,:) = dir;
                    this.v0(i,:) = tP(i,:);
                    this.r_index(i) = r_index;
                    this.path_x(i) = tP(i,1);
                    this.path_y(i) = tP(i,2);
                    this.path_z(i) = tP(i,3);
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
        
%         function plot(this, figureNr, lengthOfStray)
%             figure(figureNr)
%             hold on
%             axis equal
%             for i = 1:this.num_rays
%                 if this.num_rays == 1
%                     tempLine = [this.v0(1,:); this.v0(1,:) + this.v*lengthOfStray];
%                     plot3(tempLine(:,1),tempLine(:,2),tempLine(:,3),'--r')
%                 else
%                     plot3(this.v0(:,1),this.v0(:,2),this.v0(:,3),'r')
% %                     tempLine = [this.v0(end,:);this.v0(end,:) + this.v*lengthOfStray];
% %                     plot3(tempLine(:,1),tempLine(:,2),tempLine(:,3),'--r')
%                 end
%             end
%         end
        
    end
end