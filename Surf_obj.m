classdef Surf_obj < handle
    % Surface object
    %
    % Pre OcTree creation object. Used for visualization and measuring
    % berfore OcTree creation. 
    
    properties
        v;
        f;
        n;
        name;
        numIntersections;
        nIndexOutside;
        nIndexInside;
        maxVertDist;
        center;
        centroid;
        v_count;
        f_count;
    end
    
    methods
        function this = Surf_obj()
            [FileName,PathName] = uigetfile('*.stl','Select surfaces to load');
            filename = strcat(PathName,FileName);
            [v, f, n] = Stlread(filename);
            
            % find uniqe vertex points
            vNew = unique(v,'rows');
            faceNums = size(f,1);
            fNew = zeros(faceNums,3);
            
            for i = 1:faceNums
                vVal = [v(f(i,1),:);v(f(i,2),:);v(f(i,3),:)];
                newInd = zeros(1,3);
                newInd(1) = find(ismember(vNew,vVal(1,:),'rows'));
                newInd(2) = find(ismember(vNew,vVal(2,:),'rows'));
                newInd(3) = find(ismember(vNew,vVal(3,:),'rows'));
                fNew(i,:) = newInd;
            end
            this.v = vNew;
            this.f = fNew;
            this.v_count = size(vNew,1);
            this.f_count = faceNums;
            for i = 1:faceNums
                n(i) = n(i)./norm(n(i));
            end
            this.n = n;
            this.name = FileName(1:end-4);
        end
        
        function remakeN(this)
            for i = 1:size(this.f,1)
                u = this.v(this.f(i,2),:)-this.v(this.f(i,1),:);
                w = this.v(this.f(i,3),:)-this.v(this.f(i,1),:);
                nt = [(u(2)*w(3) - u(3)*w(2)) (u(3)*w(1) - u(1)*w(3)) (u(1)*w(2) - u(2)*w(1))];
                nt = nt./sqrt(nt(1)^2 + nt(2)^2 + nt(3)^2);
                this.n(i,:) = nt;
            end
        end
        function plot(this,figureNr,color)
            figure(figureNr)
%             hold on
            axis equal
            trisurf(this.f,this.v(:,1),this.v(:,2),this.v(:,3),'Facecolor',color,'FaceAlpha',0.1,'EdgeAlpha',0.3)
        end
        function plot_centers(this, figureNr)
            figure(figureNr)
            plot3(this.centroid(1),this.centroid(2),this.centroid(3),'ro')
            plot3(this.center(1),this.center(2),this.center(3),'go')
        end
        function plot_normals(this,figureNr,length)
            figure(figureNr)
            hold on
            for i = 1:size(this.n,1)
                op = this.v(this.f(i,1),:);
%                 tline = [op; op + this.n.*length];
                quiver3(op(1),op(2),op(3),this.n(i,1),this.n(i,2),this.n(i,3),length);
%                 quiver3(tline(:,1),tline(:,2),tline(:,3),'r--')
            end
        end
        function plot_triangles(this,tri_list)
            faces = zeros(size(tri_list,2),3);
            for i = 1:size(tri_list,2)
                faces(i,:) = this.f(tri_list(i),:);
            end
            trisurf(faces,this.v(:,1),this.v(:,2),this.v(:,3),'Facecolor','b','FaceAlpha',0.1,'EdgeAlpha',0.3)
        end
        
        function RotSurf(this, xyz, point, angle)
            % Function that rotates the surface around an axis and point in
            % space. 
            if xyz == 'x'
                R = [1 0 0; 0 cosd(angle) -sind(angle); 0 sind(angle) cosd(angle)];          
            elseif xyz == 'y'
                R = [cosd(angle) 0 sind(angle); 0 1 0; -sind(angle) 0 cosd(angle)];
            elseif xyz == 'z'
                R = [cosd(angle) -sind(angle) 0; sind(angle) cosd(angle) 0; 0 0 1];
            else
                disp('Error: not valid axis ');
                return;
            end
            this.v = bsxfun(@minus,this.v,point);
            this.v = this.v*R;
            this.v = bsxfun(@plus,this.v,point);
            this.n = this.n*R;          
        end
        
        function ScaleSurf(this, scaling)
            this.v = this.v.*scaling;
        end
        function move(this, vector)
            this.v = bsxfun(@plus,this.v,vector);
        end
        
        function find_centroid(this)
            x = mean(this.v(:,1));
            y = mean(this.v(:,2));
            z = mean(this.v(:,3));
            this.centroid = [x y z];
        end
        
        function find_center(this)
            x = min(this.v(:,1))+(max(this.v(:,1))-min(this.v(:,1)))/2;
            y = min(this.v(:,2))+(max(this.v(:,2))-min(this.v(:,2)))/2;
            z = min(this.v(:,3))+(max(this.v(:,3))-min(this.v(:,3)))/2;
            this.center = [x y z];
        end
        
            function ray = TraceRays(this, ray)
                if isempty(this.maxVertDist)
                   disp('ERROR: Max vertex distance is not defined');
                   return;
                end
%                 if ray.terminated == 1
%                     return;
%                 end
            vertIndex = vertexReduction(this,ray);
            closest = -1;
            intersectFound = 0;
            closestSurfNum = 0;
            newPoint = [0 0 0];
            
            for i = 1:size(vertIndex,2)
                Tv0 = this.v(this.f(vertIndex(i),1),:);
                Tv1 = this.v(this.f(vertIndex(i),2),:);
                Tv2 = this.v(this.f(vertIndex(i),3),:);

                [I, intersects] = RayTriIntersect(Tv0,Tv1,Tv2, this.n(vertIndex(i),:),ray.v, ray.v0);
        
                if intersects == 1
                    if closest < 0 || norm(I-ray.v0) < closest
                        closest = norm(I-ray.v0);
                        newPoint = I;
                        intersectFound = 1;
                        closestSurfNum = vertIndex(i);
                    end
                end
            end
            if intersectFound == 1
        % Adds the new points to the path matrix
                ray.v0 = [ray.v0 ; newPoint];
                [ray.v, reflected] = Snell(ray.v, this.n(closestSurfNum,:),ray.n, this.nIndexInside);
%         if reflected == true && InternalReflections == 1
%             rays(ii).colour = 'b';
%         elseif reflected == true && InternalReflections == 0
%             rays(ii).terminated = 1;
%         end
        % ray now get the now index of refraction
                ray.n = this.nIndexOutside;
%         if surface.block == 1
%             rays(ii).terminated = 1;
%         end
%         if surface.countArea == 1;
%            rays(ii).counting = 1; 
%         end
            else
                ray.terminated = 1;
            end
            end
    end 
end

