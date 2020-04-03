function ray_trace(source, volume)

orientation_matrix = [0 0 0 5 3 2; 0 0 1 6 4 0; 0 1 0 7 0 4;...
    0 2 3 8 0 0; 1 0 0 0 7 6; 2 0 5 0 8 0; 3 5 0 0 0 8;...
    4 6 7 0 0 0];

% Max step value in case of errors
max_step = 5;

oc_index = [2 6; 1 5];
oc_index(:,:,2) = [4 8; 3 7];

% Run the ray-trace in paralell
arrayfun(@trace_ray,1:source.num_rays);

    % Main ray trace for each ray
    function trace_ray(ray_index)
%         if ray_index == 75
%             disp("test")
%         end
        % Start values for the ray.
        current_step = source.steps(ray_index);
        v0 = [source.path_x(ray_index,current_step), source.path_y(ray_index,current_step), source.path_z(ray_index,current_step)];
        v = source.v(ray_index,:);
        triangle_intersect_bin = 0;
        
        outside = 0;
        if source.v0(ray_index,1) < volume.bin_boundaries(1,1) || source.v0(ray_index,1) > volume.bin_boundaries(1,4)
            outside = 1;
        elseif source.v0(ray_index,2) < volume.bin_boundaries(1,2) || source.v0(ray_index,2) > volume.bin_boundaries(1,5)
            outside = 1;
        elseif source.v0(ray_index,3) < volume.bin_boundaries(1,3) || source.v0(ray_index,3) > volume.bin_boundaries(1,6)
            outside = 1;
        else
            % find start bin
            current_bin = volume.query(v0);
            bin_depth = volume.bin_depths(current_bin);            
%             [tmin, ~] = rayBoxGPU(volume.bin_boundaries(current_bin,:),v0,v);
            bIP = v0;
            current_oc_index = determine_oc_index();
        end
        
        if outside
            current_bin = 1;
            [tmin, flag] = rayBoxGPU(volume.bin_boundaries(current_bin,:),v0,v);
            if flag == 0 % Return if missed volume
                source.status(ray_index) = 'Ray miss';
                return;
            end
            bIP = v0 + v*tmin; % Bin intersection point
            bin_depth = 0;
        end
        
        % Contiue while ray is alive
        while source.ray_alive(ray_index)
            if source.steps(ray_index) > max_step
                fprintf(strcat('ERROR: Something whent wrong with ray "',num2str(ray_index),'" To many steps \n'));
                break
            end
            
            missed = 0;
            jump_step = 0;
            run = 1;
            while run
                if missed == 1 || jump_step == 1
                    n = 3;
                    missed = 0;
                    jump_step = 0;
                else
                    n = next_move(); % Find what to do next
                end
                switch n
                    case 1
                        current_bin = find_lower_bin();
                        bin_depth = bin_depth  + 1;
                    case 2
                        [iP, intersection, vNew, intersect_flag] = find_intersect();
                        if intersection == 1
                            if isequal(intersect_flag, 'miss')
                                source.ray_alive(ray_index) = 0;
                                source.status(ray_index) = 'Ray missed';
                            else
                                v = vNew;
                                v0 = iP;
                                source.v(ray_index,:) = vNew;
                                source.steps(ray_index) = source.steps(ray_index) + 1;
                                source.path_x(ray_index,source.steps(ray_index)) = iP(1);
                                source.path_y(ray_index,source.steps(ray_index)) = iP(2);
                                source.path_z(ray_index,source.steps(ray_index)) = iP(3);
                                triangle_intersect_bin = current_bin;
                                if isequal(intersect_flag, 'blocked')
                                    source.ray_alive(ray_index) = 0;
                                    source.status(ray_index) = 'Ray ended';
                                end
                            end
                            run = 0;
                        else
                            missed = 1;
                        end
                    case 3 % No lower bins & no sufraces in current bin -> move to next bin.
                        [bIP, current_oc_index] = traverse_bin(); % find next bin
                        if current_oc_index == 0    % No more bins in same level -> move up.
                            jump_step = 1;
                            bin_depth = bin_depth  - 1;
                            current_bin = volume.bin_parents(current_bin);
                            if current_bin == 1 % Ray exit volume
                                source.status(ray_index) = 'Ray lost in space';
                                source.ray_alive(ray_index) = 0;
                                run = 0;
                            else
                                temp_inds = volume.bin_childs(:,volume.bin_parents(current_bin));
                                current_oc_index = find(temp_inds == current_bin);
                            end
                        else
                            current_bin = volume.bin_childs(current_oc_index,volume.bin_parents(current_bin));
                        end
                end
            end
        end
        function flag = next_move()
            % Determining the next step. Fist checks if there are lower bins in
            % the current one. Next checks if the bin is empty. If both are
            % empty then move to edge of bin.
            if volume.bin_childs(1,current_bin) ~= 0
                flag = 1;
            elseif volume.bin_triangles(1,current_bin) ~= 0
                if triangle_intersect_bin == current_bin
                    flag = 3;
                else
                    flag = 2;
                end
            else
                flag = 3;
            end
        end
        
        function current_oc_index = determine_oc_index() 
            xm = 1;
            ym = 1;
            zm = 2;
            xh = (volume.bin_boundaries(current_bin,4)-volume.bin_boundaries(current_bin,1))/2;
            yh = (volume.bin_boundaries(current_bin,5)-volume.bin_boundaries(current_bin,2))/2;
            zh = (volume.bin_boundaries(current_bin,6)-volume.bin_boundaries(current_bin,3))/2;
            if (bIP(1)-xh) > volume.bin_boundaries(current_bin,1)
                xm = 2;
            end
            if (bIP(2)-yh) > volume.bin_boundaries(current_bin,2)
                ym = 2;
            end
            if (bIP(3)-zh) > volume.bin_boundaries(current_bin,3)
                zm = 1;
            end
            % octreeIndex [z,x,y]
            current_oc_index = oc_index(zm,xm,ym);
        end
        function next_bin = find_lower_bin()
            % Uses current intersection point with the local logical indices to
            % determine the next lower bin.
            current_oc_index = determine_oc_index();
            next_bin = volume.bin_childs(current_oc_index,current_bin);
        end
        
        function [iP, intersection, vNew, intersect_flag] = find_intersect()
            % Creates a list of triangles in the current bin and checks for any
            % intersections. Calculate refraction if intersection is found.
            triList = nonzeros(volume.bin_triangles(:,current_bin));
            closest = inf;
            closestIp = 0;
            %         triNum = 0;
            intersection = 0;
            vNew = v;
            for i = 1:length(triList)
                Tv0 = volume.points(volume.faces(triList(i),1),:);
                Tv1 = volume.points(volume.faces(triList(i),2),:);
                Tv2 = volume.points(volume.faces(triList(i),3),:);
                            u = Tv1-Tv0;
                            w = Tv2-Tv0;
                            N = cross(u,w);
                            N = N./sqrt(N(1)^2 + N(2)^2 +N(3)^2);
%                             N = volume.face_n(triList(i));
%                 [iP, intersect] = RayTriIntersect(Tv0,Tv1,Tv2,volume.face_n(triList(i),:), v, v0);
                [iP, intersect] = RayTriIntersect(Tv0,Tv1,Tv2,N, v, v0);
                if intersect == 1
                    dist = norm(iP-v0(end,:));
                    if  dist < closest
                        closest = dist;
                        closestIp = iP;
                        triNum = i;
                        intersection = 1;
%                         triN = volume.face_n(triList(i),:);
                          triN = N;
                    end
                end
            end
            iP = closestIp;
            if intersection == 1
                if volume.face_surf_index(volume.bin_triangles(triNum,current_bin)) == volume.refract_order(source.steps(ray_index))
                    if volume.surface_blocking(source.steps(ray_index)) == 0
                        [vNew, ~] = Snell(v, triN, source.refract_index(ray_index), volume.surf_refract_index(source.steps(ray_index)));
                        source.refract_index(ray_index) = volume.surf_refract_index(source.steps(ray_index));
                        intersect_flag = 'refract';
                    else
                        intersect_flag = 'blocked';
                    end
                else
                    intersect_flag = 'miss';
                end
            else
                intersect_flag = 'miss';
            end
        end
        
        function [nextP, nextInd] = traverse_bin()
            % Moves from current point in bin to bin edge. Then cheks if there
            % are any more bins in the same level, otherwise moves up one
            % level. if this is max level, then exits the volume.
            
            %         parenBin = OT.BinParents(current_bin);
            % currentOcIndex. use it to determine next bin or parent borders.
            
            %         tx = 0; ty = 0; tz = 0;
            d = zeros(1,3);
            if v(1) >= 0
                tx = (volume.bin_boundaries(current_bin,4)-bIP(1))/v(1);
                d(1) = 3;
            else
                tx = (volume.bin_boundaries(current_bin,1)-bIP(1))/v(1);
            end
            if v(2) >= 0
                ty = (volume.bin_boundaries(current_bin,5)-bIP(2))/v(2);
                d(2) = 3;
            else
                ty = (volume.bin_boundaries(current_bin,2)-bIP(2))/v(2);
            end
            if v(3) >= 0
                tz = (volume.bin_boundaries(current_bin,6)-bIP(3))/v(3);
                d(3) = 3;
            else
                tz = (volume.bin_boundaries(current_bin,3)-bIP(3))/v(3);
            end
            [t, I] = min([tx,ty,tz]);
            try
            nextInd = orientation_matrix(current_oc_index,I+d(I));
            catch
                disp(strcat("ERROR: current_oc_index not defined in ray", string(ray_index)));
            end
            %         [nexOcInd, oob] = Orientation();
            nextP = bIP + v*t;
        end
        
        %     function [nexOcInd, oob] = Orientation()
        %         % Uses the orientation matrix to determine what the next OcTree
        %         % index and if the ray exits the parent bin. (oob = out of bounds)
        %     end
    end
end