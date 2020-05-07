% Comment to code.
% This code uses step indices. which is a fast travel algorithm. However if
% we are to calculate absolute distance traveled for each voxel to
% calculate absorption then we might aswell do absolute distance for the
% travel algorithm. 

% Error: This code will not work. if the ray starts inside the volume since
% "tmin" is only determined for outside rays. I think this might be fixed
% by just set "tmin = 0" if inside.

function absorption_values = trace_retina(receptor_volume,source, abs_coeff, volume_id)
% Function for tracing the absorption in volume with volume_id. 

% ! This function assumes that each ray only passes through a absorption
% volume once !

% Determine the bounding volume of the receptor grid.
boundingBox = [(receptor_volume.boxMin-receptor_volume.stepSize./2) (receptor_volume.boxMax+receptor_volume.stepSize./2)];
boxSize = [boundingBox(4)-boundingBox(1),boundingBox(5)-boundingBox(2),boundingBox(6)-boundingBox(3)];
boxN = round(((receptor_volume.boxMax-receptor_volume.boxMin)+receptor_volume.stepSize)./receptor_volume.stepSize);

% Create a zero array for the absorption values of the receptors
absorption_values = zeros(receptor_volume.receptor_nums,1);

% Trace the ray absorption for all rays in parallel
arrayfun(@voxel_trace,1:source.num_rays);

function voxel_trace(ray_ind)
    % Detrmine if the ray traces through the abosption volume in question. 
    if ~ismember(source.absorption(ray_ind,:),volume_id)
        return;
    end
    
    % Find the staring position of the ray. In current version we asume
    % that the ray starts att the second last position.
    
    start = [source.path_x(ray_ind,source.steps(ray_ind-1)),source.path_y(ray_ind,source.steps(ray_ind-1)),source.path_z(ray_ind,source.steps(ray_ind-1))];
    % Determine if ray starts outside or inside the voxel grid
    if start(1) < boundingBox(1) || start(1) > boundingBox(4)
        outside = 1;
    elseif start(2) < boundingBox(2) || start(2) > boundingBox(5)
        outside = 1;
    elseif start(3) < boundingBox(3) || start(3) > boundingBox(6)
        outside = 1;
    else
        outside = 0;
    end
    dir = source.v(ray_ind,:);
    
    % If outside then find first intersection point
    tmin = 0;
    if outside
        [tmin, flag] = rayBoxGPU(boundingBox,start, dir);
        if ~flag
            return
        end
        start = start + tmin*dir;
    end
    
    % This determines the volume x,y,z (voxel) indices that the ray starts
    % in. 
    x = floor(((start(1)-boundingBox(1))/boxSize(1))*boxN(1)) + 1;
    y = floor(((start(2)-boundingBox(2))/boxSize(2))*boxN(2)) + 1;
    z = floor(((start(3)-boundingBox(3))/boxSize(3))*boxN(3)) + 1;
    
    if x <= 0; x = 1;
    elseif x > boxN(1); x = boxN(1); end
    if y <= 0; y = 1;
    elseif y > boxN(2); y = boxN(2); end
    if z <= 0; z = 1;
    elseif z > boxN(3); z = boxN(3); end
    
    % Detrime the direction of the absorption rays.
    % tVoxelX etc is the fraction of voxels the ray has room to travel.
    if (dir(1)>=0)
        tVoxelX = (x)/boxN(1);
        stepX = 1;
    else
        tVoxelX = (x-1)/boxN(1);
        stepX = -1;
    end
    if (dir(2)>=0)
        tVoxelY = (y)/boxN(2);
        stepY = 1;
    else
        tVoxelY = (y-1)/boxN(2);
        stepY = -1;
    end
    if (dir(3)>=0)
        tVoxelZ = (z)/boxN(3);
        stepZ = 1;
    else
        tVoxelZ = (z-1)/boxN(3);
        stepZ = -1;
    end
    
    % Determine the fraction of volume the ray has room to travel.
    voxelMaxX  = boundingBox(1) + tVoxelX*boxSize(1);
    voxelMaxY  = boundingBox(2) + tVoxelY*boxSize(2);
    voxelMaxZ  = boundingBox(3) + tVoxelZ*boxSize(3);
    
    % ERROR. tmin only exists if the ray starts outside. 
    
    % tMaxX etc determines how long (in fractions of a sep) the ray has to
    % travel to reach the wall of the current voxel.
    
    % Using this value we can determine what index will be the next
    % traveled voxel. 
    
    if dir(1) == 0
        tMaxX = inf;
    else
        tMaxX = tmin + (voxelMaxX-start(1))/dir(1);
    end
    if dir(2) == 0
        tMaxY = inf;
    else
        tMaxY = tmin + (voxelMaxY-start(2))/dir(2);
    end
    if dir(3) == 0
        tMaxZ = inf;
    else
        tMaxZ = tmin + (voxelMaxZ-start(3))/dir(3);
    end
    
    % Absolute size of the voxels [m]
    voxelSizeX = boxSize(1)/boxN(1);
    voxelSizeY = boxSize(2)/boxN(2);
    voxelSizeZ = boxSize(3)/boxN(3);
    
    % Absolute stepsize [m/step]
    tDeltaX = voxelSizeX/abs(dir(1));
    tDeltaY = voxelSizeY/abs(dir(2));
    tDeltaZ = voxelSizeZ/abs(dir(3));
    
    % Create a max step fail check
    maxStep = (max(boxN))^2;
    step = 1;
    I = 1;
    
    % Loop runs while x,y,z voxel indices stays within boxNumber,
    while ( (x<=boxN(1)&&(x>=1)) && (y<=boxN(2))&&(y>=1) && (z<=boxN(3))&&(z>=1) )
        %     path(step,1) = x;
        %     path(step,2) = y;
        %     path(step,3) = z;
        
        % Fist check if the current indices has any receptors linked to it.
        if point_inds(x,y,z) ~= 0
            In = exp(-step*abs_coeff);
            absVal = I - In;
            I = In;
%             absorption_values(nonzeros(receptors_inds(point_inds(x,y,z),:))) = absorption_values(nonzeros(receptors_inds(point_inds(x,y,z),:))) + absVal;
            step = step + 1;
        end
        
        if (tMaxX < tMaxY)
            if (tMaxX < tMaxZ)
                x = x + stepX;
                tMaxX = tMaxX + tDeltaX;
            else
                z = z + stepZ;
                tMaxZ = tMaxZ + tDeltaZ;
            end
        else
            if (tMaxY < tMaxZ)
                y = y + stepY;
                tMaxY = tMaxY + tDeltaY;
            else
                z = z + stepZ;
                tMaxZ = tMaxZ + tDeltaZ;
            end
        end
        %     step = step+1;
    end
    % path = path(1:step-1,:);
    
end
end