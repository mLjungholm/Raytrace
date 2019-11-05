% add subfolders
% addpath(genpath(pwd))

% figure(1)
% hold on
% axis equal
% scatter3(fitted_points(:,1),fitted_points(:,2),fitted_points(:,3),'.b')
% scatter3(points(:,1),points(:,2),points(:,3),'.r')
% scatter3(oP(:,1),oP(:,2),oP(:,3),'.c')
% scatter3(mP(:,1),mP(:,2),mP(:,3),'.b')
% scatter3(aP(:,1),aP(:,2),aP(:,3),'+g')

% scatter3(oP(i,1),oP(i,2),oP(i,3),'k')
% scatter3(mP(i,1),mP(i,2),mP(i,3),'.b')
% scatter3(aP(i,1),aP(i,2),aP(i,3),'k')
% quiver3(aP(i,1),aP(i,2),aP(i,3),cone_dir(i,1),cone_dir(i,2),cone_dir(i,3),8)
% quiver3(aP(i,1),aP(i,2),aP(i,3),pvec(1),pvec(2),pvec(3),8)
% % % quiver3(mP(:,1),mP(:,2),mP(:,3),cone_dir(:,1),cone_dir(:,2),cone_dir(:,3),4)
% scatter3(points(p_ind,1),points(p_ind,2),points(p_ind,3),'r')

close all
clc

start = [-2.5 1 1];
dir = [1 -0.3 -0.4];

gridX = -1:0.005:1;
voxNum = size(gridX,2)^3;

tic
for i = 1:10000
path = voxel_trace([-1 -1 -1],[1 1 1],0.01, start, dir, 0);
end
toc

function [x,y,z] = create_grid(minP, maxP, stepSize)
stepN = floor((maxP-minP)/stepSize);
endP = minP + stepN*stepSize;
endP = [endP; endP + stepSize];
[~, closest] = min(abs(maxP-endP));
[x,y,z] = meshgrid(minP(1):stepSize:endP(closest(1),1),minP(2):stepSize:endP(closest(2),2),minP(3):stepSize:endP(closest(3),3));
end

% Voxel grid traversal algorithm
function path = voxel_trace(boxMin,boxMax, stepSize, start, dir, verbose)

boundingBox = [(boxMin-stepSize./2) (boxMax+stepSize./2)];
boxSize = [boundingBox(4)-boundingBox(1),boundingBox(5)-boundingBox(2),boundingBox(6)-boundingBox(3)];
boxN = round(((boxMax-boxMin)+stepSize)./stepSize);

if (verbose)
    figure;
    hold on;
    axis equal
    text(start(1), start(2), start(3), 'origin');
    plot3(start(1), start(2), start(3), 'k.', 'MarkerSize', 15);
    quiver3(start(1), start(2), start(3), dir(1), dir(2), dir(3), norm((boxMax-boxMin)),'r');
    
    vmin = boundingBox(1:3);
    vmax = boundingBox(4:6);
    BoxVertices = [vmax(1) vmin(2) vmin(3); vmax(1) vmax(2) vmin(3); vmin(1) vmax(2) vmin(3); vmin(1) vmax(2) vmax(3); vmin(1) vmin(2) vmax(3); vmax(1) vmin(2) vmax(3); vmin; vmax ];
    BoxFaces = [1 2 3 7; 1 2 8 6; 1 6 5 7; 7 5 4 3; 2 8 4 3; 8 6 5 4];
    h = patch('Vertices',BoxVertices,'Faces',BoxFaces,'FaceColor','yellow');
    set(h, 'FaceAlpha', 0.1);
    
    view(60,30);
    axis tight;
    xlabel('x')
    ylabel('y');
    zlabel('z');
    grid on;
end

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

% If outside then find first intersection point
if outside
    [tmin, flag] = rayBoxGPU(boundingBox,start, dir);
    if ~flag
        path = nan;
        return
    end
    start = start + tmin*dir;
end

if (verbose)
    plot3(start(1), start(2), start(3), 'r.', 'MarkerSize', 15);
end

x = floor(((start(1)-boundingBox(1))/boxSize(1))*boxN(1)) + 1;
y = floor(((start(2)-boundingBox(2))/boxSize(2))*boxN(2)) + 1;
z = floor(((start(3)-boundingBox(3))/boxSize(3))*boxN(3)) + 1;

if x <= 0; x = 1;
elseif x > boxN(1); x = boxN(1); end
if y <= 0; y = 1;
elseif y > boxN(2); y = boxN(2); end
if z <= 0; z = 1;
elseif z > boxN(3); z = boxN(3); end


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

voxelMaxX  = boundingBox(1) + tVoxelX*boxSize(1);
voxelMaxY  = boundingBox(2) + tVoxelY*boxSize(2);
voxelMaxZ  = boundingBox(3) + tVoxelZ*boxSize(3);

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

voxelSizeX = boxSize(1)/boxN(1);
voxelSizeY = boxSize(2)/boxN(2);
voxelSizeZ = boxSize(3)/boxN(3);

tDeltaX = voxelSizeX/abs(dir(1));
tDeltaY = voxelSizeY/abs(dir(2));
tDeltaZ = voxelSizeZ/abs(dir(3));

maxStep = (max(boxN))^2;
path = zeros(maxStep,3); 
step = 1;
        
while ( (x<=boxN(1)&&(x>=1)) && (y<=boxN(2))&&(y>=1) && (z<=boxN(3))&&(z>=1) )
    path(step,1) = x;
    path(step,2) = y;
    path(step,3) = z;   
    step = step+1;
    
    if (verbose)
        fprintf('Intersection: voxel = [%d %d %d] \n', [x y z]);
        
        t1 = [(x-1)/boxN(1), (y-1)/boxN(2), (z-1)/boxN(3) ];
        t2 = [  (x)/boxN(1),  (y)/boxN(2),    (z)/boxN(3) ];
        
        vmin = (boundingBox(1:3) + t1.*boxSize);
        vmax = (boundingBox(1:3) + t2.*boxSize);
        
        smallBoxVertices = [vmax(1) vmin(2) vmin(3); vmax(1) vmax(2) vmin(3); vmin(1) vmax(2) vmin(3); vmin(1) vmax(2) vmax(3); vmin(1) vmin(2) vmax(3); vmax(1) vmin(2) vmax(3); vmin; vmax ];
        smallBoxFaces    = [1 2 3 7; 1 2 8 6; 1 6 5 7; 7 5 4 3; 2 8 4 3; 8 6 5 4];
        
        h = patch('Vertices', smallBoxVertices, 'Faces', smallBoxFaces, 'FaceColor', 'blue', 'EdgeColor', 'white');
        set(h,'FaceAlpha',0.2);    
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
end
path = path(1:step-1,:);

end
