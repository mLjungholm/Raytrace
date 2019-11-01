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




% clear
close all
clc

[x,y,z] = create_grid([-1 -1 -1], [1 1 1], 0.4);
% points = [reshape(x,[],1),reshape(y,[],1),reshape(z,[],1)];
% points = single(points);


function [x,y,z] = create_grid(minP, maxP, stepSize)
stepN = floor((maxP-minP)/stepSize);
endP = minP + stepN*stepSize;
endP = [endP; endP + stepSize];
[~, closest] = min(abs(maxP-endP));
[x,y,z] = meshgrid(minP(1):stepSize:endP(closest(1),1),minP(2):stepSize:endP(closest(2),2),minP(3):stepSize:endP(closest(3),3));
end

% Voxel grid traversal algorithm
function path = voxel_trace(xGrid,yGrid,zGrid, start, dir)

path = zeros(10,3); 
path(1,:) = start;
pathInd = 2;

voxelSize = xGrid(1,2,1)-xGrid(1,1,1); % Determine step size of the uniform grid.
boundingBox = [xGrid(1,1,1),yGrid(1,1,1),zGrid(1,1,1),xGrid(1,end,1),yGrid(end,1,1),zGrid(1,1,end)];
boxSize = [boundingBox(4)-boundingBox(1),boundingBox(5)-boundingBox(2),boundingBox(6)-boundingBox(3)];
boxN = [size(xGrid,2),size(yGrid,2),size(zGrid,2)];

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
    [flag, tmin] = rayBoxGPU(boundingBox,start, dir);
    if ~flag
        path = nan;
        return
    end
    start = start + tmin*dir;
    path(pathInd,:) = start;
    pathInd = pathInd + 1;
end

x = floor(((start(1)-boundingBox(1))/boxSize(1))*boxN(1)) + 1;
y = floor(((start(2)-boundingBox(2))/boxSize(2))*boxN(2)) + 1;
z = floor(((start(3)-boundingBox(3))/boxSize(3))*boxN(3)) + 1;

if x <= 0
    x = 1;
elseif x > boxN(1)
    x = boxN(1);
end

if y <= 0
    y = 1;
elseif y > boxN(2)
    y = boxN(2);
end

if z <= 0
    z = 1;
elseif z > boxN(3)
    z = boxN(3);
end




end
