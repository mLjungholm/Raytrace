% This function is used in the raytracing algortihtm for the octreee

function [tmin, flag] = rayBoxGPU(boundaries, v0, v)
% Ray-box intersection algorithm of Smit (1998) formatted for arrayfun to 
% allow hardware acceleration:
% Smits, B. (1998). Efficiency issues for ray tracing.  Journal of Graphics
% Tools, 3(2):1–14.
% Ray box intersection is typically used to locate (axis aligned) spatial 
% bins (e.g. octree bins / regular grid) to optimise ray-tracing

% INPUT (scalar):
%    orx, ory, orz: xyz componants of the ray origin
%    Dx, Dy, Dz: xyz components of the ray directional (unit) vectors
%    minx, miny, minz: xyz componants of the box minima 
%    maxx, maxy, maxz: xyz componants of the box maxima
% OUTPUT (scalar:
%    tmin: minimum distance from from the ray-box intersection to the 
%    origin or nan if no intersection is located
%    flag: 1 if intersection / 0 if no intersection

% Use with the arrayfun as a gpuarray object
%
% Usage example: 
% Step 1: convert mx3 direction vectors, D = [Dx Dy Dz] to gpuarray object
% >> gD = gpuArray(D);
% Step 2: call rayBoxGPU using arrayfun with scalar input formatting
% where min, max are the nx3 vertex lists of the box min-max corner points
% and where or is the xyz coordinates of the origin
% >> [tmin, flag] = arrayfun(@rayBoxGPU, min(:,1)', min(:,2)', min(:,3)', ...
%                             max(:,1)', max(:,2)', max(:,3)', ...
%                             or(:,1), or(:,2), or(:,3), ...
%                             gD(:,1),gD(:,2),gD(:,3));
% Step 3: recover data
% distmin = gather(tmin);
% flagBox = gather(flag);
% Output is one mxn array containing a the distance from the ray-box
% intersection to the origin or nan if no intersection is located (distmin)
% and one mxn logical containing flags for ray-box intersections.

% Per-ray flags can be obtained from the output tmin using the following 
% method:
% >> flagB = true(size(D,1),1);
% >> flagB(sum(isnan(tmin),2) == size(min,1)) = false;
% This may save transfer time off the GPU

% Dependencies: requires Parallel Computing Toolbox
% Based upon the implementation by Jesus P. Mena-Chalco
minx = boundaries(1);
miny = boundaries(2);
minz = boundaries(3);
maxx = boundaries(4);
maxy = boundaries(5);
maxz = boundaries(6);
orx = v0(1);
ory = v0(2);
orz = v0(3);
Dx = v(1);
Dy = v(2);
Dz = v(3);

if  Dx >= 0
    tmin = (minx - orx) / Dx;
    tmax = (maxx - orx) / Dx;
elseif Dx == 0
    tmin = -inf;
    tmax = inf;
else
    tmin = (maxx - orx) / Dx;
    tmax = (minx - orx) / Dx;
end

if  Dy > 0
    tymin = (miny - ory) / Dy;
    tymax = (maxy - ory) / Dy;
elseif Dy == 0
    tymin = -inf;
    tymax = inf;
else
    tymin = (maxy - ory) / Dy;
    tymax = (miny - ory) / Dy;
end

if  tmin > tymax || tymin > tmax
    flag = false;
    tmin = nan;
    return
end

if tymin > tmin
   tmin = tymin;
end

if tymax < tmax
   tmax = tymax;
end

if  Dz > 0
    tzmin = (minz - orz) / Dz;
    tzmax = (maxz - orz) / Dz;
elseif Dz == 0
    tzmin = -inf;
    tzmax = inf;
else
    tzmin = (maxz - orz) / Dz;
    tzmax = (minz - orz) / Dz;
end

if tmin > tzmax || tzmin > tmax
   flag = false;
   tmin = nan;
    return
end

if tzmin > tmin
   tmin = tzmin;
end

% ray-box hit = true 
flag = true;

end

