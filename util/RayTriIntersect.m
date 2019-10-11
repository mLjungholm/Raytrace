function [l, intersect] = RayTriIntersect(Tv0,Tv1,Tv2,N, Rv, v0)

%Check to se if ray intersects surface.
% P0 - point in plane, N - plane normal, v - ray vecor, v0 - ray
% startpoint. 

% New Version: Removing all "dot" and "norm" functions and replacing them
% with pure math.

l = [0,0,0];   % Intersect point.
intersect = 0; % Boolean flag for intersection

u = Tv1-Tv0;
v = Tv2-Tv0;
% N = cross(u,v); % Fix this
% N = [(u(2)*v(3) - u(3)*v(2)) (u(3)*v(1) - u(1)*v(3)) (u(1)*v(2) - u(2)*v(1))];

w0 = v0 - Tv0;
a = -(N(1)*w0(1) + N(2)*w0(2) + N(3)*w0(3));
% a = -dot(N,w0);
b = N(1)*Rv(1) + N(2)*Rv(2) + N(3)*Rv(3);
% b = dot(N,Rv);
r = a/b;

if abs(b) < 10^-5  % If numericlaly close to 0 then parallel
    intersect = 0;
%       intersect = nan;
elseif abs(a) < 10^-5
    intersect = 0;
% intersect = nan;
elseif r < 0
    intersect = 0;
% intersect = nan;
else
    l = v0 + r*Rv;
    
    uu = u(1)^2 + u(2)^2 + u(3)^2;
    vv = v(1)^2 + v(2)^2 + v(3)^2;
    uv = u(1)*v(1) + u(2)*v(2) + u(3)*v(3);
    
%     uu = dot(u,u);
%     uv = dot(u,v);
%     vv = dot(v,v);
    w = l-Tv0;
    
    wu = u(1)*w(1) + u(2)*w(2) + u(3)*w(3);
%     wu = dot(w,u);
    wv = v(1)*w(1) + v(2)*w(2) + v(3)*w(3);
%     wv = dot(w,v);
    
    D = uv*uv - uu*vv;
    S = (uv*wv-vv*wu)/D;
    f = (uv*wu - uu*wv)/D;
    if S < 0 || S > 1
        intersect = 0;
% intersect = nan;
    elseif f < 0 || (S+f) > 1
        intersect = 0;
% intersect = nan;
    else
        intersect = 1;
    end
end
