% This version of spherepoints uses a set angle or frequency to place
% points evenly spaced in longitudinal and latitudinal direction.

function [pc,d_theta,d_phi] = sphere_points_2d_angular(ang, radius, half_sphere ,plot)
% full = false;
% ang = 10; % Set the angular spacing.
ang_rad = 2*pi/360*ang; % Transform into radians.
N = round((4*pi)/ang_rad^2);
r = 1; % There is no point of the radius beeing anything but 1.

points = zeros(N,3); % Allocate a matrix for the coodrinates.
% pc = zeros(N,3); % Cartesian coordinates.

Ncount = 0; % Counter.
a = 4*pi*r^2/N; 
d = sqrt(a);

M_theta = round(pi/d);
d_theta = pi/M_theta; 
d_phi = a/d_theta;

pc = zeros(N,3);
for m = 0:(M_theta-1)
    theta = pi*(m+0.5)/M_theta;
    M_phi = round(2*pi*sin(theta)/d_phi);
    for n = 0:(M_phi-1)
        phi = 2*pi*n/M_phi;
        points(Ncount+1,:) = [phi, theta, r];
%         [x,y,z] = sph2cart(theta,phi,r);
        x = radius*sin(theta)*cos(phi);
        y = radius*sin(theta)*sin(phi);
        z = radius*cos(theta);
        pc(Ncount+1,:) = [x,y,z];
        Ncount = Ncount +1;
    end
end

if half_sphere
    inds = pc(:,1) < 0;
    pc = pc(inds,:);
end

if plot
    figure(1)
    scatter3(pc(:,1),pc(:,2),pc(:,3),'.')
    axis equal
end
end
