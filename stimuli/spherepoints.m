% Code to generate roughly equaly spaced point on a sphere. 

% Badly comented
% Right now the coordinates ens up on the positive side of the x-axis
function [sphCord, numP] = spherepoints(sphereNumPoints, radius)

%N = 4000; % Number of points
N = sphereNumPoints;
%r = 1;  % Radius.
r = radius;


points = zeros(N,3);

Ncount = 0; % Counter.
a = 4*pi*r^2/N; 
d = sqrt(a);

M_theta = round(pi/d);
d_theta = pi/M_theta; 
d_phi = a/d_theta;

for m = 0:(M_theta-1)
    theta = pi*(m+0.5)/M_theta;
    M_phi = round(2*pi*sin(theta)/d_phi);
    for n = 0:(M_phi-1)
        phi = 2*pi*n/M_phi;
        points(Ncount+1,:) = [phi, theta, r];
        Ncount = Ncount +1;
    end
end


cartCord = zeros(1,3);
Pcount = 1;
% This part sorts points on the sphere to select points within a sertain
% range
for i=1:N
    theta = points(i,2);
    phi = points(i,1);
    r = points(i,3);
    x = r*sin(theta)*cos(phi);
    y = r*sin(theta)*sin(phi);
    z = r*cos(theta);
%     if (theta >= 10*pi/180 && theta <= 150*pi/180 && phi >= 0*pi/180 && phi <= 80*pi/180)
%         cartCord(Pcount,:) = [x,y,z];
%         Pcount = Pcount+1;
%     end
%     if (theta >= 10*pi/180 && theta <= 150*pi/180 && phi >= 280*pi/180 && phi <= 360*pi/180)
%         cartCord(Pcount,:) = [x,y,z];
%         Pcount = Pcount+1;
%     end
    if (theta >= 0.01*pi/180 && theta <= 180*pi/180 && phi >= 0*pi/180 && phi <= 90*pi/180)
        cartCord(Pcount,:) = [x,y,z];
        Pcount = Pcount+1;
    end
    if (theta >= 0.01*pi/180 && theta <= 180*pi/180 && phi >= 270*pi/180 && phi <= 360*pi/180)
        cartCord(Pcount,:) = [x,y,z];
        Pcount = Pcount+1;
    end
end

sphCordT = cartCord;
numP = size(sphCordT,1);
for i = 1:numP
    [az,el,r] = cart2sph(cartCord(i,1),cartCord(i,2),cartCord(i,3));
    sphCordT(i,:) = [az, el, r];
end
sphCord = sortrows(sphCordT,2);

plot3(cartCord(:,1),cartCord(:,2),cartCord(:,3),'.')
axis equal
end

