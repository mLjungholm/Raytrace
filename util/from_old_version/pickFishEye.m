% Function to pick out values at seartain spherical coordinates from a
% fisheye lens image. Returns pixel coordinets [row,cool]. 

% theta is the azimutal angle 
function [x,y] = pickFishEye(phi,theta, size)
% x = (theta+1.0153)*r/pi;
% y = (phi+1.3963)*2*sqrt(r^2-x^2)/pi;
r=(size-1)/2;
x = (pi/2-theta)*2*r/pi;
d = sqrt(r^2-(r-x)^2);
y = (pi/2-phi)*2*d/pi + (r-d);
x = round(x) + 1;
y = round(y) + 1;
% x = x +1;
% y = y +1;
end