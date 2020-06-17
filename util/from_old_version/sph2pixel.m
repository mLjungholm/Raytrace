function [yi,xi] = sph2pixel(sphCord,imSize, diam)
% First find the angle theta betwen the vector OP and y-axis


[xp,yp,zp] = sph2cart(sphCord(1), sphCord(2),1);
% OP = [xp,yp,zp];
theta = acosd(dot([xp,yp,zp],[1,0,0]));

% Use theta to find the ri-radius
p = [-2.39756826754156e-05 0.0131523270675485 -0.00163513838631874];
    
ri = (p(1)*theta^2 + p(2)*theta + p(3))*diam/2;

% Find the rotation angle delta betwen OP and the xy-plane
if yp == 0 && zp == 0
    delta = 0;
else
% delta = acos(dot([0,yp,zp],[xp,yp,0])/(norm([0,yp,zp])*norm([xp,yp,0])));
delta = acos(dot([zp,yp],[0,1])/norm([yp,zp]));
end

xs = ri*cos(delta);
ys = ri*sin(delta);

center = zeros(1,2);
if mod(imSize(2),2)
    center(2) = imSize(2)/2 +0.5;
else
    center(2) = imSize(2)/2;
end

if mod(imSize(1),2)
    center(1) = imSize(1)/2 +0.5;
else
    center(1) = imSize(1)/2;
end

xi = round(center(2)+xs);
if sphCord(2) > 0
    yi = round(center(1)-ys); % Becaus lower pixel index means higher upp in the image
else
    yi = round(center(1)+ys);
end
end