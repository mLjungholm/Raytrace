point_a = [1 1 1];
point_b = [3 2 2];
r = 10;
r2 = r/2;

resol = 100;

dir = point_b - point_a;
l = sqrt(dir(1)^2 + dir(2)^2 + dir(3)^2); 
% l = norm(dir);
dir_n = dir./l;

step = 2*pi/resol;
theta = 0:step:2*pi-step;

x2d = r*cos(theta);
y2d = r*sin(theta);
p3d = [x2d' y2d' zeros(resol,1)];
% plot(y,x)


ph_dir = acos(dir(3)/l);
if ph_dir == 0 || ph_dir == pi
    th_dir = 0;
else
    th_dir = atan(dir(2)/dir(1));
end



% rotX = [1 0 0; 0 cos(ph_dir) -sin(ph_dir); 0 sin(ph_dir) cos(ph_dir)];
rotY =[cos(ph_dir) 0 sin(ph_dir); 0 1 0; -sin(ph_dir) 0 cos(ph_dir)];
rotZ = [cos(th_dir) -sin(th_dir) 0; sin(th_dir) cos(th_dir) 0; 0 0 1];
rotT = rotY*rotZ;

p3dR = p3d*rotT; 
p3dRhalf = p3dR.*0.5;

% p3dR = p3dR+point_a;
% p3dRhalf = p3dRhalf+point_a + dir;


figure(1)
hold on
axis equal
quiver3(0,0,0,dir_n(1),dir_n(2),dir_n(3),l)
plot3(point_a(1),point_a(2),point_a(3),'ro')
plot3(point_b(1),point_b(2),point_b(3),'bo')
% scatter3(p3d(:,1),p3d(:,2),p3d(:,3),'r.')
scatter3(p3dR(:,1),p3dR(:,2),p3dR(:,3),'b.')
scatter3(p3dRhalf(:,1),p3dRhalf(:,2),p3dRhalf(:,3),'r.')

function [a,b] = cone_points(point_a,point_b,r)
resol = 10;

dir = point_b - point_a;
l = norm(dir);
dir_n = dir./l;

step = 2*pi/resol;
theta = 0:step:2*pi-step;

x2da = r*cos(theta);
y2da = r*sin(theta);
p3d = [x2da' y2da' zeros(resol,1)];


ph_dir = acos(dir(3)/l);
if ph_dir == 0 || ph_dir == pi
    th_dir = 0;
else
    th_dir = atan(dir(2)/dir(1));
end

rotY =[cos(ph_dir) 0 sin(ph_dir); 0 1 0; -sin(ph_dir) 0 cos(ph_dir)];
rotZ = [cos(th_dir) -sin(th_dir) 0; sin(th_dir) cos(th_dir) 0; 0 0 1];
rotT = rotZ*rotY;

p3dR = p3d*rotT; 
p3dRhalf = p3dR.*0.5;
end