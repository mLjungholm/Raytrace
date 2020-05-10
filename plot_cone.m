
% point_a = [1 -5 5];
% point_b = [1 -4 1];
% cone_points(point_a,point_b,1,1,'c');
% pt = [p3dR; p3dRhalf];
% 
% resol = 10;
% dir = point_b - point_a;
% l = sqrt(dir(1)^2 + dir(2)^2 + dir(3)^2); 
% step = 2*pi/resol;
% theta = 0:step:2*pi-step;
% x2d = r*cos(theta);
% y2d = r*sin(theta);
% p3d = [x2d' y2d' zeros(resol,1)];
% ph_dir = -acos(dir(3)/l);
% if ph_dir == 0 || ph_dir == pi
%     th_dir = 0;
% else
%     th_dir = -atan(dir(2)/dir(1));
% end
% rotY =[cos(ph_dir) 0 sin(ph_dir); 0 1 0; -sin(ph_dir) 0 cos(ph_dir)];
% rotZ = [cos(th_dir) -sin(th_dir) 0; sin(th_dir) cos(th_dir) 0; 0 0 1];
% rotT = rotY*rotZ;
% p3dR = p3d*rotT; 
% p3dRhalf = p3dR.*0.5;
% a = p3dR+point_a;
% b = p3dRhalf+point_a + dir;
% 

function plot_cone(point_a,point_b,r, figurnr,color)
dir = point_b - point_a;
l = sqrt(dir(1)^2 + dir(2)^2 + dir(3)^2); 

p3d = [1,0,0;0.809016994374948,0.587785252292473,0;...
    0.309016994374947,0.951056516295154,0;...
    -0.309016994374947,0.951056516295154,0;...
    -0.809016994374947,0.587785252292473,0;...
    -1,1.22464679914735e-16,0;...
    -0.809016994374948,-0.587785252292473,0;...
    -0.309016994374948,-0.951056516295154,0;...
    0.309016994374947,-0.951056516295154,0;...
    0.809016994374947,-0.587785252292473,0];

p3d = p3d*r;

ph_dir = -acos(dir(3)/l);
if ph_dir == 0 || ph_dir == pi
    th_dir = 0;
else
    th_dir = -atan(dir(2)/dir(1));
end
rotY =[cos(ph_dir) 0 sin(ph_dir); 0 1 0; -sin(ph_dir) 0 cos(ph_dir)];
rotZ = [cos(th_dir) -sin(th_dir) 0; sin(th_dir) cos(th_dir) 0; 0 0 1];
rotT = rotY*rotZ;
p3dR = p3d*rotT; 
p3dRhalf = p3dR.*0.5;
p3dR = p3dR+point_a;
p3dRhalf = p3dRhalf+point_a + dir;

faces = [1,2,11;2,3,12;3,4,13;4,5,14;5,6,15;6,7,16;...
    7,8,17;8,9,18;9,10,19;10,11,20;11,12,2;12,13,3;...
    13,14,4;14,15,5;15,16,6;16,17,7;17,18,8;18,19,9;...
    19,20,10;1,11,10];

pt = [p3dR; p3dRhalf];

figure(figurnr)
hold on
axis equal
plot3(point_a(1),point_a(2),point_a(3),'ro')
plot3(point_b(1),point_b(2),point_b(3),'bo')
trisurf(faces,pt(:,1),pt(:,2),pt(:,3),'facecolor',color,'FaceAlpha',0.3,'EdgeAlpha',0.3,'edgecolor',color)
end

% function [a,b] = cone_points_2(point_a,point_b,r)
% resol = 10;
% dir = point_b - point_a;
% l = sqrt(dir(1)^2 + dir(2)^2 + dir(3)^2); 
% step = 2*pi/resol;
% theta = 0:step:2*pi-step;
% x2d = r*cos(theta);
% y2d = r*sin(theta);
% p3d = [x2d' y2d' zeros(resol,1)];
% ph_dir = -acos(dir(3)/l);
% if ph_dir == 0 || ph_dir == pi
%     th_dir = 0;
% else
%     th_dir = -atan(dir(2)/dir(1));
% end
% rotY =[cos(ph_dir) 0 sin(ph_dir); 0 1 0; -sin(ph_dir) 0 cos(ph_dir)];
% rotZ = [cos(th_dir) -sin(th_dir) 0; sin(th_dir) cos(th_dir) 0; 0 0 1];
% rotT = rotY*rotZ;
% p3dR = p3d*rotT; 
% p3dRhalf = p3dR.*0.5;
% a = p3dR+point_a;
% b = p3dRhalf+point_a + dir;
% 
% faces = [1,2,11;2,3,12;3,4,13;4,5,14;5,6,15;6,7,16;...
%     7,8,17;8,9,18;9,10,19;10,11,20;11,12,2;12,13,3;...
%     13,14,4;14,15,5;15,16,6;16,17,7;17,18,8;18,19,9;...
%     19,20,10;1,11,10];
% end