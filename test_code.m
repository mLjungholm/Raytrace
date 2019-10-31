% add subfolders
% addpath(genpath(pwd))
% close all
% clc

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



% minP = min([mP;oP],[],1);
% maxP = max([mP;oP],[],1);
minP = [-1 -1 -1];
maxP = [1 1 1];
step = 0.3;



[x,y,z] = create_grid(minP, maxP, step);
points = [reshape(x,[],1),reshape(y,[],1),reshape(z,[],1)];
points = single(points);



function [x,y,z] = create_grid(minP, maxP, step_size)
[x,y,z] = meshgrid(minP(1):step_size:maxP(1),minP(2):step_size:maxP(2),minP(3):step_size:maxP(3));
end