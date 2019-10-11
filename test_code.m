% add subfolders
addpath(genpath(pwd))
close all
clc




% surf1 = Surf_obj;
% surf2 = Surf_obj;
% 
% surf2.move([10 0 0]);
% 
surf_tree = Surf_tree(surf1);
% figure(1)
% hold on
% axis equal
% surf_tree.plot_surface('all')
surf_tree.surf_refract_index = [1.4 1];
surf_tree.surface_blocking = [1 1];
% surf_tree1 = Surf_tree(surf1);
% surf_tree1.surf_refract_index = [1.4 1.4];
% surf_tree1.surface_blocking = [0 0];

s = Source([1 0 0],[-20 0 0],5, 114, 1);
tic
ray_trace2(s,surf_tree);
toc
% 
% figure(1)
% axis equal
% hold on
% surf1.plot(1,'b')
% surf2.plot(1,'g')
% s.plot_grid(1)
% s.plot(1)
% s.plot_stray(1,10)
% surf_tree.plot_surface('all')

% a = randi(10,20,5);
% [row,col] = find(a==10);

% type = 3;
% gpuFun = @testfun;
% C = gpuArray(A);
% 
% % gputime = gputimeit(@() arrayfun(gpuFun, C(:,1),C(:,2),C(:,3)));
% fprintf(strcat('ERROR:',string(type),' is not an accepted type. \n'))
%                     fprintf(string(type))
%                     fprintf(' is not an accepted type. \n')


% a = zeros(100,3,'single');
% dd = randi(10,100,3);
% 
% cl = cell(100,1,'uint16');
% 
% for i = 1:100
%     cl{i} = dd(i,:);
% end
% 
% b = cell2mat(cl);
% 
% for i =1:100
%     a(i,:) = dd(i,:);
% end
