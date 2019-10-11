% add subfolders
addpath(genpath(pwd))
close all
clc


% surf1 = Surf_obj;
% surf2 = Surf_obj;
% 
% surf2.move([10 0 0]);
% 
% surf_tree = Surf_tree(surf1,surf2);
% surf_tree.surf_refract_index = [1.4 1.2];
% surf_tree.surface_blocking = [0 1];

s = Source([1 0 0],[-20 0 0], 5, 5, 1);
ray_trace(s,surf_tree);

figure(1)
axis equal
hold on
surf1.plot(1,'b')
surf2.plot(1,'g')
s.plot_grid(1)



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
