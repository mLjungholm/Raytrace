% This is the main raytracing loop. The rays will be run on parallel cores
% if posible to speed up the simulation but might still take quite some
% time depending on complexity of the model. 

% The files needed for the simulation are the generated "surf_tree" ie. the
% space partitioned volume of all surfaces. You will also need a list of
% point source origins "pc". 

% The program will tell you to select a folder for saving the resulting
% trace files. One file for each points source. 

% add subfolders
addpath(genpath(pwd))
close all
clc

folder_name = uigetdir('','Select save folder');
mkdir(folder_name,'tracedRays');
% mkdir(folder_name,'absorbedRays');

% [pc,d_theta,d_phi] = sphere_points_2d_angular(4.535,150, 1 ,0);

f = waitbar(0,'1','Name','Tracing rays...',...
    'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
setappdata(f,'canceling',0);

source_nums = size(pc,1);

for source_num = 1:source_nums
    if getappdata(f,'canceling')
        break
    end
    waitbar(source_num/source_nums,f,sprintf('tracing...'))
    
    start = pc(source_num,:);
    dir = -start;
    s = Source(dir,start,60,114,1);  % 114 size = 10000rays, s81=5022 rays. 
    s.source_id = source_num;
    
    % Raytracing
    ray_trace_single(s,cornia_tree);
    s.unblock;
    ray_trace_single(s,lens_tree);    
    
    % Absorption
    receptor_space.absorbed_val(:) = 0;
    absorption_trace(receptor_space,s);
    s.absorption = receptor_space.absorbed_val;
    
    savename_raytrace = strcat(folder_name,'\tracedRays\',mat2str(start,4),'_traced.mat');
    save(savename_raytrace,'s')
end
delete(f)