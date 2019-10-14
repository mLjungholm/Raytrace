% add subfolders
addpath(genpath(pwd))
close all
clc


folder_name = uigetdir('','Select save folder');
mkdir(folder_name,'tracedRays');
mkdir(folder_name,'absorbedRays');

[pc,d_theta,d_phi] = sphere_points_2d_angular(4.535,150, 1 ,0);


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
    s = Source(dir,start,70, 81, 1);  % 114 size = 10000rays, s81=5022 rays.
    ray_trace(s,surf_tree);
    savename = strcat(folder_name,'\tracedRays\',mat2str(start,4),'_traced.mat');
    save(savename,'s')
end
delete(f)