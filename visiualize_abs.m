function visiualize_abs(receptor_id,abs_struct,filename,path,receptor_space)
% temp_var = load(strcat(path, char(filename(abs_struct.file_id(receptor_id)))),'s');
% s = temp_var.s;
[Az, El] = meshgrid(-pi:0.01:pi,-pi/2:0.01:pi/2);
[az,el,r] = cart2sph(-abs_struct.source_coords(:,1),abs_struct.source_coords(:,2),abs_struct.source_coords(:,3));
% % interpolate nonuniformly spaced points


% r1 = ones(size(el,1),1)*120; 
R = griddata(az,el,r,Az,El);
C = griddata(az,el,abs_struct.absorption_mat(:,receptor_id),Az,El).*10^6;
% % convert to cart
[x, y, z] = sph2cart(Az,El,R);
% 
figure(1)
hold on
colormap(viridis)
% axis equal off vis3d
axis equal
grid on
receptor_space.plot_cone(receptor_id,1,'g')
surface(-x,y,z,C,'edgealpha',0.05)
quiver3(1.5,0,0,0.6,0,0,'r','linewidth',1.2)
quiver3(0,1.3,0,0,0.6,0,'b','linewidth',1.2)
quiver3(0,0,1.1,0,0,0.5,'b','linewidth',1.2)
end