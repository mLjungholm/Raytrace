% cornia_surface.rotate('y', [0,0,0], 180)
% lens_surface.rotate('y', [0,0,0], 180)
% retina_surface.rotate('y', [0,0,0], 180)

% cornia_surface.plot(1,'g')
% lens_surface.plot(1,'g')
% retina_surface.plot(1,'g')


% st = Surf_tree(cornia_surface,lens_surface,retina_surface);
% st.surf_refract_index = [1.35,1.4,1.37];
% st.surface_blocking = [0 0 1];
% st.refract_order = [1,2,2,3];
% st.surf_absorbing = [0 0 1 0];


% rs = Receptor_space(oP,mP,[4,0.7]);
% rs.allocate_space([20,20,20]);
% rs.fit_receptors2bins;
% rs.volume_id = 3;
% rs.absorption_coeff = 0.005;


start = [-20,-10,0];
dir = -start;
dir = dir./norm(dir);
s = Source(dir, start,17,25,1);
ray_trace(s,st);
% 
% s.plot_grid(1)
s.plot(1)
% % s.plot_stray(1,20)
st.plot_surface('all')
% 
absorption_trace(rs,s);


%%
absVals = rs.absorbed_val;
absVals = absVals./max(absVals);
[Az, El] = meshgrid(-pi:0.01:pi,-pi/2:0.01:pi/2);
% interpolate nonuniformly spaced points
[az,el,r] = cart2sph(rs.base_pos(:,1),rs.base_pos(:,2),rs.base_pos(:,3));
C = griddata(az,el,absVals,Az,El);
R = griddata(az,el,r,Az,El);
C = C.*1000;
% convert to cart
[x, y, z] = sph2cart(Az,El,R);

figure(1)
hold on
s.plot(1)
% s.plot_stray(1,20)
st.plot_surface('all')
% colormap(inferno)
% axis equal off vis3d
axis equal
surface(x,y,z,C,'edgealpha',0.05)
% retina_surface.plot(1,'g')
%%

absVals = rs.absorbed_val;
absVals = absVals./max(absVals);

[X,Y]= meshgrid(-21:0.5:21,-21:0.5:21);
C = griddata(rs.base_pos(:,1),rs.base_pos(:,2),absVals,X,Y);
imagesc(C)
% 
%%
% 
absVals = rs.absorbed_val;
absVals = absVals./max(absVals);
bas_sph = zeros(rs.receptor_nums,2);
for i = 1:rs.receptor_nums 
[az,el,~] = cart2sph(rs.base_pos(i,1),rs.base_pos(i,2),rs.base_pos(i,3));
bas_sph(i,:) = [az,el];
end
bas_sph = bas_sph.*180/pi;
[AZ,EL]= meshgrid(-100:0.5:100,-100:0.5:100);
C = griddata(bas_sph(:,1),bas_sph(:,2),absVals,AZ,EL);
figure(1)
imagesc(C)
% colormap(plasma)
% v = [0.015:0.01:0.20];
% contourf(AZ,EL,C,v,'ShowText','off')
% % title('Azimuthal resolution [cycles/deg]')
% xlabel('Azimuthal angle [deg]')
% ylabel('Elevation angle [deg]')
% grid on
% axis equal