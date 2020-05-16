% cornia_surface.rotate('y', [0,0,0], 180)
% lens_surface.rotate('y', [0,0,0], 180)
% retina_surface.rotate('y', [0,0,0], 180)

% ro = 90;
% cornia.rotate('z',[0,0,0],ro);
% lens.rotate('z',[0,0,0],ro);
% retina.rotate('z',[0,0,0],ro);

% cornia.plot(1,'g')
% lens.plot(1,'g')
% retina.plot(1,'g')

% cornia_tree = Surf_tree(cornia);
% lens_tree = Surf_tree(lens);
% retina_tree = Surf_tree(retina);

% cornia_tree.surf_refract_index = 1.4;
% cornia_tree.surface_blocking = 1;
% cornia_tree.refract_order = 1;
% cornia_tree.surf_absorbing = 0;

% lens_tree.surf_refract_index = [1.44 1.36];
% lens_tree.surface_blocking = [0 0];
% lens_tree.refract_order = [1 1];
% lens_tree.surf_absorbing = [0 1];

% st = Surf_tree(cornia,lens,retina);
% st.surf_refract_index = [1.4,1.44,1.36];
% st.surface_blocking = [0 0 0 1];
% st.refract_order = [1,2,2,3];
% st.surf_absorbing = [0 0 1 0];

% st.plot_surface('all')
% s = Source([1,0,0],[-100,0,0],60,357,1);
% s.plot_grid(1)
% tic
% ray_trace(s,st);
% toc
% s.plot(1);

% sg = Source_definedSourceGrid(retina.v, 1, [0,0,0]);
% ray_trace(sg,st);
% mP = sg.export_endpoints;
% oP = retina.v;
% 
% figure(1)
% hold on
% axis equal
% scatter3(mP(:,1),mP(:,2),mP(:,3),'.');
% scatter3(oP(:,1),oP(:,2),oP(:,3),'.');

% rs = Receptor_space(oP,mP,[4,0.7]);
% rs.allocate_space([20,20,20]);
% rs.fit_receptors2bins;
% rs.volume_id = 3;
% rs.absorption_coeff = 0.005;


%% Trace from source to end of lens
% cornia_tree.surface_blocking = 1;
% cornia_tree.surf_refract_index = 1.4; % 1.4
% lens_tree.surf_refract_index = [0 1.44 1.36]; %1.36
% lens_tree.surface_blocking = [0 0 1];
% lens_tree.refract_order = [1 1 1];
% lens_tree.surf_absorbing = [0 0 1];
% 
start = [-100,0,50];
dir = -start;
s = Source(dir,start,60,25,1);
ray_trace_single(s,cornia_tree);
s.unblock;
ray_trace_single(s,lens_tree);

% figure(1)
% cornia.plot(1,'g')
% s.plot(1)


%% Find receptor endpoints
% oP = retina.v;
% key = oP(:,1) > 40;
% oP = oP(key,:);
% retina_source = Source_definedSourceGrid(oP, 1, [0,0,0]);
% lens_tree.surf_refract_index = 1;
% lens_tree.surface_blocking = 1;
% lens_tree.refract_order = 1;
% lens_tree.surf_absorbing = 0;
% ray_trace_single(retina_source,lens_tree);
% mP = retina_source.export_endpoints;

% figure(1)
% retina.plot(1,'g')
% lens.plot(1,'b')
% cornia.plot(1,'y')
% retina_source.plot(1)

%% Check size of receptor base  ~ 4.8

% dists = zeros(retina.f_count,1);
% for i = 1:retina.f_count
%     p1 = retina.v(retina.f(i,1),:);
%     p2 = retina.v(retina.f(i,2),:);
%     p3 = retina.v(retina.f(i,3),:);
%     d1 = norm(p2-p1);
%     d2 = norm(p3-p1);
%     d3 = norm(p3-p2);
%     dists(i) = (d1+d2+d3)/3;
% end
% hist(dists,100)

%% Create retina space

% receptor_space = Receptor_space(oP,mP,[4.7,0.7]);
% receptor_space.allocate_space([100,100,100]);
% receptor_space.fit_receptors2bins;
% figure(1)
% retina.plot(1,'g')
% lens.plot(1,'b')
% receptor_space.plot_cone(1500,1,'r')

%% Test absorption

% receptor_space.absorbed_val(:) = 0;
% receptor_space.absorption_coeff = 0.0067;
% receptor_space.volume_id = 3;
% absorption_trace(receptor_space,s);
% absVals = receptor_space.absorbed_val;
% absVals = absVals./max(absVals);
% [Az, El] = meshgrid(-pi:0.01:pi,-pi/2:0.01:pi/2);
% % interpolate nonuniformly spaced points
% [az,el,r] = cart2sph(receptor_space.base_pos(:,1),receptor_space.base_pos(:,2),receptor_space.base_pos(:,3));
% C = griddata(az,el,absVals,Az,El);
% R = griddata(az,el,r,Az,El);
% C = C.*1000;
% % convert to cart
% [x, y, z] = sph2cart(Az,El,R);
% % 
% figure(1)
% hold on
% s.plot(1)
% cornia.plot(1,'y')
% lens.plot(1,'b')
% retina.plot(1,'g')
% % s.plot_stray(1,20)
% % st.plot_surface('all')
% % colormap(inferno)
% % axis equal off vis3d
% axis equal
% surface(x,y,z,C,'edgealpha',0.05)

%% Speed test trace.

cornia_tree.surface_blocking = 1;
cornia_tree.surf_refract_index = 1.4; % 1.4
lens_tree.surf_refract_index = [0 1.44 1.36]; %1.36
lens_tree.surface_blocking = [0 0 1];
lens_tree.refract_order = [1 1 1];
lens_tree.surf_absorbing = [0 0 1];

start = [-100,0,0];
dir = -start;
s = Source(dir,start,60,11,1);
tic
ray_trace_single(s,cornia_tree);
s.unblock;
ray_trace_single(s,lens_tree);
toc

% s2 = Source(dir,start,60,25,1);
% tic
% ray_trace(s2,st);
% toc
tic
receptor_space.absorbed_val(:) = 0;
receptor_space.absorption_coeff = 0.0067;
receptor_space.volume_id = 3;
absorption_trace(receptor_space,s);
toc
