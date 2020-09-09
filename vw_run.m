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

cornia.surf_refract_index = 1.4;
cornia.surface_blocking = 1;
cornia.refract_order = 1;
cornia.surf_absorbing = 0;

lens.surf_refract_index = [1.44 1.36];
lens.surface_blocking = [0 0];
lens.refract_order = [1 1];
lens.surf_absorbing = [0 1];

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
cornia.surface_blocking = 1;
cornia.surf_refract_index = 1.4; % 1.4
cornia.refract_order = 1;
cornia.surf_absorbing = 0;
lens.surf_refract_index = [0 1.44 1.36]; %1.36
lens.surface_blocking = [0 0 1];
lens.refract_order = [1 1 1];
lens.surf_absorbing = [0 0 0];
retina.surface_blocking = [0 0 0 1];
retina.surf_refract_index = [0 0 0 1.36];
retina.refract_order = [0 0 0 1];
retina.surf_absorbing = [0 0 0 0];
% % 
% start = [-40,0,-60];
% start = [-80,0,0];
d = 80;
start = [-d*cosd(55),0,-d*sind(55)];
dir = -start;
s = Source(dir,start,55,50,1);
ray_trace_single(s,cornia);
s.unblock;
ray_trace_single(s,lens);
s.unblock;
ray_trace_single(s,retina);

% d = 80;
% start2 = [-d*cosd(55),0,-d*sind(55)];
% dir2 = -start2;
% s2 = Source(dir2,start2,50,15,1);
% ray_trace_single(s2,cornia);
% s2.unblock;
% ray_trace_single(s2,lens);
% s2.unblock;
% ray_trace_single(s2,retina);
% 
figure(1)
hold on
axis equal
% cornia2.plot(1,'g')
% lens_surface.plot(1,'b')
% retina_surface.plot(1,'w')
% pigment.plot(1,'w')
trisurf(cornia2.f,cornia2.v(:,1),cornia2.v(:,2),cornia2.v(:,3),'Facecolor','b','FaceAlpha',0.2,'EdgeAlpha',0.1)
trisurf(lens_surface.f,lens_surface.v(:,1),lens_surface.v(:,2),lens_surface.v(:,3),'Facecolor','b','FaceAlpha',0.3,'EdgeAlpha',0.1)
trisurf(retina_surface.f,retina_surface.v(:,1),retina_surface.v(:,2),retina_surface.v(:,3),'Facecolor','y','FaceAlpha',0.1,'EdgeAlpha',0.1)
trisurf(pigment.f,pigment.v(:,1),pigment.v(:,2),pigment.v(:,3),'Facecolor','w','FaceAlpha',0.1,'EdgeAlpha',0.1)
s.plot_spec(1,'r',0.5)
% s2.plot_spec(1,'b',0.5)
% s.plot_stray(1,40)

%%
figure(1)
hold on
% axis equal
axis equal off vis3d
% 
% cornia2.ScaleSurf(2);
% lens_surface.ScaleSurf(2);
% retina_surface.ScaleSurf(2);
% pigment.ScaleSurf(2);
trisurf(cornia_surface.f,cornia_surface.v(:,1),cornia_surface.v(:,2),cornia_surface.v(:,3),'Facecolor','b','FaceAlpha',0.2,'EdgeAlpha',0.3)
trisurf(lens_surface.f,lens_surface.v(:,1),lens_surface.v(:,2),lens_surface.v(:,3),'Facecolor','b','FaceAlpha',0.3,'EdgeAlpha',0.1)
trisurf(retina_surface.f,retina_surface.v(:,1),retina_surface.v(:,2),retina_surface.v(:,3),'Facecolor','y','FaceAlpha',0.1,'EdgeAlpha',0.2)
trisurf(pigment.f,pigment.v(:,1),pigment.v(:,2),pigment.v(:,3),'Facecolor','w','FaceAlpha',0.1,'EdgeAlpha',0.1)

quiver3(-50, 110, 0, -15, 0, 0,3, 'LineWidth',2, 'color','r')
quiver3(-50, 110, 0, 0, 15, 0,3, 'LineWidth',2, 'color','k')
quiver3(-50, 110, 0, 0, 0, -15,3, 'LineWidth',2, 'color','b')

a = -pi : pi/2 : pi;                                % Define Corners
ph = pi/4;                                          % Define Angular Orientation (‘Phase’)
x = [cos(a+ph); cos(a+ph)]/cos(ph).*6 -56;
y = [sin(a+ph); sin(a+ph)]/sin(ph).*6 +116;
z = [-ones(size(a)); ones(size(a))].*6 -6;
% figure
surf(x, y, z, 'FaceColor','k', 'facealpha', 0.1)                      % Plot Cube
% hold on
% patch(x', y', z', 'k', 'facealpha)                              % Make Cube Appear Solid
% hold off
% axis([ -1  1    -1  1    -1  1]*1.5)
% grid on
receptor_space.plot_cone(1043,1,'r')
receptor_space.plot_cone(1041,1,'r')
receptor_space.plot_cone(1073,1,'r')


%% Find receptor endpoints
oP = retina_surface.v;
key = oP(:,1) > 25;
oP = oP(key,:);
retina_source = Source_definedSourceGrid(oP, 1, [0,0,0]);
lens.surf_refract_index = 1;
lens.surface_blocking = 1;
lens.refract_order = 1;
lens.surf_absorbing = 0;
ray_trace_single(retina_source,lens);
mP = retina_source.export_endpoints;


figure(1)
retina_surface.plot(1,'g')
lens_surface.plot(1,'b')
cornia_surface.plot(1,'y')
retina_source.plot(1)

%% Check size of receptor base  ~ 4.8

dists = zeros(retina_surface.f_count,1);
for i = 1:retina_surface.f_count
    p1 = retina_surface.v(retina_surface.f(i,1),:);
    p2 = retina_surface.v(retina_surface.f(i,2),:);
    p3 = retina_surface.v(retina_surface.f(i,3),:);
    d1 = norm(p2-p1);
    d2 = norm(p3-p1);
    d3 = norm(p3-p2);
    dists(i) = (d1+d2+d3)/3;
end
hist(dists,100)

%% Create retina space

receptor_space = Receptor_space(oP,mP,[3.5,0.7]);
% receptor_space.allocate_space([200,200,200]);
% receptor_space.fit_receptors2bins;
figure(1)
retina_surface.plot(1,'g')
lens_surface.plot(1,'b')
receptor_space.plot_cone(2588,1,'r')
% receptor_space.plot_cone(2576,1,'r')
% receptor_space.plot_cone(2546,1,'r')
receptor_space.plot_cone(3519,1,'r')
receptor_space.plot_cone(1657,1,'r')

%% Test absorption

receptor_space.absorbed_val(:) = 0;
receptor_space.absorption_coeff = 0.0067;
receptor_space.volume_id = 3;
%%
absorption_trace(receptor_space,s);
absVals = receptor_space.absorbed_val;
absVals = absVals./max(absVals);
[Az, El] = meshgrid(-pi:0.01:pi,-pi/2:0.01:pi/2);
% interpolate nonuniformly spaced points
[az,el,r] = cart2sph(receptor_space.base_pos(:,1),receptor_space.base_pos(:,2),receptor_space.base_pos(:,3));
C = griddata(az,el,absVals,Az,El);
R = griddata(az,el,r,Az,El);
C = C.*1000;
% convert to cart
[x, y, z] = sph2cart(Az,El,R);
% 
figure(1)
hold on
s.plot(1)
cornia.plot(1,'y')
lens.plot(1,'b')
retina.plot(1,'g')
% s.plot_stray(1,20)
% st.plot_surface('all')
% colormap(inferno)
% axis equal off vis3d
axis equal
surface(x,y,z,C,'edgealpha',0.05)

%% Speed test trace.

cornia.surface_blocking = 1;
cornia.surf_refract_index = 1.4; % 1.4
lens.surf_refract_index = [0 1.44 1.36]; %1.36
lens.surface_blocking = [0 0 1];
lens.refract_order = [1 1 1];
lens.surf_absorbing = [0 0 1];

start = [-100,0,0];
dir = -start;
s = Source(dir,start,60,11,1);
tic
ray_trace_single(s,cornia);
s.unblock;
ray_trace_single(s,lens);
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
