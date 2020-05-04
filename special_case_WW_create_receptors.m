% Speceial case of creating the receptors for the velvet worm simulation.

% ------------------------------------------------------------------------
% Define the origin points  (oP) for the receptorss using the back retina grid
% coordinates. Cut off coodrinates if needed.

oP = retina_surface.v;

% figure(1)
% hold on
% scatter3(oP(:,1),oP(:,2),oP(:,3))
% lens_surface.plot(1,'b')
% axis equal


% ------------------------------------------------------------------------

% Create special source for tracing.
% 
% key = (oP(:,3) < 5);
% oP = oP(key,:);
% key = (oP(:,3) > 4);
% oP = oP(key,:);
% 
% key = (oP(:,1) < 5);
% oP = oP(key,:);
% key = (oP(:,1) > 0);
% oP = oP(key,:);
% 
% key = (oP(:,2) > -20);
% oP = oP(key,:);
% key = (oP(:,2) < -15);
% oP = oP(key,:);
% oP = oP(oP(:,3) > -5);
S = Source_definedSourceGrid(oP,1,[0,0,0]);


% ------------------------------------------------------------------------

% Trace source from back of retina to lens to find the intesection points
% and extract those.

surf_tree.refract_order = [2,2];
surf_tree.surface_blocking = [1,1];
ray_trace(S,surf_tree);


% ------------------------------------------------------------------------

% retina_surface.plot(1,'g')
lens_surface.plot(1,'b')
S.plot(1)
% surf_tree.plot_tri([1031;1069]);
% surf_tree.plot_bin(1074)
% surf_tree.plot_bin(947)
% surf_tree.plot_bin(917) % parent
% surf_tree.plot_bin(915) % parent
% surf_tree.plot_bin(946)
% surf_tree.plot_bin(948)
% S.plot_stray(1,10)

% figure(1)
% surf_tree.plot_bin(2)
% surf_tree.plot_bin(3)
% 
% surf_tree.plot_bin(4)
xlabel('X')
ylabel('Y')
zlabel('Z')
% surf_tree.plot_bin(5)
% 
% surf_tree.plot_bin(6)
% 
% surf_tree.plot_bin(7)
% 
% surf_tree.plot_bin(8)