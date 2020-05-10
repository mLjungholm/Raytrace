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
% oP = [oP(60,:);oP(62,:)];
S = Source_definedSourceGrid(oP,1,[0,0,0]);
% 
% 
% % ------------------------------------------------------------------------
% 
% % Trace source from back of retina to lens to find the intesection points
% % and extract those.
% 
surf_tree.refract_order = [2,2];
surf_tree.surface_blocking = [1,1];
ray_trace(S,surf_tree);

mP = [S.path_x(:,2),S.path_y(:,2),S.path_z(:,2)];

% ------------------------------------------------------------------------
% % tri_list [371;1527;518;1522;517;1638;131]
% retina_surface.plot(1,'g')
% surf_tree.plot_bin(507)
% surf_tree.plot_tri([371;1527;518;1522;517;1638;131])
% % lens_surface.plot(1,'b')
% S.plot_stray(1,10)
% xlabel('X')
% ylabel('Y')
% zlabel('Z')
