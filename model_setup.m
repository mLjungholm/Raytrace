% Flowchart for importing 3d files and setting up the model that will be
% used in the ray tracing simulation.

%------------------------------------------------------------------------

% Step 1. Importing 3d models. This assumes that the model of the eye has
% been separated into ".stl-files". These are imported by creating and 
% naming a Surf_obj. 
%------------------------------------------------------------------------

% Step 2. Checking and orienting surfaces
% Check the position and orientation of all surfaces. This can be done by
% using the built in funtions of the Surf_obj. Check Surf_obj documentation
% for further information on utility funtions.

%                               ! Note !
%       Here it is importaint to think about where the light source is
%       supposed to originate from. Matlab math conversions makes the x 
%       axis most suitable for the purpose of main propagation axis. 
%       Traditionaly in optics the z axis is ususaly chosen as the main 
%       propagation axis. This might however cause problems in the code.

%------------------------------------------------------------------------

% Step 3. Create surface partition. Use the Surf_tree object to partition
% the 3d space adding all the relevant surfaces. This object will later be
% used in the actuall ray tracing. 

% Here you will also need to add refractive index values to the surfaces,
% define if a surface is a blocking surface and the refraction order if
% any. See Surf_tree object for further documentation (unfinished)

%                           ! IMPORTANT !
%       Note that this is computational heavy and might take upp to 30min
%       depending on computer and size of the models. Do not do this part 
%       before you have finnished part 2. Remember to save the file after
%       compleation to avoid having to reapeat this process.

%------------------------------------------------------------------------

% Step 4. Save the created surf_tree file. This will be used in the main
% raytracing simulation. 

%------------------------------------------------------------------------
%   Ceated by Mikael Ljungholm 
%   2020-03-02

%   Dependencies:
%       Surf_obj.m
%       Surf_tree.m


%% Step 1. Import 3d models

% Example 
% ------------------------------------------------------------
% Exmple ".stl" surfaces can be found in /data/test_data/
% cornia_surface = Surf_obj;
% lens_surface = Surf_obj;
% retina_surface = Surf_obj;
% pigment = Surf_obj;
% ------------------------------------------------------------


%% Step 2. Checking and orienting surfaces

% Example
% ------------------------------------------------------------
% cornia2.plot(1,'r')
% lens_surface.plot(1,'b')
% retina_surface.plot(1,'g')
% pigment.plot(1,'k')

% v = [-213, -173,0];
% cornia2.move(v);
% 
% cornia2.plot(1,'r')
% lens_surface.plot(1,'b')
% retina_surface.plot(1,'g')
% pigment.plot(1,'k')

% lens_surface.move(v);
% retina_surface.move(v)
% pigment.move(v);
% cornia2.rotate('x', [0 0 0], 180) % 90 degree rotation around the x axis at position [0 0 0]
% lens_surface.rotate('x', [0 0 0], 180)
% retina_surface.rotate('x', [0 0 0], 180)
% 
% cornia2.plot(1,'r')
% lens_surface.plot(1,'b')
% retina_surface.plot(1,'g')
% pigment.plot(1,'k')
% ------------------------------------------------------------

%% Step 3. Create the surface partition

% Example
%------------------------------------------------------------------------
cornia = Surf_tree(cornia2);
lens = Surf_tree(lens_surface);
retina = Surf_tree(retina_surface);

% trisurf(cornia2.f,cornia2.v(:,1),cornia2.v(:,2),cornia2.v(:,3),'Facecolor','b','FaceAlpha',0.2,'EdgeAlpha',0.1)
% trisurf(lens_surface.f,lens_surface.v(:,1),lens_surface.v(:,2),lens_surface.v(:,3),'Facecolor','b','FaceAlpha',0.3,'EdgeAlpha',0.1)
% trisurf(retina_surface.f,retina_surface.v(:,1),retina_surface.v(:,2),retina_surface.v(:,3),'Facecolor','y','FaceAlpha',0.1,'EdgeAlpha',0.1)
% trisurf(pigment.f,pigment.v(:,1),pigment.v(:,2),pigment.v(:,3),'Facecolor','w','FaceAlpha',0.1,'EdgeAlpha',0.1)

% surf_tree = Surf_tree(cornia_surface, lens_surface, retina_surface);
% surf_tree.surf_refract_index = [1.35 1.4 1.37]; % Refractive index of the surfaces in the order they where entered.
% surf_tree.surface_blocking = [0 0 1]; % Last surface (retina in this case) is set to be a blocking surface.
% surf_tree.refract_order = [1 2 2 3]; % In this case the known order of refraction is (Cornia -> Lens -> Lens again -> Retina).
%------------------------------------------------------------------------

