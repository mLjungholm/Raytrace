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
% ------------------------------------------------------------


%% Step 2. Checking and orienting surfaces

% Example
% ------------------------------------------------------------
% cornia_surface.plot(1,'r')
% ------------------------------------------------------------

%% Step 3. Create the surface partition

% Example
%------------------------------------------------------------------------
% surf_tree = Surf_tree(cornia_surface, lens_surface, retina_surface);
% surf_tree.surf_refract_index = [1.35 1.4 1.37]; % Refractive index of the surfaces in the order they where entered.
% surf_tree.surface_blocking = [0 0 1]; % Last surface (retina in this case) is set to be a blocking surface.
% surf_tree.refract_order = [1 2 2 3]; % In this case the known order of refraction is (Cornia -> Lens -> Lens again -> Retina).
%------------------------------------------------------------------------

