% This is a flow of how to set up one or more pointsources and testrun
% single ray simulations to do a initial check of the model.

%-------------------------------------------------------------------------
% Step 1. Load the surf_tree object (This shloud have been saved when
% creating the surface_tree).

%-------------------------------------------------------------------------
% Step 2. Visiualize the surfaces to determine a good position for a light
% source.

%-------------------------------------------------------------------------
% Step 3. Testrun a single point source simulation and trace the result to
% verify the model and position of the source. Repeat this for aditional
% positions and directions as is needed for the simulation. 

%                               ! Note !
%       Do not use more than ~1000 rays for the test simulation as this
%       will slow down the graphical rendering.

%-------------------------------------------------------------------------
% Step 4. Create a list of all pointsoruces to be included in the
% simulation. Use the vizualization functions to enshure that the
% positioning looks correct.

%-------------------------------------------------------------------------
% Step 5. Save the point source coordinate list to be used in the main
% raytracing run.

%-------------------------------------------------------------------------
%   Created by Mikael Ljungholm
%   2020-03-02

%   Dependencies:

%-------------------------------------------------------------------------


%% Step 2 Visualize surfaces

%-------------------------------------------------------------------------
% Example
% cornia_surface.plot(1,'r')
% lens_surface.plot(1,'b')
% retina_surface.plot(1,'g')
%-------------------------------------------------------------------------

%% Step 3 

%-------------------------------------------------------------------------
% Example
dir = [0 0 1]; % Source direction in positive z-axis.
origin = [0 0 -20]; % Origin of source.
% hsize = Radius of the source. make it big enough to encompass the entire 
% aperature but not unessesarily big which will result in unessesary 
% computations
hsize = 17;
raysPerAxis = 15; % Density of rays.
refract_index = 1; % Refractive index of the initial medium for the source.
test_source = Source(dir, origin, hsize, raysPerAxis, refract_index);

cornia_surface.plot(1,'r')
lens_surface.plot(1,'b')
retina_surface.plot(1,'g')
test_source.plot_grid(1) % Plots the initial points for all rays.

%% Step 4
