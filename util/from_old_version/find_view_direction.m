% This function takes the absorption values for a photoreceptor and
% inerpolates that data onto all alngles in a range. It then returns the
% weighted center of the angle dependent absorption function. 

% absorption_values - list of absorption values for each source index.
% source_coordinates - list of coordinates (in spherical) for each source
% index.
% rangeAz, rangeEl - 

function view_direction = find_view_direction(absorption_values, source_coordinates, rangeAz, rangeEl)

[aza,ela] = meshgrid(rangeAz, rangeEl); % Create a grid that spanns the whole field.
C = griddata(source_coordinates(:,1),source_coordinates(:,2),absorption_values,aza,ela); % Interpolate the data points onto the grid.
% C = imrotate(C,90); % Rotate the image to compensate for the orientation of the 3d model of the eye.
C = flipud(C); % Flip the matrix since matlab plots images with inverted y axis
C(isnan(C)) = 0; % Remove all NaN values from the matrix.
% imagesc(C)
% Find the weighted center of the image
C = C./max(max(C)); % Normalize the image.
cutoff = 0.5; % Create a cutoff intensity value ( 50% in this case)
C(C < cutoff) = 0; % Remove values below cutoff.
view_direction = find_weighted_centroid(C,0); % Get the weighted center coordinates (in pixels).
dy = (rangeAz(end)-rangeAz(1))/size(C,1); % Find the deg/pixel scale.
dx = (rangeEl(end)-rangeEl(1))/size(C,2);
view_direction(1) = (view_direction(1) + rangeAz(1))*dx; % Scale the coordinates
view_direction(2) = (view_direction(2) + rangeEl(1))*dy;
% view_direction = fliplr(view_direction);
end