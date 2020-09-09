% Script for generating horizontal bar stimuli to be used with image filter.

im_size = 256;      % Size of stimuli image
target_width = 20;  % Width of stimuli (in degrees of visual field at horizon)
type = 'sine';      % Type of stimmuli
horizontal_shift = 0; % Shift of the stimulis centre (in degrees)

arena_dist = 45;    % Distance to arena wall
arena_height = 36;  % Arena wall height
% Set to 1 to compensate for the perspective difference between bottom and 
% top of the arena wall
perspective_compensation = 0;       
projection_file = 'dan_equisolid';  % Name of the projection used.

% Step 1. Find the corresponding spherical coordinates for all pixels
% coirdinates in the output image (stimuli image). This depends on the
% projection type used.

stimuli_image = zeros(im_size);
pixel_coords = create_pixel_coord_list(im_size);
for pixel_ind = 1:im_size^2
    if pixel_coords(pixel_ind,1) == 0 || pixel_coords(pixel_ind,2) == 0
        continue;
    end
    if perspective_compensation
%         pixel_val
    else
        pixel_val = 
    end
%     [min_angle,~] = find_min_angle(pixel_coords(pixel_ind,1:2), point_list);
%     point_val = stimuli_value(min_angle*180/pi,separation_angle,type);
%     stimuli(pixel_coords(pixel_ind,3),pixel_coords(pixel_ind,4)) = point_val;
end

function sin_val = get_sin_val(target_width,pixel_angle, horizontal_shift)
sine_width = (target_width/0.6667)*2; % Get period of the sine function from the FWHM target.
if pixel_angle < (-sine_width*3/4 + horizontal_shift) || pixel_angle > (sine_width*3/4 + horizontal_shift)
    sin_val = 0;
else
    

end

