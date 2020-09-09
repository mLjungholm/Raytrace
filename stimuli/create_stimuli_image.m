function stimuli = create_stimuli_image(type, point_list, separation_angle, image_size)
% Function for creation of a dotted 3d stimmulus pattern with a fishey
% projection.
pixel_coords = create_pixel_coord_list(image_size);
stimuli = zeros(image_size);
for pixel_ind = 1:image_size^2
    if pixel_coords(pixel_ind,1) == 0 || pixel_coords(pixel_ind,2) == 0
        continue;
    end
    [min_angle,~] = find_min_angle(pixel_coords(pixel_ind,1:2), point_list);
    point_val = stimuli_value(min_angle*180/pi,separation_angle,type);
    stimuli(pixel_coords(pixel_ind,3),pixel_coords(pixel_ind,4)) = point_val;
end

end