function pixel_sph_coords = create_pixel_coord_list(image_size)
% Creating a list of the pixel2sph coordinates.
radOr = (3937-48)/4032;
diam = round(image_size*radOr);
pixel_sph_coords = zeros(image_size^2,4);
pixel_ind = 1;
for x_ind = 1:image_size
    for y_ind = 1:image_size
        [az,el] = pixel2sph([x_ind,y_ind],[image_size,image_size],diam);
        pixel_sph_coords(pixel_ind,:) = [az,el,x_ind,y_ind];
        pixel_ind = pixel_ind + 1;
    end
end
end