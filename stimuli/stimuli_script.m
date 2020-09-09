% Script for creating stimuli pattern.
% test
image_size = 512;
radOr = (3937-48)/4032;
diam = round(image_size*radOr);

[az,el] = pixel2sph([125,50],[image_size,image_size],diam);
% [az,el].*180/pi;
[x,y,z] = sph2cart(az,el,50);
figure(1)
hold on
axis equal
% plot3([0;50],[0;0],[0;0],'r')
quiver3(0,0,0,50,0,0,'r')
quiver3(0,0,0,0,50,0,'b')
quiver3(0,0,0,0,0,50,'g')
plot3(-x,y,z,'o')
grid on

%%
% figure(1)
a_in = 1/0.11;
[point_list, cycles_out] = sphere_points_hexagon(a_in, 1, 0, 'spherical');
a_out = 1/cycles_out
% a_out_control = acosd(dot(point_list(1,:),point_list(2,:)));


% a_out_control = zeros(size(point_list,1),1);
% for point_ind_1 = 1:size(point_list,1)
%     temp_list = zeros(1,size(point_list,1));
%     for point_ind_2 = 1:size(point_list,1)
%         temp_list(point_ind_2) = acosd(dot(point_list(point_ind_1,:),point_list(point_ind_2,:)));
%     end
%     temp_list = sort(temp_list);
%     a_out_control(point_ind_1) = min(temp_list(2:end));
% end

%%

stimuli1 = create_stimuli_image('tophat', point_list, cycles_out, 256);
filename = strcat('C:\Users\Mikael\Dev\Ray_tracing\data\test_data\stimuli_images\tophat',num2str(a_out,'%.3f'),'.tiff');
imwrite(stimuli1,filename)
figure(1)
imshow(stimuli1)

stimuli2 = create_stimuli_image('cosine', point_list, cycles_out, 256);
figure(2)
imshow(stimuli2)
filename = strcat('C:\Users\Mikael\Dev\Ray_tracing\data\test_data\stimuli_images\cosine',num2str(a_out,'%.3f'),'.tiff');
imwrite(stimuli2,filename)
% figure(2)
% % stimuli = imcomplement(stimuli);
% imshow(stimuli)

%%
% filename = fprintnum2str(k,'%02d')
filename = strcat('stimuli_images\tophat',num2str(a_out,'%.3f'),'.tiff');
imwrite(stimuli,filename)


%%
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


function stimuli = create_stimuli_image(type, point_list, separation_angle, image_size)
% Function for creation of a dotted 3d stimmulus pattern with a fishey
% projection.
pixel_coords = create_pixel_coord_list(image_size);
stimuli = zeros(image_size);
for pixel_ind = 1:image_size^2
%     if pixel_coords(pixel_ind,1) == 0 || pixel_coords(pixel_ind,2) == 0
    if isnan(pixel_coords(pixel_ind,1))
        continue;
    end
    [min_angle,point_index] = find_min_angle(pixel_coords(pixel_ind,1:2), point_list);
    point_val = stimuli_value(min_angle*180/pi,separation_angle,type);
    stimuli(pixel_coords(pixel_ind,3),pixel_coords(pixel_ind,4)) = 1-point_val;
end

end
