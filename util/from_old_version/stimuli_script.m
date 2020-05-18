% Script for creating stimuli pattern.
% test


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
[point_list, cycles_out] = sphere_points_hexagon(0.01, 1, 0, 'spherical');
%%
stimuli = create_stimuli_image('tophat', point_list, 1/cycles_out, 256).*255;
%%
imshow(stimuli)




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

function [min_angle,point_index] = find_min_angle(point, point_list)
min_angle = inf;
point_index = 0;
az1 = point(1);
el1 = point(2);
for ind = 1:size(point_list,1)
    el2 = point_list(ind,2);
    az2 = point_list(ind,1);
    da = acos(cos(el1)*cos(el2)*(cos(az1)*cos(az2) + sin(az1)*sin(az2)) + sin(el1)*sin(el2));
 %     da = acos( cos(theta_1)*cos(theta_2) + sin(theta_1)*sin(theta_2)*cos(phi_1-phi_2));
    if da < min_angle
        min_angle = da;
        point_index = ind;
    end
end
end

function val = stimuli_value(angle,stimulis_angle,type)
% Determine pixel value [0:1] for a stimulis pattern. type alterniatives
% are 'tophat' or 'cosine'.
switch type
    case 'tophat'
        if angle < stimulis_angle/2
            val = 1;
        else
            val = 0;
        end
        
    case 'cosine'
        val = cosd(angle*90/stimulis_angle);
        
    case 'gaussian'
        disp('gaussian, not implemented. use cosine instead')
        val = nan;
    otherwise
        dip('Unsupported type of stimuli')
        val = nan;
end
end


function stimuli = create_stimuli_image(type, point_list, separation_angle, image_size)
% Function for creation of a dotted 3d stimmulus pattern with a fishey
% projection.
pixel_coords = create_pixel_coord_list(image_size);
stimuli = zeros(image_size);
for pixel_ind = 1:image_size^2
    if pixel_coords(pixel_ind,1) == 0 || pixel_coords(pixel_ind,2) == 0
        continue;
    end
%     closest_point = 0;
    [min_angle,point_index] = find_min_angle(pixel_coords(pixel_ind,1:2), point_list);
    point_val = stimuli_value(min_angle,separation_angle,type);
    stimuli(pixel_coords(pixel_ind,3),pixel_coords(pixel_ind,4)) = point_val;
end

end
