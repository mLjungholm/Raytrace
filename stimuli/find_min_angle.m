function [min_angle,point_index] = find_min_angle(point, point_list)
% finds the closest point to reference point and the separation angle.
% Coordinates shoud be in radians and in the format [az,el]

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