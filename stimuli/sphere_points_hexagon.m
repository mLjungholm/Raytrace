% Octagon_Points creates a set of points on a sphere where the points are
% spaced att the center of hexagons fitted on the surface of the sphere.

% It is not mathematicly possible to fit only hexagons on a sphere [Euler's
% theorem]. It is however posible to fit hexagons onto the faces of a
% Platonic solid and project the hexagons onto a sphere. This will create a
% pentagon in each vertex (corner) of the triangles in the platonic solid.
% The solid that is the most spherical is the icosahedron (This will cause
% the least amount of warping of the hexagons) will requiere a total of 12
% pentagons in the hexagon sphere. 

% Since it is only possible to fit a certain number of hexagons evenly in
% each triangle; the angular spacing (Dihedral angle) is limited to a set series of numbers.
% The program will try to achieve as a close fit as possible to the input
% desirered cycles/deg. Note that verry low spatial frequincies will result
% in a high step between each allowed angle. 

% v1.0
% Written by: Mikael Ljungholm,  2019-01-17

function [point_list, angle_out] = sphere_points_hexagon(separation_deg, half_sphere, plot, type)

% Find the closet poisble number of hexagon divisions. 
dh_angle_in = separation_deg;
divs = round(63.4349/separation_deg)-1;
angle_out = 63.4349 /(divs + 1);

% dh_angle = 138.19 / (divs + 1); % Dihedral angle /  angle between nodal points.

% Define the points on the triangle.
tri_points = create_triangle(divs);
% Define the points and triangle vectors of the icosahedron
[ico_points, ico_triangs] = create_ico();

% Map all the points for each triangle on to the sphere.
point_list = zeros(size(tri_points,1)*20,3);
p_count = 1;
% for tri_ind = 1:20
for tri_ind = 1:20
    inds = ico_triangs(tri_ind,:);
    v1 = ico_points(inds(1),:);
    v2 = ico_points(inds(2),:);
    v3 = ico_points(inds(3),:);
    for p_ind = 1:size(tri_points,1)
        p = tri_points(p_ind,:);
        p3d = mapGridpoint2Sphere(p,v1,v2,v3);
        point_list(p_count,:) = p3d;
        p_count = p_count + 1;
    end
end
% point_list = sortrows(point_list);
point_list = uniquetol(point_list,10^(-4),'ByRows',true);

if half_sphere
    inds = point_list(:,1) >= 0;
    point_list = point_list(inds,:);
end

if plot
    figure(1)
    hold on
    axis equal
    scatter3(point_list(:,1),point_list(:,2),point_list(:,3),'.')
end

point_n = 10;
if point_n < size(point_list,1)
    point_n = size(point_list,1);
end

% This code used to test if the asigned angles where correct.
% min_angles = zeros(10,10);
% for ind1 = 1:10
%     temp_list = zeros(10,1);
%     for ind2 = 1:10
%         temp_list(ind2) = mydot(point_list(ind1+1,:),point_list(ind2+1,:));
%     end
%     min_angles(ind1,:) = sort(temp_list);
% end
% angle_out = median(sort(min_angles(:,2)));

switch type
    case 'cartesian'
        
    case 'spherical'
        new_list = zeros(size(point_list,1),2);
        for ind = 1:size(point_list,1)
            [az,el,r] = cart2sph(point_list(ind,1),point_list(ind,2),point_list(ind,3));
            new_list(ind,:) = [az,el];
        end
        point_list = new_list;
    otherwise
end

    % Function for creating the hexagon points on a triangle to be mapped
    % on the sphere. The riangle is deffined by [-1/2,0; 1/2,0;
    % 0,1/2*sqrt(3)].
    function tri_points = create_triangle(divs)
        base_triangle = [-1/2,0; 1/2,0; 0,1/2*sqrt(3)]; % Base triangle.
        nodal_baseP = 2 + divs; % Number of nodal basepoints.
        nodal_N = (divs+2)*(divs+3)/2; % Number of nodal points for one triangle surface. (Triangular numbers)

        % The vertex code can be used if the borders of the hexagons are of
        % interest and not only the centers.
        % vertex_baseP = 1 + divs; % Number of vertex basepoints.
        % vertex_N = (divs+1)*(divs+2)/2; % Number of vertex points in one triangle surface.
        % vertex_inds = reshape([1:(vertex_N*4)],[],4); % List of vertex indices.
        % vertex_points = zeros(vertex_N*4,2); % List of vertex points.

        tri_points = zeros(nodal_N,2);   % List of nodal points.
        dx_spacing = (base_triangle(2,1) - base_triangle(1,1)) / (divs + 1);
        dy_spacing = (base_triangle(3,2) - base_triangle(1,2)) / (divs + 1);
        
        % Assign coordinates to the nodal points.
        h = 0;
        y = base_triangle(1,2); x = base_triangle(1,1);
        n = 1;
        for ind = 1:nodal_baseP
            for j = 1:(nodal_baseP-h)
                tri_points(n,:) = [x,y];
                x = x + dx_spacing;
                n = n + 1;
            end
            h = h + 1;
            x = base_triangle(1,1) + (dx_spacing / 2)*ind;
            y = y + dy_spacing;
        end
    end

    % Creating the coordinates for the icosahedron.
    function [ico_points, ico_triangs] = create_ico()
        s = 2/sqrt(5);
        c = 1/sqrt(5);
        top_points = zeros(5,3);
        for ind = 0:4
            top_points(ind+1,:) = [s*cos(ind*2*pi/5), s*sin(ind*2*pi/5), c];
        end
        bottom_points = top_points.*[-1,1,-1];
        ico_points = [ [0,0,1]; top_points; [0,0,-1]; bottom_points];
        
        ico_triangs = zeros(20,3);
        for ind = 1:5
            ico_triangs(ind,:) = [1,ind+1,mod(ind,5)+2];
            ico_triangs(ind+5,:) = [7,ind+7,mod((ind),5)+8];
            ico_triangs(ind+10,:) = [ind+1,mod((ind),5)+2,mod((8-ind),5)+8];
            ico_triangs(ind+15,:) = [ind+1,mod((8-ind),5)+8,mod((9-ind),5)+8];
        end
    end

    % Converitng ponts on the triangle [-1/2,0; 1/2,0; 0,1/2*sqrt(3)] to
    % barycentric coordinates.
    function [l1,l2,l3] = barycentric_coords(p)
        x = p(1);
        y = p(2);
        l3 = y.*2/sqrt(3); % l3*sqrt(3)/2 = y
        l2 = x + 0.5*(1 - l3); % 0.5*(l2 - l1) = x
        l1 = 1 - l2 - l3; % l1 + l2 + l3 = 1
    end

    % Mapping the trianhgle points onto the sphere
    function out = mapGridpoint2Sphere(p,s1,s2,s3)
        [l1,l2,l3] = barycentric_coords(p);
        if norm(l3-1) < 1e-10
            out = s3;
            return;
        end
        l2s = l2/(l1+l2);
        p12 = my_slerp(s1,s2,l2s);
        out = my_slerp(p12,s3,l3);
    end

    % Spherical linear interpolation
    function map = my_slerp(p0,p1,t)
        ang0_cos = dot(p0,p1)/dot(p0,p0);
        ang0_sin = sqrt(1 - ang0_cos*ang0_cos);
        ang0 = atan2(ang0_sin,ang0_cos);
        l0 = sin((1-t)*ang0);
        l1 = sin(t*ang0);
        map = (l0.*p0 + l1.*p1)./ang0_sin;
    end

function val = mydot(v,u)
    if u(1) == v(1) && u(2) == v(2) && u(3) == v(3)
        val = 0;
    else
        av = sqrt(v(1)^2 + v(2)^2 + v(3)^2);
        au = sqrt(u(1)^2 + u(2)^2 + u(3)^2);
        d1 = v(1)*u(1) + v(2)*u(2) + v(3)*u(3);
        d2 = d1/(av*au);
        val = acosd(d2);
    end
end

end
