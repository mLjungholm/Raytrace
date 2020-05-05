% script for creating the points used in the receptor_tree

% Note! 
%   Function for creating cones is missing

% Step 1. Create the receptor shapes -> make lists of variables defining
% the space

% Step 2. Create a grid space with the desiered resolution.

% Step 3. Match all grid points with the receptor volumes and make list of
% the receptors that match each point.

% Step 4. Create an octree for the receptor points.

% Step 5. Trace the rays though the octree and add absorbed values to the
% respective receptors.


%% Generate gridpoints 
minP = min([mP;oP],[],1);
maxP = max([mP;oP],[],1);
[x,y,z] = meshgrid(minP(1):1:maxP(1),minP(2):1:maxP(2),minP(3):1:maxP(3));
points = [reshape(x,[],1),reshape(y,[],1),reshape(z,[],1)];
points = single(points);


%% Mach points with receptor volumes

cone_dir = single(zeros(size(mP,1),3));
aP = single(zeros(size(mP,1),3));
h = single(zeros(size(mP,1),1));
alpha = single(zeros(size(mP,1),1));
beta = single(zeros(size(mP,1),1));
for i = 1:size(mP,1)
    
    v = oP(i,:)-mP(i,:);
    h(i) = norm(v);
    cone_dir(i,:) = v./h(i);
    aP(i,:) = (mP(i,:)-2*v);
    alpha(i) = atand(4/(3*h(i)));
    beta(i) = atand(0.7/(3*h(i)));
end

hWaitBar = waitbar(0, 'Fitting retina Volume', 'CreateCancelBtn', ...
                   @(src, event) setappdata(gcbf(), 'Cancelled', true));
setappdata(hWaitBar, 'Cancelled', false);

cones_per_points = create_receptor_grid(points,aP,cone_dir,h,alpha,beta,hWaitBar);
% cpp = cones_per_points';

not_empty_points = any(cones_per_points,2); % Map of unempty points
max_fit = max(sum(cones_per_points~=0,2));  % Max cones for any point
fitted_points = points(not_empty_points,:); % New list of unempty points
cones_in_points = cones_per_points(not_empty_points,1:max_fit);  % Cone list / point
% cones_in_points = cones_in_points(:,1:max_fit); % remove unessesary zeros
% 

function cones_per_point = create_receptor_grid(points,cone_ap, cone_dir,h, alpha, beta,hWaitBar)
coneN = size(cone_ap,1);
pointN = size(points,1);

numComplete = 0;

% Initiate fitting
pcount = 1:pointN;
pcount = pcount';
cones_per_point = cell2mat(arrayfun(@point_cone_fit, pcount,'UniformOutput',false));

delete(hWaitBar);


% Refom matrix and remove all points that have zero matches.

    function match = point_cone_fit(p_ind)
        match = single(zeros(1,100)); % Guessing no more than hundred cones for a point
        n = 1;
        for i = 1:coneN
            l = sqrt((points(p_ind,1)-cone_ap(i,1))^2 + (points(p_ind,2)-cone_ap(i,2))^2 + (points(p_ind,3)-cone_ap(i,3))^2);
            if  l > h(i)*2 && l < h(i)*3
                pvec = (points(p_ind,:)-cone_ap(i,:))./norm(points(p_ind,:)-cone_ap(i,:));
                ang = acosd(pvec(1)*cone_dir(i,1) + pvec(2)*cone_dir(i,2) + pvec(3)*cone_dir(i,3));
%                 if  ang <= beta(i) && ang >= alpha(i)
                if  ang <= alpha(i)
                    match(n) = i;
                    n = n + 1;
                end
            end
        end
        numComplete = numComplete + 1;
        fractionComplete = numComplete / pointN;
        waitbar(fractionComplete, hWaitBar);
    end
end


% oc_vol

% ptest = [reshape(x,[],1) reshape(y,[],1) reshape(z,[],1)];
% figure(1)
% scatter3(points(:,1),points(:,2),points(:,3))
% axis equal