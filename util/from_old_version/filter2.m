% This is a new try to create a TrueSim filter.

% Process order: Make a map (traceCords) of all sph-cords used to create the
% simmulation. Make a map (receptCords) for the coordinates of all the
% receptors. 

% Calibrate the pixel2sph func for relevant immage size. Make a 3d matix
% for the response of every pixel point. 

% Run an immage though the filter.


%% Making the traceCords and receptCords map

traceCords = round(sphcord(:,1:2).*1000);
traceCords = sortrows(traceCords,[2,1]);

receptCords = zeros(4063,3);
for i = 1:4063
    [az,el,r] = cart2sph(coneGrid(i).basep(1),coneGrid(i).basep(2),coneGrid(i).basep(3));
    receptCords(i,:) = [az,el,r];
end

receptCords = sortrows(receptCords,[2,1]);

receptCartCords = zeros(4063,3);
for ii = 1:4063
    [x,y,z] = sph2cart(receptCords(ii,1),receptCords(ii,2),receptCords(ii,3));
    receptCartCords(ii,:) = [x,y,z];
end

%% Calibrate the pixel2sph func

% [filename, path] = uigetfile('*');
% filePath = strcat(path,filename);
% I = imread(filePath);
% I = rgb2gray(I);   % Transform into gray-scale.

B=I;

% B = B(35:4792,:);
% B = B(:,745:5503);
B = imresize(B,[90 90]);
% B = imrotate(B,-90);

imshow(B)
axis equal
