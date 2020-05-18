%% Complimet analysis and vizualisation to the filter creation. 

% This part selects one source point and plots the absorption for that
% source on the back of the retina volume. It is usefull for visualization
% of the absorrption and determining directions.

% scN - source point index
% scN = 3687; % From below(model orientation). Or the forward (real world orientation)
% scN = 1; % Almoast straight
scN = 230; % Positive Y -> from Above (real world)
rcN = 2945; % rcN - receptor index

% Pots the response function of a source point onto the back of the retina.

% [Az, El] = meshgrid(-pi:0.01:pi,-pi:0.01:pi);
% [az,el,r] = cart2sph(receptor_coords(:,1),receptor_coords(:,2),receptor_coords(:,3));
% % interpolate nonuniformly spaced points
% C = griddata(az,el,abs_mat_2(:,scN),Az,El);
% R = griddata(az,el,r,Az,El);
% % convert to cart
% [x, y, z] = sph2cart(Az,El,R);

% Visualization
figure(2)
hold on
% colormap(inferno)
axis equal on vis3d
% surface(x,y,z,C,'edgealpha',0.05)
% plot3(source_coords(scN,1),source_coords(scN,2),source_coords(scN,3),'+')
scatter3(receptor_coords(:,1), receptor_coords(:,2), receptor_coords(:,3),'.')
plot3(receptor_coords(rcN,1),receptor_coords(rcN,2),receptor_coords(rcN,3),'o')

%% This part creates angular response maps for individual photoreceptors.

% rcN = 3939; % Receptor index. (Back of retina)
% rcN = 4041; % central
% rcN = 251; % central down (looking upp)
% rcN = 291; % central upp (looking down)
rcN = 2777; % central back (looking forward)

%Convert the radians of the sources to angular for eayy visualization.
[az,el,r] = cart2sph(-source_coords(:,1),source_coords(:,2),-source_coords(:,3));
source_angular = [az,el].*180/pi;

[aza,ela] = meshgrid(-100:0.2:100, -100:0.2:100); % Create a grid that spanns the whole field.
C = griddata(source_angular(:,1),source_angular(:,2),abs_mat_2(rcN,:),aza,ela); % Interpolate the data points onto the grid.
C = imrotate(C,90); % Rotate the image to compensate for the orientation of the 3d model of the eye.
C = flipud(C); % Flip the matrix since matlab plots images with inverted y axis
C(isnan(C)) = 0; % Remove all NaN values from the matrix. 

% Find the weighted center of the image
I = C./max(max(C)); % Normalize the image.
cutoff = 0.5; % Create a cutoff intensity value ( 50% in this case)
I(I < cutoff) = 0; % Remove values below cutoff.
range = [-100,100,-100,100]; % Define the range of the image (same as the meshgrid)
center_coords = find_weighted_centroid(I,0); % Get the weighted center coordinates (in pixels).
dx = (range(2)-range(1))/size(C,1); % Find the deg/pixel scale.
dy = (range(4)-range(3))/size(C,2);
center_coords(1) = (center_coords(1) + range(1))*dx; % Scale the coordinates
center_coords(2) = (center_coords(2) + range(3))*dy;

C = C./max(max(C));
C(250,500:650) = 0.7;
% Visualization.
figure(2)
colormap(inferno)
hold on
imagesc([-100 100],[-100 100],C)
% plot([500 150],[750 750],'w','LineWidth',8)
% plot(center_coords(1), center_coords(2), 'r+', 'LineWidth', 1, 'MarkerSize', 10);
xlabel('Horizontal angle [deg]')
ylabel('Elevation angle [deg]')
axis equal
xtickformat('degrees')
ytickformat('degrees')


%% Creating angular coordinates for the receptors

[az,el,] = cart2sph(receptor_coords(:,1),receptor_coords(:,2),receptor_coords(:,3));
receptor_angular = [az,el].*180/pi;

%% Creating a list of viewing directions for all the receptors.
% rcN = 3939; % Receptor index. (Back of retina)
% rcN = 4041; % central
% rcN = 251; % central down (looking upp)
% rcN = 291; % central upp (looking down)
% rcN = 2777; % central back (looking forward)
view_dir = zeros(4063,2);
for ind = 1:4063
    view_dir(ind,:) = find_view_direction(abs_mat_2(ind,:),source_angular,-100:1:100,-100:1:100);
end

%     view_dir(3,:) = find_view_direction(abs_mat_2(rcN,:),source_angular,-100:1:100,-100:1:100);

% Converting to carteesian coordinates
view_dir_cart = zeros(4063,3);
for ind = 1:4063
    [x,y,z] = sph2cart(view_dir(ind,1)/180*pi,view_dir(ind,2)/180*pi,100);
    view_dir_cart(ind,:) = [x,y,z];
end
%% Visualizatrion of the receptors, their direction and absorption function
% rcN = 3939; % Receptor index. (Back of retina)
rcN = 4041; % central
% rcN = 251; % central down (looking upp)
% rcN = 291; % central upp (looking down)
% rcN = 2777; % central back (looking forward)
% view_dir = find_view_direction(abs_mat_2(rcN,:),source_angular,-100:1:100,-100:1:100);
% [x,y,z] = sph2cart(view_dir(1)/180*pi,view_dir(2)/180*pi,100);
% view_dir = [x,y,z];


figure(1)
hold on
% scatter3(-view_dir_cart(:,1),view_dir_cart(:,2),view_dir_cart(:,3),'.b')
plot3(receptor_coords(:,1), receptor_coords(:,2), receptor_coords(:,3),'.g','MarkerSize', 1)
plot3(receptor_coords(rcN,1),receptor_coords(rcN,2),receptor_coords(rcN,3),'o','LineWidth', 2, 'MarkerSize', 10)
plot3(-view_dir_cart(rcN,1),view_dir_cart(rcN,2),view_dir_cart(rcN,3),'r*','LineWidth', 2, 'MarkerSize', 10)
% plot3(source_coords(:,1).*1.3,source_coords(:,2).*1.3,source_coords(:,3).*1.3,'.b','MarkerSize', 1)
plot3(-view_dir_cart(:,1),view_dir_cart(:,2),view_dir_cart(:,3),'.b','MarkerSize', 1)
axis equal 
grid on


[Az, El] = meshgrid(-pi:0.01:pi,-pi:0.01:pi);
[az,el,r] = cart2sph(-source_coords(:,1),source_coords(:,2),source_coords(:,3));
% interpolate nonuniformly spaced points
C = griddata(az,el,abs_mat_2(rcN,:)',Az,El);
R = griddata(az,el,r,Az,El);
% convert to cart
[x, y, z] = sph2cart(Az,El,R);

% Visualization
figure(2)
hold on
% colormap(inferno)
axis equal on vis3d
surface(-x,y,z,C,'edgealpha',0.05)
plot3(receptor_coords(:,1), receptor_coords(:,2), receptor_coords(:,3),'.g','MarkerSize', 4)
plot3(receptor_coords(rcN,1),receptor_coords(rcN,2),receptor_coords(rcN,3),'o','LineWidth', 2, 'MarkerSize', 10)
plot3(-view_dir_cart(rcN,1),view_dir_cart(rcN,2),view_dir_cart(rcN,3),'r*','LineWidth', 2, 'MarkerSize', 10)

%%  Testing the fishey projection on a few viewing directions.

% Original image size and new image size
iSizeOr = 4032;
radOr = (3937-48)/4032;
diamMod = radOr;
iSize = 512;

rcN = 4041; % central
% rcN = 251; % central down (looking upp)
% rcN = 291; % central upp (looking down)
% rcN = 2777; % central back (looking forward)

receptor_inds = [251,291,2777];
I = zeros(iSize,iSize);
pointList = zeros(3,2);

for ind = 3:3
     az = -view_dir(receptor_inds(ind),1)/180*pi;
     el = view_dir(receptor_inds(ind),2)/180*pi;
    [xi,yi] = sph2pixel([az,el],[iSize iSize],round(iSize*diamMod));
    I(yi,xi) = 1;
    pointList(ind,:) = [xi,yi];
end

% figure(1)
imagesc(I)
% axis equal
%% map the source coordinates onto the fishey projection.

% Original image size and new image size
iSizeOr = 4032;
radOr = (3937-48)/4032;
iSize = 256;
source_fishey_coords = zeros(3900,2);

for source_ind = 1:3900
    [xi,yi] = sph2pixel([source_coordinates(source_ind,1)/1000,source_coordinates(source_ind,2)/1000],[iSize iSize],round(iSize*radOr));
    source_fishey_coords(source_ind,:) = [xi,yi];
end

% Question: Should the interpolation be done in 3d space or in the 2D fisheye projection? Will it make a diffenence?

%% Creating one absorption image for one receptor
rcN = 4041; % central
% rcN = 251; % central down (looking upp)
% rcN = 291; % central upp (looking down)
% rcN = 2777; % central back (looking forward)

% abs_mat_2D_single = single(zeros(iSize,iSize,6));
iSize = 256;
abs_mat_2D_single = zeros(iSize,iSize);

[yq,xq] = meshgrid(1:iSize, 1:iSize);
interpolated_data = scatteredInterpolant(source_fishey_coords(:,1),source_fishey_coords(:,2),abs_mat_2(rcN,:),yq,xq);

colormap(viridis)
axis equal
imagesc(interpolated_data)

%% Creating voronoi cells
radOr = (3937-48)/4032;
iSize = 256;

view_dir_fisheye_coords = zeros(4063,2);
for ind = 1:4063
    [xi,yi] = sph2pixel([view_dir(ind,1)/180*pi,view_dir(ind,2)/180*pi],[iSize iSize],round(iSize*radOr));
    view_dir_fisheye_coords(ind,:) = [xi,yi];
end
voronoiList = VoronoiAllocation(view_dir_fisheye_coords,[iSize,iSize]);

%% Testing the voronoi allocation by plotting one cell and its coresponding view_direction.
rcN = 4041;

I = voronoiList;
I(I ~= rcN) = 0;
I(I > 1) = 1;
figure(1)
hold on
imagesc(I)


%% Creating an absorptionMatrix for all photoreceptors

% Workaround for the sizelimit of the matrix.
a = single(zeros(iSize,iSize,2000));
b = single(zeros(iSize,iSize,2063));
absMat2D = cat(3,a,b);

% absMat2D = single(zeros(iSize,iSize,receptorNum));
% coordTransferMat = zeros(4063,2);
% iSizeOr = 4032;
% radOr = (3937-48)/4032;

% for i = 1:sourceNum
%     [yi,xi] = sph2pixel([-coords(i,2)/100,coords(i,1)/100],[iSize iSize],round(iSize*diamMod));
%     coordTransferMat(i,:) = [yi,xi];
% end

[yq,xq] = meshgrid(1:iSize, 1:iSize);

h = waitbar(0,'Initializing waitbar...');
perc = 0;
step = 100/4063;
for i = 1:4063
    vq = griddata(source_fishey_coords(:,1),source_fishey_coords(:,2),abs_mat_2(i,:),yq,xq);
    absMat2D(:,:,i) = vq;
    
    perc = perc + step;
    tempPerc = round(perc);
    waitbar(tempPerc/100,h,sprintf('%d%% done...',tempPerc))
end
absMat2D(isnan(absMat2D)) = 0;

%% White field calibration

iSize = 256;

receptor_total_abs = zeros(4063,1);
for ind = 1:4063
%     receptor_total_abs(ind) = sum(absorptionMat(ind,:));
    receptor_total_abs(ind) = sum(sum(absMat2D(:,:,ind)));
end

white_filter = single(zeros(iSize,iSize));

for receptor_ind = 1:4063
    white_filter(voronoiList == receptor_ind) = 1/receptor_total_abs(receptor_ind);
end

%% Edge filter 

edge_filter = single(zeros(iSize,iSize));
for ind = 1:4063
    edge_filter(voronoiList == ind) = receptor_total_abs(ind);
end

%% test

% testIm = single(zeros(iSize,iSize));
% testIm = edge_filter.*white_filter;
% imagesc(testIm)
im2 = uint8(filterIm.*255);


%% Loading one Image

iSize = 256;
[filename, path] = uigetfile('*');
filePath = strcat(path,filename);
I = imread(filePath);
I = I(:,:,2);
% I = rgb2gray(I);   % Transform into gray-scale.
% 
B=I;
B = B(35:4792,:);
B = B(:,745:5503);
B = imresize(B,[iSize iSize]);
C = im2single(B);


figure(2)
imshow(C)

%% Running the filter on one image
% [filename, path] = uigetfile('*');
% filePath = strcat(path,filename);
% I = imread(filePath);
% 
% B=I;
% B = ones(256);
% C = im2single(B);
% figure(2)
% imshow(C+BW)
% C = flip(C,2);
C = flip(C,1);

filterIm = single(zeros(iSize,iSize));
for ind = 1:4063
%     if list(i) < 85
    absVal = sum(sum(absMat2D(:,:,ind).*C));
    filterIm(voronoiList == ind) = absVal;
%     else
%         filterIm(voronoiList == i) = 0;
%     end
end
% filterIm = filterIm.*inCalibIm.*lightFilter;
filterIm = filterIm.*white_filter;
filterIm = filterIm./max(max(filterIm)).*255;
% % filterIm = flip(filterIm,2);
fI2 = uint8(filterIm)-edge_filter;


comIm = C(:,:,[1 1 1]).*255;
comIm = flip(comIm,1);
comIm(:,:,1) = comIm(:,:,1) + BW;
comIm = cat(2,uint8(comIm),fI2(:,:,[1 1 1]));
comIm = insertMarker(comIm,[192,124],'plus','color','red','size',5);
% 
% comIm = insertMarker(C.*255,[192,124],'plus','color','red','size',5);
% comIm(:,:,1) = comIm(:,:,1) + BW;
% comIm = cat(2,uint8(comIm),fI2(:,:,[1 1 1]));
% t1(:,:,1) = t1(:,:,1) + BW;
figure(1)
imshow(comIm)

% filename2 = strcat('D:\RayTrace_Program\beta0.3\veletwormFilter\newSim2\filtered_images\','filtered_',filename(1:end-10),'_',filename(end-9:end-5),'cycles.tiff');
filename2 = strcat('D:\RayTrace_Program\beta0.3\veletwormFilter\newSim2\filtered_images\','filtered_',filename,'.tiff');
imwrite(comIm,filename2)
% imshow(filterIm2)
%
% figure(3)
% hist(reshape(filterIm2,1,[]))
% ! The images needs to be flipped because I made an error with the
% coordinates of the phooreceptors vs source points. !
% filterIm = flip(filterIm);
% filterIm = fliplr(filterIm);

% figure(1)
% imshow(filterIm)

%% Some image modifications

BW = edge(edge_filter);
BW = single(BW).*255;

radOr = (3937-48)/4032;
iSize = 256;

phi = 0;
theta = -50/180*pi;
% [x,y] = pickFishEye(phi,theta, round(iSize*radOr));
[x,y] = pickFishEye(0,-0.8727, round(246.9206));

t1 = im2single(B);
t1 = insertMarker(t1,[x,y],'plus','color','red','size',5);
t1(:,:,1) = t1(:,:,1) + BW;
% t1 = t1 + crossim;

figure(2)
imshow(t1)

%% Control tests
test_abs = zeros(1,4063); 
for i = 1:4063
     test_abs(i) = sum(sum(absMat2D(:,:,i)));
end

%%
[x,y] = pickFishEye(0,-0.8727, round(246.9206));
