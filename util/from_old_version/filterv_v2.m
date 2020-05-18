% This script reads all abosorbed receptor grids and creates a matrix with
% an absorption value for each source coordinate [m x n] matrix.

% It also takes the names of the receptor grids to create a list of
% coordinates for the sources and receptors.

% The filter_v2_compliment script includes code for further visualization and
% analysis

%%  This part loads all receptor grids. It is verry time consuming. (run only if needed)

receptorNum = 4063;
sourceNum = 3900;
absorptionMat = zeros(sourceNum,receptorNum);

[filename, path] = uigetfile('*','MultiSelect','on');
fileNums = size(filename,2);
source_coordinates = zeros(fileNums,2);

h = waitbar(0,'Initializing waitbar...');
perc = 0;
step = 100/fileNums;

for j = 1:fileNums
    %     Load the receptor grid.
    
    rG = LoadFile(strcat(path, char(filename(j))),'coneGrid');
    for i = 1:receptorNum
        absorptionMat(j,i) = rG(i).absValue;
    end
    temp_name = filename{j};
    temp_coords = temp_name(11:end-4);
    temp_coords  = strsplit(temp_coords ,'Y');
    temp_coords_1  = str2num(temp_coords{1});
    temp_coords_2  = str2num(temp_coords{2});
    source_coordinates(j,:) = [temp_coords_2,temp_coords_1];
    
    perc = perc + step;
    tempPerc = round(perc);
    waitbar(tempPerc/100,h,sprintf('%d%% done...',tempPerc))
end

%% This converts the loaded coordinates into catreesian coordinates.

source_coords = zeros(3900,3);
for i = 1:3900
    sc = source_coordinates./1000;
    [x,y,z] = sph2cart(sc(i,1),sc(i,2),100);
    source_coords(i,:) = [x,y,-z].*(-1);
end

%% Collecting all coordinates for the receptors in one list
% This requiers you to load one receptor grid manualy.

receptor_Nr = 4063; % Number of receptors
receptor_coords = zeros(receptor_Nr,3); % Coordinate list for the receptors

for i = 1:receptor_Nr
    receptor_coords(i,:) = coneGrid(i).basep;
end

%% Visualization for the receptors and sources

figure(1)
hold on
scatter3(source_coords(:,1),source_coords(:,2),source_coords(:,3),'.')
scatter3(receptor_coords(:,1), receptor_coords(:,2), receptor_coords(:,3),'o')
axis equal
grid on