% This function loads all the traced and absourbed files and returns a
% struct cointaining the abosption matrix and coordinates.

% if for some reason the source object in the saved file was not named 's'
% then this needs to be changed in the code.

function absorption_results = create_absorption_matrix()
[filename, path] = uigetfile('*','MultiSelect','on');
fileNums = size(filename,2);

temp_var = load(strcat(path, char(filename(1))),'s');
receptor_nums = size(temp_var.s.absorption,1);
absorption_mat = zeros(fileNums,receptor_nums);
source_coords = zeros(fileNums,3);
file_id = zeros(fileNums,1);
% file_name = cell(fileNums,1);

h = waitbar(0,'Initializing waitbar...');
perc = 0;
step = 100/fileNums;

for file_ind = 1:fileNums
    loadstr = load(strcat(path, char(filename(file_ind))),'s');
    source_id = loadstr.s.source_id;
    source_coords(source_id,:) = loadstr.s.origin;
%     abs_vals = loadstr.s.absorption;
    absorption_mat(source_id,:) = loadstr.s.absorption';
    file_id(source_id) = file_ind;
%     file_name{source_id} = filename(file_ind);
    
    perc = perc + step;
    tempPerc = round(perc);
    waitbar(tempPerc/100,h,sprintf('%d%% done...',tempPerc))
end
delete(h)
absorption_results = struct('receptor_nums',receptor_nums,'source_nums',fileNums,'source_coords',source_coords,'absorption_mat',absorption_mat, 'file_id',file_id);

end