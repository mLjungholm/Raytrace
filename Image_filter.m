classdef Image_filter < handle
    
    properties
        filter_mat
        image_size
        voronoi_map
        edge_filter
        luminosity_filter
        coord_transfer_mat
        source_nums
        receptor_nums
        weights
        down_supersampling
        receptor_view_center
        diamMod = 1.3423e+03/1080;
%         diamMod = (3937-48)/4032;
    end

%                 iSizeOr = 4032;
%             radOr = (3937-48)/4032;
              
    
    methods
        function this = Image_filter()
        end
        
        % This function creates the image filter matrix used to create the
        % filtered images. It takes three input variables. The
        % absorption_result object, the output imageSize. The last input
        % variable is the optics calibration file. 
        
        % The reference frame shuld specify what the reference axis are. In
        % order ["forward" "right" "up"] as clumn vectors. An example from the onycohperan
        % would be:
        %          f    r   u
        %        [-1    0   0 ]
        %        [ 0    1   0 ]
        %        [ 0    0   1 ]
        % Ie the eye is oriented so that the light travels in the +x
        % direction and the eye "looks" in the -x direction.
        function create_filter(this,absorption_results, imageSize)
%             function create_filter(this,absorption_results, optics_calibration_file, imageSize)
            this.source_nums = absorption_results.source_nums;
            this.receptor_nums = absorption_results.receptor_nums;
            this.image_size = imageSize;            
            % List of pixel coordinates for the center of each receptor.
            this.receptor_view_center = zeros(absorption_results.receptor_nums,2);
            
            
            
            for receptor_ind = 1:absorption_results.receptor_nums
%                 receptor_view_center(receptor_ind,:) = sph2pixel(absorption_results.view_dir(receptor_ind,:)./180*pi,[imageSize imageSize],optics_calibration_file);
                [yi,xi] = sph2pixel(absorption_results.view_dir(receptor_ind,:)./180*pi,[imageSize imageSize],this.diamMod*imageSize);
                this.receptor_view_center(receptor_ind,:) = [yi,xi];
            end
            % Create a voronoi cell map using the center pixel coordinates
            % for the receptors
            this.voronoi_map = VoronoiAllocation(this.receptor_view_center,[imageSize,imageSize]);
            % Create a map of where all source coodrinates land in the
            % image plane.
            h = waitbar(0,'Building source to pixel coordinate transfer map...');
            perc = 0;
            step = 100/this.source_nums;
            this.coord_transfer_mat = zeros(this.source_nums,2);
            for source_ind = 1:this.source_nums
                [az,el] = cart2sph(-absorption_results.source_coords(source_ind,1),absorption_results.source_coords(source_ind,2),absorption_results.source_coords(source_ind,3));
                [yi,xi] = sph2pixel([az,el],[imageSize imageSize],round(imageSize*this.diamMod));
                this.coord_transfer_mat(source_ind,:) = [yi,xi];
                perc = perc + step;
                tempPerc = round(perc);
                waitbar(tempPerc/100,h,sprintf('%d%% done...',tempPerc))
            end
            delete(h)
            % Downsampling if there are multiple sourcepoints pointing
            % towards the same pixel. For ease, use a simple multisampling
            % aproach. For maximum quality, (unessesary in my opinion) use a
            % sinc function (Lanczos).
            
            h = waitbar(0,'Checking for source to pixel coordinate overlap...');
            perc = 0;            
            [~,uniq_coords,duplicates] = unique(this.coord_transfer_mat,'rows');
            logic_dups = true(size(duplicates,1),1);
            logic_dups(uniq_coords) = false;
            I = find(logic_dups);
            if any(I)
                step = 100/size(I,1);
                this.weights = zeros(this.source_nums,1);
                for dup_ind = 1:size(I,1)
                    tind = duplicates(I(dup_ind));
                    if this.weights(tind)
                        continue
                    end
                    tindS = duplicates == tind;
                    this.weights(tindS) = 1/sum(tindS);
                    
                    perc = perc + step;
                    tempPerc = round(perc);
                    waitbar(tempPerc/100,h,sprintf('%d%% done...',tempPerc))
                end
            end
            this.down_supersampling = any(this.weights);
            if this.down_supersampling
                this.weights(this.weights == 0) = 1;
            end
            delete(h)
            
%             Creating the receptor filter stack. Each receptor creates one
%             image-filter using the coord_transfer_map and weights vector
%             and a linear interpolation for the empy pixels.
%             We also reduce to single pecision for reduced memmory usage
            h = waitbar(0,'Building receptor filter stack...');
            perc = 0;
            step = 100/this.receptor_nums;
            this.filter_mat = single(zeros(imageSize,imageSize,this.receptor_nums));
            [y_grid,x_grid] = meshgrid(1:imageSize, 1:imageSize);
            for receptor_ind = 1:this.receptor_nums
                if this.down_supersampling
                    interpolated_image = griddata(this.coord_transfer_mat(:,1),this.coord_transfer_mat(:,2),absorption_results.absorption_mat(:,receptor_ind).*this.weights,x_grid,y_grid);
                else
                    interpolated_image = griddata(this.coord_transfer_mat(:,1),this.coord_transfer_mat(:,2),absorption_results.absorption_mat(:,receptor_ind),x_grid,y_grid);
                end
                this.filter_mat(:,:,receptor_ind) = interpolated_image;
                perc = perc + step;
                tempPerc = round(perc);
                waitbar(tempPerc/100,h,sprintf('%d%% done...',tempPerc))
            end
            this.filter_mat(isnan(this.filter_mat)) = 0;
            delete(h)
            
            % Creating the luminosity compensation filter.
            this.luminosity_filter = zeros(imageSize,imageSize);
            for receptor_ind = 1:this.receptor_nums
                abs_val = sum(sum(this.filter_mat(:,:,receptor_ind)));
                this.luminosity_filter(this.voronoi_map == receptor_ind) = abs_val;                
            end
            this.luminosity_filter = this.luminosity_filter./max(max(this.luminosity_filter));
            this.luminosity_filter(this.luminosity_filter == 0) = nan;
            this.luminosity_filter = 1./this.luminosity_filter;
            this.luminosity_filter(this.luminosity_filter == inf) = 0;
            
            % Creating edge filter
            theshold_val = 0.1;
            this.edge_filter = griddata(this.coord_transfer_mat(:,1),this.coord_transfer_mat(:,2),absorption_results.total_abs,y_grid,x_grid);
            this.edge_filter = this.edge_filter./max(max(this.edge_filter));
            this.edge_filter(isnan(this.edge_filter)) = 0;
            this.edge_filter(this.edge_filter > theshold_val) = 1;        
        end
        
        function filter_image_from_file(this)
            [filename, path] = uigetfile('*');
            filePath = strcat(path,filename);
            I = imread(filePath);
            I = rgb2gray(I);   % Transform into gray-scale.
            % %
            B=I;
            B = B(35:4792,:);
            B = B(:,745:5503);
            C = imresize(B,[this.image_size, this.image_size]);
            C = im2single(C);
            filtered_I = zeros(this.image_size, this.image_size);
            for receptor_ind = 1:this.receptor_nums
                tempIm = C.*this.filter_mat(:,:,receptor_ind);
                imVal = sum(sum(tempIm));
                filtered_I(this.voronoi_map == receptor_ind) = imVal;                
            end
            filtered_I = filtered_I./(max(max(filtered_I)))*255;
            filtered_I = uint8(filtered_I.*this.luminosity_filter.*this.edge_filter);
            comp_im = [uint8(C), filtered_I];
            figure(1)
            imshow(comp_im)
            % imwrite(im3,strcat(filePath(1:end-4),'_Filtered.tif'));            
        end
        
        function filter_image(this, input_image)
            I = rgb2gray(input_image);   % Transform into gray-scale.
            B=I;
            B = B(35:4792,:);
            B = B(:,745:5503);
            C = imresize(B,[this.image_size, this.image_size]);
            CI = im2single(C);
            filtered_I = zeros(this.image_size, this.image_size);
            for receptor_ind = 1:this.receptor_nums
                tempIm = CI.*this.filter_mat(:,:,receptor_ind);
                imVal = sum(sum(tempIm));
                filtered_I(this.voronoi_map == receptor_ind) = imVal;                
            end
            % .*this.luminosity_filter.*this.edge_filter
            filtered_I = filtered_I./(max(max(filtered_I)))*255;
            filtered_I = uint8(filtered_I.*this.edge_filter.*this.luminosity_filter);
            comp_im = [C, filtered_I];
            figure(1)
            imshow(comp_im)
            % imwrite(im3,strcat(filePath(1:end-4),'_Filtered.tif'));            
        end
        
        function filter_stimmuli(this)
            [filename, path] = uigetfile('*');
%             selpath = uigetdir('C:\Users\Mikael\Dev\Ray_tracing','Select save folder')
            filePath = strcat(path,filename);
            I = imread(filePath);
            C = im2single(I);
            filtered_I = zeros(this.image_size, this.image_size);
            for receptor_ind = 1:this.receptor_nums
                tempIm = C.*this.filter_mat(:,:,receptor_ind);
                imVal = sum(sum(tempIm));
                filtered_I(this.voronoi_map == receptor_ind) = imVal;                
            end
            filtered_I = filtered_I./(max(max(filtered_I)))*255;
            filtered_I = uint8(filtered_I.*this.luminosity_filter.*this.edge_filter); %
%             comp_im = [uint8(I), filtered_I];
%             figure(1)
%             imshow(filtered_I)
            imwrite(filtered_I,strcat(filePath(1:end-4),'_Filtered.tif'));
%             imshow(comp_im)
%             imwrite(comp_im,strcat(filePath(1:end-4),'_Filtered.tif'));            
        end
        
        function show_single_filter(this, receptor_num)
            I1 = this.filter_mat(:,:,receptor_num);
            I1 = I1./(max(max(I1)));
            I2 = single(zeros(this.image_size));
            I2(this.voronoi_map == receptor_num) = 1;
            comp_im = [I1, I2];
            figure(1)
            imagesc(comp_im)
        end
        
%         function filter_video(this,input_video)
%         end
    end
end