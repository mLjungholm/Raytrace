classdef Absorption_results < handle
    
    properties
        receptor_nums
        source_nums
        source_coords
        absorption_mat
        total_abs
        file_ids
        cF_horizontal
        cF_vertical
        view_dir
    end
    
    methods
        function this = Absorption_results()
        end
        
        % This function loads the traced files and loads the data into the
        % corresponding variable. This function does not do any
        % calulations.
        function create_absorption_mat(this)
            [filename, path] = uigetfile('*','MultiSelect','on');
            fileNums = size(filename,2);
            
            temp_var = load(strcat(path, char(filename(1))),'s');
            this.receptor_nums = size(temp_var.s.absorption,1);
            this.absorption_mat = zeros(fileNums,this.receptor_nums);
            this.source_coords = zeros(fileNums,3);
            this.file_ids = zeros(fileNums,1);

            h = waitbar(0,'Initializing waitbar...');
            perc = 0;
            step = 100/fileNums;
            
            for file_ind = 1:fileNums
                loadstr = load(strcat(path, char(filename(file_ind))),'s');
                source_id = loadstr.s.source_id;
                this.source_coords(source_id,:) = loadstr.s.origin;
                this.absorption_mat(source_id,:) = loadstr.s.absorption';
                this.file_ids(source_id) = file_ind;
               
                perc = perc + step;
                tempPerc = round(perc);
                waitbar(tempPerc/100,h,sprintf('%d%% done...',tempPerc))
            end
            delete(h)
        end
        
        % This function calculates the coutoff frequency and viewing
        % direction for each receptor
        function calculate_cutoff_frequency(this)
            h = waitbar(0,'Calculating cutoff frequencies...');
            perc = 0;
            step = 100/fileNums;
            
            [Az, El] = meshgrid(-90:1:90,-90:1:90);
            [az,el,~] = cart2sph(-this.source_coords(:,1),this.source_coords(:,2),this.source_coords(:,3));
            az = az.*180/pi;
            el = el.*180/pi;
            ind2angle = @(matInd,matSize) (matInd/matSize)*180-90;
            this.cF_horizontal = zeros(this.receptor_nums,1);
            this.cF_vertical = zeros(this.receptor_nums,1);
            this.view_dir = zeros(this.receptor_nums,2);            
            for r_ind = 1:this.receptor_nums
                % interpolate nonuniformly spaced points
                C = griddata(az,el,this.absorption_mat(:,r_ind),Az,El);
                C(isnan(C)) = 0;
                cM = find_centroid(C,0.5,0);
                this.view_dir(r_ind,:) = [ind2angle(cM(1),181),ind2angle(cM(2),181)];
                this.cF_horizontal(r_ind) = CalculateCutoffFrequency(C,'h',0);
                this.cF_vertical(r_ind) = CalculateCutoffFrequency(C,'v',0);
                perc = perc + step;
                tempPerc = round(perc);
                waitbar(tempPerc/100,h,sprintf('%d%% done...',tempPerc))
            end
            delete(h)
        end
        
        function calculate_total_abs(this) 
            this.total_abs = zeros(this.source_nums,1);
            this.total_abs = sum(this.absorption_mat,2);
            this.total_abs = this.total_abs./max(this.total_abs);
        end
        
        % Function for visiualization of the cut_off frequency in 3d
        function visualize_cutoff(this,direction)
            arguments
                this
                direction {mustBeMember(direction,{'v','h'})}
            end
            if isempty(this.absorption_mat)
                error('absorption matrix is emtpy')
            end
            [Az, El] = meshgrid(-pi:0.01:pi,-pi/2:0.01:pi/2);
            r = ones(this.receptor_nums,1)*120;
            R = griddata(this.view_dir(:,1)./180*pi,this.view_dir(:,2)./180*pi,r,Az,El);
            if isequal(direction,'h')
                C = griddata(this.view_dir(:,1)./180*pi,this.view_dir(:,2)./180*pi,this.cF_horizontal,Az,El);
            elseif isequal(direction, 'v')
                C = griddata(this.view_dir(:,1)./180*pi,this.view_dir(:,2)./180*pi,this.cF_vertical,Az,El);
            else
                error('unsupported direction')
            end
            v = 0.015:0.01:0.15;
            C2 = getContour(C,v);
            [x, y, z] = sph2cart(Az,El,R);
            figure(1)
            hold on
            % colormap(viridis)
            % colormap(inferno)
            axis equal off vis3d
            surface(-x,y,z,C2,'edgealpha',0.05)
            function outMat = getContour(matIn, scale_in)
                [yi,xi] = size(matIn);
                outMat = zeros(yi,xi);
                for y_ind = 1:yi
                    for x_ind = 1:xi
                        if ~isnan(matIn(y_ind,x_ind))
                            t = scale_in - matIn(y_ind,x_ind);
                            [~,I] = min(abs(t));
                            outMat(y_ind,x_ind) = scale_in(I);
                        end
                    end
                end
            end
        end
        
        function visualize_total_abs(this,threshold)
            arguments
                this
                threshold {mustBeNumeric, mustBeGreaterThanOrEqual(threshold,0), mustBeLessThan(threshold,1)}          
            end
            [Az, El] = meshgrid(-pi/2:0.01:pi/2,-pi/2:0.01:pi/2);
            [az,el,r] = cart2sph(-this.source_coords(:,1),this.source_coords(:,2),this.source_coords(:,3));
            R = griddata(az,el,r,Az,El);
            C = griddata(az,el,this.total_abs,Az,El);
            C(C<threshold) = nan;
            [x, y, z] = sph2cart(Az,El,R);
            figure(1)
            hold on
            % colormap(viridis)
            colormap(inferno)
            axis equal off vis3d
            surface(-x,y,z,C,'edgealpha',0.05)
        end
        
        % This is a supportfunction to create a 3d contour plot.
        
    end
end
