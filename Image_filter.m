classdef Image_filter < handle
    
    properties
        filter_mat
        image_size
        voronoi_map
        edge_filter
        luminosity_filter
    end
    
    methods
        function this = Image_filter()
        end
        
        % This function creates the image filter matrix used to create the
        % filtered images. It takes three input variables. The
        % absorption_result object, the output imageSize. The last input
        % variable is the optics calibration file. 
        function create_filter(this,absorption_results, optics_calibration_file, imageSize)
            this.image_size = imageSize;
%             iSizeOr = 4032;
%             radOr = (3937-48)/4032;
%             diamMod = 1.3423e+03/1080;
%             iSize = 512;
            
            I = zeros(iSize,iSize);
            pointList = zeros(receptorNum,2);
            for j = 1:receptorNum
                [az,el,r] = cart2sph(rG(j).basep(1),rG(j).basep(2),rG(j).basep(3));
                [yi,xi] = sph2pixel([az,el],[iSize iSize],round(iSize*radOr));
                pointList(j,:) = [xi,yi];
            end
            
            voronoiList = VoronoiAllocation(pointList,[iSize,iSize]);
            
%             ! This code is alternative so to not need to flip the end image
            
            Original image size and new image size
            iSizeOr = 4032;
            radOr = (3937-48)/4032; % I do not realy know what this is? Diameter
            factor or something?
            receptorNum = 4063;
            sourceNum = 3900;
            diamMod = 1.3423e+03/1080;
            iSize = 512;
            
            I = zeros(iSize,iSize);
            pointList = zeros(receptorNum,2);
            for j = 1:receptorNum
                [az,el,r] = cart2sph(rG(j).basep(1),rG(j).basep(2),rG(j).basep(3));
                [yi,xi] = sph2pixel([-az,-el],[iSize iSize],round(iSize*diamMod));
                pointList(j,:) = [xi,yi];
            end
            
            voronoiList = VoronoiAllocation(pointList,[iSize,iSize]);
        end
        
        function filter_image(this,input_image)
        end
        
        function filter_video(this,input_video)
        end
    end
end