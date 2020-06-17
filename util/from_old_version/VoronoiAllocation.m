% Creates voronoi cells for the given points and returns a pointList index
% for each pixel in the image determined by the imSize

% pointList = [k x 2] matrix of points [x,y]
% imSize = [m , n] size of the image

% pointIndex = [m x n] matrix indicating which point each pixel belongs to


function voronoi_map = VoronoiAllocation(pointList,imSize)
voronoi_map = zeros(imSize(1),imSize(2));
for xi = 1:imSize(2)
    for yi = 1:imSize(1)
        closestInd = nan;
        closestDist = inf;
        for pi = 1:size(pointList,1)
            dist = sqrt((pointList(pi,1)-yi)^2 + (pointList(pi,2)-xi)^2);
            if dist < closestDist
                closestDist = dist;
                closestInd = pi;
            end
        end
        voronoi_map(yi,xi) = closestInd;
    end
end

% Old version
% pixelCount = imSize(1)*imSize(2);
% imMat = zeros(pixelCount,2);
% k = 1;
% for i = 1:imSize(1)
%     for ii = 1:imSize(2)
%         imMat(k,:) = [i-0.5, ii-0.5];
%         k = k+1;
%     end
% end
% 
% pointIndex = zeros(pixelCount,1);
% for i = 1:pixelCount
%     closestInd = nan;
%     closestDist = inf;
%     for ii = 1:size(pointList,1)
%         dist = sqrt((imMat(i,1)-pointList(ii,1))^2+(imMat(i,2)-pointList(ii,2))^2);
%         if dist < closestDist
%             closestDist = dist;
%             closestInd = ii;
%         end
%     end
%     pointIndex(i) = closestInd;
% end
% 
% pointIndex = reshape(pointIndex,[imSize(1),imSize(2)]);
end