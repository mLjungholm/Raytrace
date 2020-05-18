% Returns the coordinate for the weighted Centroid of the image. 

function centerOfMass = find_weighted_centroid(grayImage, plot_flag)
binaryImage = true(size(grayImage));
% labeledImage = bwlabel(binaryImage);
measurements = regionprops(logical(binaryImage), grayImage, 'WeightedCentroid');
centerOfMass = measurements.WeightedCentroid;
% imSize = size(grayImage)
% dx = (range(2)-range(1))/size(grayImage,1);
% dy = (range(4)-range(3))/size(grayImage,2);
% center_coords = centerOfMass;
% center_coords(1) = (center_coords(1) + range(1))*dx;
% center_coords(2) = (center_coords(2) + range(3))*dy;


if plot_flag
    figure(1)
    hold on
    imagesc(grayImage)
    plot(centerOfMass(1), centerOfMass(2), 'r+', 'LineWidth', 1, 'MarkerSize', 10);
    colormap(inferno)
end
end