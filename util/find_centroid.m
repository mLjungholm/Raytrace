function centerOfMass = find_centroid(grayImage, cutoff ,plot_flag)

cutoff_val = max(max(grayImage))*cutoff;
grayImage(isnan(grayImage)) = 0;
grayImage(grayImage < cutoff_val) = 0;
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
    plot(centerOfMass(1), centerOfMass(2), 'r+', 'LineWidth', 1, 'MarkerSize', 20);
    colormap(inferno)
end
end