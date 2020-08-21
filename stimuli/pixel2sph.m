function [az,el,r] = pixel2sph(pointImg, imSize, diam)

center = zeros(1,2);
if mod(imSize(2),2)
    center(2) = imSize(2)/2 +0.5;
else
    center(2) = imSize(2)/2;
end

if mod(imSize(1),2)
    center(1) = imSize(1)/2 +0.5;
else
    center(1) = imSize(1)/2;
end

rImg = sqrt((pointImg(1)-center(1))^2 + (pointImg(2)-center(2))^2);
% if rImg > diam/2
%     az = 0; el = 0; r = 0;
%     return
% end

pind = [18.0452291332390 72.6996894778249 0.353902671216221];
rImg2 = rImg/diam*2;
ang = pind(1)*rImg2^2 + pind(2)*rImg2 +pind(3);
if ang > 90
%     az = 0; el = 0; r = 0;
     az = nan; el = nan; r = nan;
    return
end

% Find the angle beta that the trueV vector should be rotated around trueY
% axis.
yRot = -acos((pointImg(2)-center(2))/rImg);
if pointImg(1) > center(1)
    yRot = -yRot;
elseif pointImg(1) == center
    if pointImg(2) < center
        yRot = 3.1426;
    else
        yRot = 0;
    end
end

% Find cart (x,y,z) coords before rotation. zi is assumed to be 0 before
% rotation.
xi = sind(ang); yi = cosd(ang);

% Rotation matrix around y axis;
% Ry = [cos(yRot) 0 sin(yRot); 0 1 0; -sin(yRot) 0 cos(yRot)];
% trueCords = Ry*[xi;yi;0];

xf = cos(yRot)*xi;
% yf = yi;
zf = -sin(yRot)*xi;

% Convert to spherical coordinates;
% [az,el,r] = cart2sph(trueCords(1),trueCords(2),trueCords(3));
% az = atan(trueCords(2)); el = acos(trueCords(3)); %r = sqrt(trueCords(1)^2 + trueCords(2)^2 + trueCords(3)^2);
az = atan(xf/yi); el = asin(zf);
r = 1;
end

