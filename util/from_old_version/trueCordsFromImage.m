function [az,el,r] = trueCordsFromImage(pointImg, imSize, diam)

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
    az = 0; el = 0; r = 0;
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

% center = [2016 2016]; % Remember that image matrix is [row, col] ie [y,x]
% % pointImg = [zim,xim]; % Needs to be an actual point
% 
% % centerCord = pointImg - center;
% % rImg = norm(centerCord); % Rewrite this as math if it is to slow.
% % rImg = sqrt(centerCord(1)^2 + centerCord(2)^2);
% rImg = sqrt((pointImg(1)-center(1))^2 + (pointImg(2)-center(2))^2);
% if rImg > 1953
%     az = 0; el = 0; r = 0;
%     return
% end
% 
% pind = [4.55864436973203e-06 0.0370942598652501 0.159856449173671];
% ang = pind(1)*rImg^2 +pind(2)*rImg + pind(3); % This is the angle between cart-y axis and trueV vector.
% 
% % Find the angle alpha if the trueV would lie in the tueX-trueY plane
% % alpha = 90-ang; % Not needed. 
% 
% % Find the angle beta that the trueV vector should be rotated around trueY
% % axis.
% yRot = -acos((pointImg(2)-2016)/rImg);
% if pointImg(1) > 2016
%     yRot = -yRot;
% elseif pointImg(1) == 2016
%     yRot = 0;
% end
% 
% % Find cart (x,y,z) coords before rotation. zi is assumed to be 0 before
% % rotation.
% xi = sind(ang); yi = cosd(ang);
% 
% % Rotation matrix around y axis;
% % Ry = [cos(yRot) 0 sin(yRot); 0 1 0; -sin(yRot) 0 cos(yRot)];
% % trueCords = Ry*[xi;yi;0];
% 
% xf = cos(yRot)*xi;
% % yf = yi;
% zf = -sin(yRot)*xi;
% 
% % Convert to spherical coordinates;
% % [az,el,r] = cart2sph(trueCords(1),trueCords(2),trueCords(3));
% % az = atan(trueCords(2)); el = acos(trueCords(3)); %r = sqrt(trueCords(1)^2 + trueCords(2)^2 + trueCords(3)^2);
% az = atan(xf/yi); el = asin(zf);
% r = 1;