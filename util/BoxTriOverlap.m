function overlap = BoxTriOverlap(boxEgdes, triVert)
% If all tests are negative, then the triangle intersects the box.
overlap = 1;

tri2 = sort(triVert);
% xt = sort(triVert(:,1));
% yt = sort(triVert(:,2));
% zt = sort(triVert(:,3));

% if ~BBIntersection()
%     overlap = 0;
%     return;
% end
    if tri2(1,1) > boxEgdes(4) || tri2(3,1) < boxEgdes(1)
        overlap = 0;
        return;
    elseif tri2(1,2) > boxEgdes(5) || tri2(3,2) < boxEgdes(2)
        overlap = 0;
        return;
    elseif tri2(1,3) > boxEgdes(6) || tri2(3,3) < boxEgdes(3)
        overlap = 0;
        return;
    end



[c, hb] = BoxCenter();
v = zeros(3,3);
v(1,:) = SUB(triVert(1,:),c);
v(2,:) = SUB(triVert(2,:),c);
v(3,:) = SUB(triVert(3,:),c);

f = zeros(3,3);
f(1,:) = SUB(triVert(2,:),triVert(1,:));
f(2,:) = SUB(triVert(3,:),triVert(2,:));
f(3,:) = SUB(triVert(1,:),triVert(3,:));

% fex = norm(f(1,1));
% fey = norm(f(1,2));
% fez = norm(f(1,3));
fex = sqrt(f(1,1)^2);
fey = sqrt(f(1,2)^2);
fez = sqrt(f(1,3)^2);
if ~AxTestX01(f(1,3),f(1,2),fez,fey)
    overlap = 0;
    return;
elseif ~AxTestY02(f(1,3),f(1,1),fez,fex)
    overlap = 0;
    return;
elseif ~AxTestZ12(f(1,2),f(1,1),fey,fex)
    overlap = 0;
    return;
end

% fex = norm(f(2,1));
% fey = norm(f(2,2));
% fez = norm(f(2,3));
fex = sqrt(f(2,1)^2);
fey = sqrt(f(2,2)^2);
fez = sqrt(f(2,3)^2);
if ~AxTestX01(f(2,3),f(2,2),fez,fey)
    overlap = 0;
    return;
elseif ~AxTestY02(f(2,3),f(2,1),fez,fex)
    overlap = 0;
    return;
elseif ~AxTestZ0(f(2,2),f(2,1),fey,fex)
    overlap = 0;
    return;
end

% fex = norm(f(3,1));
% fey = norm(f(3,2));
% fez = norm(f(3,3));
fex = sqrt(f(3,1)^2);
fey = sqrt(f(3,2)^2);
fez = sqrt(f(3,3)^2);
if ~AxTestX2(f(3,3),f(3,2),fez,fey)
    overlap = 0;
    return;
elseif ~AxTestY1(f(3,3),f(3,1),fez,fex)
    overlap = 0;
    return;
elseif ~AxTestZ12(f(3,2),f(3,1),fey,fex)
    overlap = 0;
    return;
end




if ~planeBoxIntersect()
    overlap = 0;
    return
end




%%% X-Tests%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function flag = AxTestX01(a,b,fa,fb)
        p0 = a*v(1,2) - b*v(1,3);			       	   
        p2 = a*v(3,2) - b*v(3,3);			       	   
        minp = min(p0,p2);
        maxp = max(p0,p2);
        rad = fa *hb(2) + fb * hb(3);   
        if minp > rad || maxp < -rad
            flag = 0;
            return;
        end
        flag = 1;
    end
    function flag = AxTestX2(a,b,fa,fb)
        p0 = a*v(1,2) - b*v(1,3);			       	   
        p1 = a*v(2,2) - b*v(2,3);			       	   
        minp = min(p0,p1);
        maxp = max(p0,p1);
        rad = fa * hb(2) + fb * hb(3);   
        if minp > rad || maxp < -rad
            flag = 0;
            return;
        end
        flag = 1;
    end
%%% Y-Tests%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function flag = AxTestY02(a,b,fa,fb)
        p0 = -a*v(1,1) + b*v(1,3);			       	   
        p2 = -a*v(3,1) + b*v(3,3);			       	   
        minp = min(p0,p2);
        maxp = max(p0,p2);
        rad = fa * hb(1) + fb * hb(3);   
        if minp > rad || maxp < -rad
            flag = 0;
            return;
        end
        flag = 1;
    end
    function flag = AxTestY1(a,b,fa,fb)
        p0 = -a*v(1,1) + b*v(1,3);			       	   
        p1 = -a*v(2,1) + b*v(2,3);			       	   
        minp = min(p0,p1);
        maxp = max(p0,p1);
        rad = fa * hb(1) + fb * hb(3);   
        if minp > rad || maxp < -rad
            flag = 0;
            return;
        end
        flag = 1;
    end
%%% Z-Tests%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function flag = AxTestZ12(a,b,fa,fb)
        p1 = a*v(2,1) - b*v(2,2);			       	   
        p2 = a*v(3,1) - b*v(3,2);			       	   
        minp = min(p1,p2);
        maxp = max(p1,p2);
        rad = fa * hb(1) + fb * hb(2);   
        if minp > rad || maxp < -rad
            flag = 0;
            return;
        end
        flag = 1;
    end
    function flag = AxTestZ0(a,b,fa,fb)
        p0 = a*v(1,1) - b*v(1,2);			       	   
        p1 = a*v(2,1) - b*v(2,2);			       	   
        minp = min(p0,p1);
        maxp = max(p0,p1);
        rad = fa * hb(1) + fb * hb(2);   
        if minp > rad || maxp < -rad
            flag = 0;
            return;
        end
        flag = 1;
    end

    function dest = SUB(p1,p2)
        dest = zeros(1,3);
        dest(1) = p1(1)-p2(1);
        dest(2) = p1(2)-p2(2);
        dest(3) = p1(3)-p2(3);
    end

    function [center, h] = BoxCenter()
        % Function for finding the box center and size.
        h = [(boxEgdes(4)-boxEgdes(1))/2, (boxEgdes(5)-boxEgdes(2))/2, (boxEgdes(6)-boxEgdes(3))/2];
        center = [boxEgdes(1)+h(1), boxEgdes(2)+h(2), boxEgdes(3)+h(3)];
    end

    function flag = BBIntersection()
        % Checks if the triangle boundary box intersects the Bin-box.
%         triXmin = min(triVert(:,1));
%         triYmin = min(triVert(:,2));
%         triZmin = min(triVert(:,3));
%         triXmax = max(triVert(:,1));
%         triYmax = max(triVert(:,2));
%         triZmax = max(triVert(:,3));
%         triXmin = xt(1);
%         triYmin = yt(1);
%         triZmin = zt(1);
%         triXmax = xt(3);
%         triYmax = yt(3);
%         triZmax = zt(3);

        
%         if triXmin > boxEgdes(4) || triXmax < boxEgdes(1)
%             flag = 0;
%             return;
%         elseif triYmin > boxEgdes(5) || triYmax < boxEgdes(2)
%             flag = 0;
%             return;
%         elseif triZmin > boxEgdes(6) || triZmax < boxEgdes(3)
%             flag = 0;
%             return;
%         else 
%             flag = 1;
%         end
        if tri2(1,1) > boxEgdes(4) || tri2(1,3) < boxEgdes(1)
            flag = 0;
            return;
        elseif tri2(2,1) > boxEgdes(5) || tri2(2,3) < boxEgdes(2)
            flag = 0;
            return;
        elseif tri2(3,1) > boxEgdes(6) || tri2(3,3) < boxEgdes(3)
            flag = 0;
            return;
        else 
            flag = 1;
        end
    end

    function flag = planeBoxIntersect()
        n = cross(f(1,:),f(2,:));
        vmin = zeros(1,3);
        vmax = zeros(1,3);
        for i = 1:3
            vert = v(1,i);
            if n(i) > 0
                vmin(i) = -hb(i) - vert;
                vmax(i) = hb(i) - vert;
            else
                vmin(i) = hb(i) - vert;
                vmax(i) = -hb(i) - vert;
            end 
        end
        if dot(n,vmin) > 0
            flag = 0;
            return;
        elseif dot(n,vmax) >= 0
            flag = 1;
            return
        end
        flag = 0;
    end
end