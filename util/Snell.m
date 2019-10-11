% Function for determine Recraction of a ray. Or reflection in the case of
% total internal reflection. Could be expanded to include power. 
function [vNew, reflected] = Snell(v, N, n1, n2)

% v and N needs to be unit vectors. This shuld be confirmd in higher code
% to minimize number of calculations. 
% N = N./norm(N);
% v = v./norm(v);

a = dot(v,N);
if a > 0
    N = -N;
end

c = -1*dot(N,v);
r = n1/n2;

vNew = r*v + (r*c - sqrt(1- r^2*(1-c^2)))*N;

if isreal(vNew) == 1    % If vNew is imagenary then total internal reflection occured. 
%     vRefracted = vNew;
    reflected = false;
else
    vNew = v - 2*dot(v,N).*N;   % Code for reflection.
%     vRefracted = vNew;
    reflected = true;
end

end
