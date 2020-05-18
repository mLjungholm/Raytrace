n = 2000;
fi = (sqrt(5)+1)/2 - 1;
K = fi * pi;
fib_points = zeros(n,3);
for i = 1:n
    long = K*i;
    lat = asin(-1 + 2*i/n);
    [x,y,z] = sph2cart(long,lat,1);
    fib_points(i,:) = [x,y,z];
end

half_sphere = 1;
if half_sphere
    inds = fib_points(:,1) > 0;
    fib_points = fib_points(inds,:);
end

endnr = n/2;
figure(4)
scatter3(fib_points(:,1),fib_points(:,2),fib_points(:,3),'.')
axis equal