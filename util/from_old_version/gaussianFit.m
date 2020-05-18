% Code to fit a 2d Gaussian to an immage. Thos one WOrks!!!!

[h,w] = size(C2);

[X,Y] = meshgrid(1:h,1:w);
X = X(:); Y=Y(:); Z = C2(:);
figure(1); clf; scatter3(X,Y,Z);

% 2D gaussian fit object
gauss2 = fittype( @(a1, sigmax, sigmay, x0,y0, x, y) a1*exp(-(x-x0).^2/(2*sigmax^2)-(y-y0).^2/(2*sigmay^2)),...
'independent', {'x', 'y'},'dependent', 'z' );

a1 = max(C2(:)); % height, determine from image. may want to subtract background
sigmax = 2; % guess width
sigmay = 2; % guess width
x0 = w/2; % guess position (center seems a good place to start)
y0 = h/2; 


% compute fit
sf = fit([X,Y],double(Z),gauss2,'StartPoint',[a1, sigmax, sigmay, x0,y0]);
figure(6); clf; 
% plot(sf,[X,Y],Z);
plot(sf);

figure(2)
surface(C2,'edgealpha',0.2)
% colormap('jet')

% sf.x0 and sf.y0 is the center of gaussian.
% sf.sigmax etc will get you the other parameters. 