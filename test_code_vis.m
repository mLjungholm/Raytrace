
% absorption_results = create_absorption_matrix();
% [filename, path] = uigetfile('*','MultiSelect','on');
% fileNums = size(filename,2);

close all
cornia.plot(1,'y');
lens.plot(1,'b');
retina.plot(1,'g')
visiualize_absorption(319,absorption_results,receptor_space);

%%

% input: absorption_results.

% Center ceceptor 1064
% Looking upp 319
% looking down 1766
% forward 1070
% back 1042



% receptor_id = 319;
% 
% [Az, El] = meshgrid(-90:1:90,-90:1:90);
% [az,el,~] = cart2sph(-absorption_results.source_coords(:,1),absorption_results.source_coords(:,2),absorption_results.source_coords(:,3));
% az = az.*180/pi;
% el = el.*180/pi;
% % % interpolate nonuniformly spaced points
% C = griddata(az,el,absorption_results.absorption_mat(:,receptor_id),Az,El);

close all

cM = find_centroid(C,0.5,0);
xLine = C(round(cM(1)),:);
xLine(isnan(xLine)) = 0;

figure(1)
plot(-90:1:90,xLine)

xLineF = fft(xLine);

Fs = 1000;            % Sampling frequency                    
T = 1/Fs;             % Sampling period       
L = 1500;             % Length of signal
t = (0:L-1)*T;        % Time vector

% plot(1000*t(1:50),X(1:50))
% title('Signal Corrupted with Zero-Mean Random Noise')
% xlabel('t (milliseconds)')
% ylabel('X(t)')

% Y = fft(X);
% Compute the two-sided spectrum P2. Then compute the single-sided spectrum P1 based on P2 and the even-valued signal length L.

P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
Define the frequency domain f and plot the single-sided amplitude spectrum P1. The amplitudes are not exactly at 0.7 and 1, as expected, because of the added noise. On average, longer signals produce better frequency approximations.

f = Fs*(0:(L/2))/L;
figure(2)
plot(xLineF)



% figure(1)
% colormap(viridis)
% ax = axes();
% ax.YDir = 'reverse';
% ax.XLim = [-90,90];
% ax.YLim = [-90,90];
% imagesc(ax,C)
% ax.YDir = 'reverse';
% % imagesc(ax,C);
% % % convert to cart
% [x, y, z] = sph2cart(Az,El,R);
% 
% figure(1)
% hold on
% colormap(viridis)
% 



%%
% maxP = [343.6 36, 800];
% steps = [3 5 10];
% stepS = (maxP-minP)./(steps-1);
% % [x,y,z] = meshgrid(minP(1):stepS(1):maxP(1),minP(2):stepS(2):maxP(2),minP(3):stepS(3):maxP(3));
% X = minP(1):stepS(1):maxP(1);
% 
% % The cell array is created with the y coordinate fist so that the matrix
% % indexing mathces the physical carteesian coordinate structure.
% receptor_grid = cell(steps(2),steps(1),steps(3)); 
% 
% receptor_grid{1,1,1} = [receptor_grid{1,1,1} 1];
% 
% receptor_grid{1,1,1} = [receptor_grid{1,1,1} 2];
% 
% 
% receptor_grid{1,1,1}
% receptor_grid{steps(2),steps(1),steps(3)} = []; 

% points = [reshape(x,[],1),reshape(y,[],1),reshape(z,[],1)];
% points = single(points);



%% Speed test

% point = [this.x_grid(x_ind) this.y_grid(y_ind) this.z_grid(z_ind)];
% point = rand([1,3]);
% aP = rand([2,3]);
% p_ind = 1;
% cone_dir = rand([1,3]);
% 
% tic
% for i = 1:100000
% l = sqrt((point(1)-aP(p_ind,1))^2 + (point(2)-aP(p_ind,2))^2 + (point(3)-aP(p_ind,3))^2);
% end
% toc
% 
% tic
% for i = 1:100000
% % Vector from cone apex to grid point
% pvec = (point-aP(p_ind,:))./norm(point-aP(p_ind,:));
% ang = acosd(pvec(1)*cone_dir(p_ind,1) + pvec(2)*cone_dir(p_ind,2) + pvec(3)*cone_dir(p_ind,3));
% 
% end
% toc

%%
% figure(1)
% retina_surface.plot(1,'r')
% % receptor_space.plot([200,2],'line',1)
% % ep = S.export_endpoints;
% % scatter3(ep(:,1),ep(:,2),ep(:,3))
% % receptor_space.plot('all','end',1)
% S.plot(1)


