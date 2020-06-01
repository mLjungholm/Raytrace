
close all
[Az, El] = meshgrid(-pi:0.01:pi,-pi/2:0.01:pi/2);

r = ones(absR.receptor_nums,1)*120; 
R = griddata(absR.view_dir(:,1)./180*pi,absR.view_dir(:,2)./180*pi,r,Az,El);
C = griddata(absR.view_dir(:,1)./180*pi,absR.view_dir(:,2)./180*pi,absR.horizontal_Cf,Az,El);

v = 0.015:0.01:0.15;
C2 = getContour(C,v);

% % convert to cart
[x, y, z] = sph2cart(Az,El,R);
figure(1)
hold on
% colormap(viridis)
% colormap(inferno)
axis equal off vis3d
surface(-x,y,z,C2,'edgealpha',0.05)
% cornia.plot(1,'y');
% lens.plot(1,'b');
% retina.plot(1,'g')


function outMat = getContour(matIn, scale_in)
[yi,xi] = size(matIn);
outMat = zeros(yi,xi);
for y_ind = 1:yi
    for x_ind = 1:xi
        if ~isnan(matIn(y_ind,x_ind))
            t = scale_in - matIn(y_ind,x_ind);
            [~,I] = min(abs(t));
            outMat(y_ind,x_ind) = scale_in(I);
        end
    end
end
end


%%
% absorption_results = create_absorption_matrix();
% [filename, path] = uigetfile('*','MultiSelect','on');
% fileNums = size(filename,2);

% close all
% cornia.plot(1,'y');
% lens.plot(1,'b');
% retina.plot(1,'g')
% visiualize_absorption(319,absorption_results,receptor_space);




%%

% input: absorption_results.

% Center ceceptor 1064
% Looking upp 319
% looking down 1766
% forward 1070
% back 1042
% forward is negative angles.



% receptor_id = 1064;
% % 
% close all
% [Az, El] = meshgrid(-90:1:90,-90:1:90);
% [az,el,~] = cart2sph(-absR.source_coords(:,1),absR.source_coords(:,2),absR.source_coords(:,3));
% az = az.*180/pi;
% el = el.*180/pi;
% % ind2angle = @(matInd,matSize) (matInd/matSize)*180-90;
% % C = griddata(az,el,absR.absorption_mat(:,receptor_id),Az,El);
% % C(isnan(C)) = 0;
% % figure(1)
% % cM = find_centroid(C,0.5,0);
% % angs = [ind2angle(cM(1),181),ind2angle(cM(2),181)];
% % figure(2)
% % imagesc(flipud(C))
% ind2angle = @(matInd,matSize) (matInd/matSize)*180-90; 
% absR.horizontal_Cf = zeros(absR.receptor_nums,1);
% absR.view_dir = zeros(absR.receptor_nums,2);
% for r_ind = 1:absR.receptor_nums 
% % interpolate nonuniformly spaced points
% C = griddata(az,el,absR.absorption_mat(:,r_ind),Az,El);
% C(isnan(C)) = 0;
% cM = find_centroid(C,0.5,0);
% absR.view_dir(r_ind,:) = [ind2angle(cM(1),181),ind2angle(cM(2),181)];
% absR.horizontal_Cf(r_ind) = CalculateCutoffFrequency(C,0,0);
% end



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


