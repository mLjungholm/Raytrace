% Code for creation of the filter matrix

% Create a matrix for all photoreceptors (n) and all points in space (k).
% absorptionMat = [k*n]
% absorptionMat = zeros(3900,4063);

% I will first do a testrun doing the process for only a few photoreceptors
% to see if it works as intended. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load all coordinate files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Creating the abosrption matrix

% receptorNum = 2104;
% sourceNum = 824;
% absorptionMat = zeros(sourceNum,receptorNum);
% 
% [filename, path] = uigetfile('*','MultiSelect','on');
% fileNums = size(filename,2);
% % 
% % % loop through all grids and get the absorption for all photoreceptors
% % k = [1500,2500,3000,2700,3800,3915];
% % 
% h = waitbar(0,'Initializing waitbar...');
% perc = 0;
% % % onestep = 0;
% step = 100/fileNums;
% % 
% 
% for j = 1:fileNums
% %     Load the receptor grid.
%     s = LoadFile(strcat(path, char(filename(j))),'coneGrid'); 
%     for i = 1:receptorNum
%         absorptionMat(j,i) = s(i).absValue;
%     end
%     
%     perc = perc + step;
%     tempPerc = round(perc);
%     waitbar(tempPerc/100,h,sprintf('%d%% done...',tempPerc))
% end

%% Visualizing the 3D -> fisheye projection
% iSizeOr = 4032;
% radOr = (3937-48)/4032;

% % iSize = 128;
% iSize = 512;
% % 
% 
% I = zeros(iSize,iSize);
% k = [1500,2500,3000,2700,3800,3915];
% pointList = zeros(6,2);
% % nums = 4063-3600;
% % nums = 4063-0;
% % pointList = zeros(nums,2);
% for j = 1:6
%     [az,el,r] = cart2sph(coneGrid(k(j)).basep(1),coneGrid(k(j)).basep(2),coneGrid(k(j)).basep(3));
% %     [az,el,r] = cart2sph(coneGrid(j).basep(1),coneGrid(j).basep(2),coneGrid(j).basep(3));
% %     [az,el,r] = cart2sph(coneGrid(j+3600).basep(1),coneGrid(j+3600).basep(2),coneGrid(j+3600).basep(3));
%     [yi,xi] = sph2pixel([az,el],[iSize iSize],round(iSize*diamMod));
%     I(yi,xi) = 1;
%     pointList(j,:) = [xi,yi];
% end
% 
% figure(1)
% imagesc(I)
% axis equal


%% Creating an absorption matrix for one photoreceptor
% iSizeOr = 4032;
% radOr = (3937-48)/4032;
% iSize = 256;
% I = zeros(iSize,iSize);
% 
% 
% for j = 1:3900
% %     [az,el,r] = cart2sph(coords);
%     [yi,xi] = sph2pixel([coords(j,2)/100,coords(j,1)/100],[iSize iSize],round(iSize*radOr));
% %     [yi,xi] = sph2pixel([el,az],[iSize iSize],round(iSize*radOr));
%     I(yi,xi) = absorptionMat(j,2);
% end
% 
% figure(1)
% imagesc(I)

% Question: should the interpolation be done in 3d space or in the 2D
% fisheye projection? Will it make a diffenence?

%% Creating voronoi cells

% ! Load one receptor grid to get all the basepoints !

% Original image size and new image size
% iSizeOr = 4032;
% radOr = (3937-48)/4032;
% diamMod = 1.3423e+03/1080;
% iSize = 512;
% 
% I = zeros(iSize,iSize);
% pointList = zeros(receptorNum,2);
% for j = 1:receptorNum
%     [az,el,r] = cart2sph(rG(j).basep(1),rG(j).basep(2),rG(j).basep(3));
%     [yi,xi] = sph2pixel([az,el],[iSize iSize],round(iSize*radOr));
%     pointList(j,:) = [xi,yi];
% end
% 
% voronoiList = VoronoiAllocation(pointList,[iSize,iSize]);

% ! This code is alternative so to not need to flip the end image

% Original image size and new image size
% iSizeOr = 4032;
% radOr = (3937-48)/4032; % I do not realy know what this is? Diameter
% factor or something? 
% receptorNum = 4063;
% sourceNum = 3900;
% diamMod = 1.3423e+03/1080;
% iSize = 512;
% 
% I = zeros(iSize,iSize);
% pointList = zeros(receptorNum,2);
% for j = 1:receptorNum
%     [az,el,r] = cart2sph(rG(j).basep(1),rG(j).basep(2),rG(j).basep(3));
%     [yi,xi] = sph2pixel([-az,-el],[iSize iSize],round(iSize*diamMod));
%     pointList(j,:) = [xi,yi];
% end
% 
% voronoiList = VoronoiAllocation(pointList,[iSize,iSize]);


%% Visualizing the senistivity function of the photoreceptors
% t = 2500; % Receptor nr
% [Az, El] = meshgrid(-pi:0.01:pi,-pi:0.01:pi);
% % interpolate nonuniformly spaced points
% az = -coords(:,2)/100;
% el = coords(:,1)/100;
% r1 = ones(size(el,1),1)*120; 
% R = griddata(az,el,r1,Az,El);
% C = griddata(az,el,absorptionMat(:,t),Az,El).*10^6;
% % convert to cart
% [x, y, z] = sph2cart(Az,El,R);
% 
% figure(1)
% hold on
% colormap(viridis)
% % axis equal off vis3d
% axis equal
% grid on
% surface(-x,y,z,C,'edgealpha',0.05)
% quiver3(1.5,0,0,0.6,0,0,'r','linewidth',1.2)
% quiver3(0,1.3,0,0,0.6,0,'b','linewidth',1.2)
% quiver3(0,0,1.1,0,0,0.5,'b','linewidth',1.2)

%% Visualizing the back retina and selected photoreceptors

% figure(1)
% hold on
% axis equal
% k = [1500,2500,3000,2700,3800,3915];
% for i = t:t
% %     plot3(rG(k(i)).basep(1),rG(k(i)).basep(2),rG(k(i)).basep(3),'ro','linewidth',2)
%     plot3(rG(i).basep(1),rG(i).basep(2),rG(i).basep(3),'ro','linewidth',2)
% end
% trisurf(eye(1).f,eye(1).v(:,1),eye(1).v(:,2),eye(1).v(:,3),0.05,'facealpha',0.3,'edgealpha',0.5)
% trisurf(eye(2).f,eye(2).v(:,1),eye(2).v(:,2),eye(2).v(:,3),'Facecolor','b','facealpha',0.2,'edgealpha',0.2)
% trisurf(eye(3).f,eye(3).v(:,1),eye(3).v(:,2),eye(3).v(:,3),'Facecolor','b','facealpha',0.2,'edgealpha',0.2)
% trisurf(eye(4).f,eye(4).v(:,1),eye(4).v(:,2),eye(4).v(:,3),'Facecolor','g','facealpha',0.2,'edgealpha',0.2)
% trisurf(eye(5).f,eye(5).v(:,1),eye(5).v(:,2),eye(5).v(:,3),0.14,'facealpha',0.8,'edgealpha',0.2)

%% Visualizing the absorption function for one photoreceptor in the fisheye projection

% absMat2D = single(zeros(iSize,iSize,6));
% coordTransferMat = zeros(3900,2);
% iSizeOr = 4032;
% radOr = (3937-48)/4032;
% for i = 1:3900
%     [yi,xi] = sph2pixel([-coords(i,2)/100,coords(i,1)/100],[iSize iSize],round(iSize*radOr));
%     coordTransferMat(i,:) = [yi,xi];
% end
% 
% [yq,xq] = meshgrid(1:iSize, 1:iSize);
% for i = 1:6
%     vq = griddata(coordTransferMat(:,2),coordTransferMat(:,1),absorptionMat(:,i),yq,xq);
%     absMat2D(:,:,i) = vq;
% end

% figure(2)
% % hold on
% colormap(viridis)
% axis equal
% imagesc(absMat2D(:,:,t))
% imagesc(Icoords)

%% Loading an image and filtering it for the 6 testpoints

% [filename, path] = uigetfile('*');
% filePath = strcat(path,filename);
% I = imread(filePath);
% I = rgb2gray(I);   % Transform into gray-scale.
% 
% B=I;
% B = B(35:4792,:);
% B = B(:,745:5503);
% B = imresize(B,[iSize iSize]);
% C = im2single(B);
% % 
% % figure(1)
% % imshow(I)
% % 
% % figure(2)
% % imshow(B)
% 
% % Filter part
% 
% % Variables needed: absMat2D, voronoiList
% % absMat2D2 = absMat2D;
% % absMat2D2(isnan(absMat2D2)) = 0;
% 
% filterIm = zeros(iSize,iSize);
% 
% % imtest = absMat2D2(:,:,1).*C;
% for i = 1:6
%     absVal = sum(sum(absMat2D2(:,:,i).*C));
%     filterIm(voronoiList == i) = absVal;
% end
% filterIm = filterIm./max(max(filterIm));
% imshow(filterIm)

%% Creating an absorptionMatrix for all photoreceptors

% Workaround for the sizelimit of the matrix.
% a = single(zeros(iSize,iSize,2000));
% b = single(zeros(iSize,iSize,2063));
% absMat2D = cat(3,a,b);
% 
% % absMat2D = single(zeros(iSize,iSize,receptorNum));
% coordTransferMat = zeros(sourceNum,2);
% iSizeOr = 4032;
% radOr = (3937-48)/4032;
% 
% for i = 1:sourceNum
%     [yi,xi] = sph2pixel([-coords(i,2)/100,coords(i,1)/100],[iSize iSize],round(iSize*diamMod));
%     coordTransferMat(i,:) = [yi,xi];
% end
% 
% [yq,xq] = meshgrid(1:iSize, 1:iSize);
% 
% h = waitbar(0,'Initializing waitbar...');
% perc = 0;
% step = 100/receptorNum;
% for i = 1:receptorNum
%     vq = griddata(coordTransferMat(:,2),coordTransferMat(:,1),absorptionMat(:,i),yq,xq);
%     absMat2D(:,:,i) = vq;
%     
%     perc = perc + step;
%     tempPerc = round(perc);
%     waitbar(tempPerc/100,h,sprintf('%d%% done...',tempPerc))
% end
% absMat2D(isnan(absMat2D)) = 0;

%% Doing a complete filtering of an image with all photoreceptor

% iSize = 512;
% % [filename, path] = uigetfile('*');
% filePath = strcat(path,filename);
% I = imread(filePath);
% I = rgb2gray(I);   % Transform into gray-scale.
% % 
% B=I;
% B = B(35:4792,:);
% B = B(:,745:5503);
% B = imresize(B,[iSize iSize]);
% C = im2single(B);
% % % 
% % % figure(1)
% % % imshow(I)
% % % 
% 
% figure(2)
% imshow(C)
% 
% % 
% % Filter part
% 
% % Variables needed: absMat2D, voronoiList

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filter calibration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% receptorNum = 4063;
% sourceNum = 3900;


% calibIm = zeros(iSize,iSize);
% for i = 1:receptorNum
% %     if list(i) < 85
%     absVal = sum(sum(absMat2D(:,:,i)));
%     calibIm(voronoiList == i) = absVal;
% %     else
% %         calibIm(voronoiList == i) = 0;
% %     end
% end
% calibIm = calibIm./max(max(calibIm));
% 
% calibIm(calibIm == 0) = nan;
% inCalibIm = 1./calibIm;
% inClaibIm(inCalibIm == Inf) = 0;
% % 
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % Filter creation
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 
% filterIm = zeros(iSize,iSize);
% % % 
% % % % imtest = absMat2D2(:,:,1).*C.*inCalibIm;
% % tI = testIm;
% 
% C = imrotate(C,-90);
% C = flip(C,2);
% for i = 1:receptorNum
% %     if list(i) < 85
%     absVal = sum(sum(absMat2D(:,:,i).*C));
%     filterIm(voronoiList == i) = absVal;
% %     else
% %         filterIm(voronoiList == i) = 0;
% %     end
% end
% % filterIm = filterIm.*inCalibIm.*lightFilter;
% filterIm = filterIm.*inCalibIm;
% filterIm = filterIm./max(max(filterIm)).*255;
% filterIm = flip(filterIm,2);
% fI2 = uint8(imrotate(filterIm,90))-edgeFilter;
% % ! The images needs to be flipped because I made an error with the
% % coordinates of the phooreceptors vs source points. !
% % filterIm = flip(filterIm);
% % filterIm = fliplr(filterIm);
% 
% % figure(1)
% % imshow(filterIm)
% 
% figure(3)
% imshow(fI2)
% im3 = [B fI2];
% figure(4)
% imshow(im3)
% imwrite(im3,strcat(filePath(1:end-4),'_Filtered.tif'));

%% Edge of field of wiev calibration

% import the total absorbation file

% [yq,xq] = meshgrid(1:iSize, 1:iSize);
% vq = griddata(coordTransferMat(:,2),coordTransferMat(:,1),absVals,yq,xq);
% vq = vq./max(max(vq));
% vq(isnan(vq)) = 0;
% vq(vq > 0.35) = 1;
% lightFilter = vq./max(max(vq));
% figure(1)
% imagesc(vq)
% axis equal
% colormap(viridis)
%%
% tic
% fI2(fI2 == 0) = B(1,1);
% toc
% im3 = [B fI2];
% figure(4)
% imshow(im3)