function cf = CalculateCutoffFrequency(projactedImage,direction,showPlot)
% Function for calculating the cutoff frequency in the choosen direction

% direction = 0 -> azhimutal direction
% direction = 1 -> elevation direction
cf = 0;
I = projactedImage;
I(isnan(I)) = 0; % Fill all empty spots with zero

if direction == 1
    I = imrotate(I,90);
end

% Reshape image from spacial coordinates to angular coordinates in the
% azimuthal direction.... On second thought. This is better done before the
% image is extrapolated.


stepSize = 1; % Set stepsize for the x-sampling grid
If = fft2(I);
If = fftshift(abs(If));
hX = round(size(I,2)/2);
% hY = round(size(I,1)/2);

maxVal = max(max(If));
[rowsOfMax, colsOfMax] = find(If == maxVal);

Ifl = If(rowsOfMax,:); 
% Ifl = Ifl(hX+1:end)./max(Ifl);
Ifl = Ifl(colsOfMax:end)./max(Ifl);
% f = (0:hX)/stepSize/size(I,2); % Determine frequency step
f = (0:colsOfMax)/stepSize/size(I,2);

cf = findCf(f,Ifl);

if showPlot == 1


    figure(1)
    plot(If(rowsOfMax,:))
    title('fourier transform in x-direction')
    xlabel('deg')
    
    figure(2)
    plot(Ifl)
    title('fourier transform in x-direction')
    xlabel('deg')
    
figure(3)
imagesc(I)
title('Projected absorption')

figure(4)
imagesc(If)
title('Forier transform image')

figure(5)
% plot(f(1:round(hX/2)),Ifl(1:round(hX/2)),'r')
plot(f(1:colsOfMax/2),Ifl(1:colsOfMax/2),'r')
title('Foirer transorm in the azimuthal direction')
xlabel('cycles per deg')

    figure(6)
    [rowsOfMax, colsOfMax] = find(I == max(max(I)));
    plot(I(rowsOfMax,:))
    title('absorption function before fourier transform')

end

function Cf = findCf(xGrid,f)
halfMax = max(f)/2;
% oldf = f(1); newf = 0; oldx = xGrid(1); newx = 0;
for i = 2:size(f,2)
    newf = f(i);
    newx = xGrid(i);
    oldf = f(i-1); oldx = xGrid(i-1);
    if newf <= halfMax
%         df = newf-oldf;
%         dx = newx-oldx;
        a = (newf-oldf)/(newx-oldx);
        b = newf - (a*newx);
        Cf = (halfMax - b)/a*2;
        return;
    else
        Cf = 0;
    end
end
disp('ERROR: no value found');
end

end