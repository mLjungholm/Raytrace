% anglesRes = [41 65];
% attAngles = [12 16];
% 
% plot(attAngles, anglesRes)
% grid on

sp = [12.5 10 11 13 18 23 24 26 27 28 28 30 30 30 28 25 20 18 15 11.5 10.5 9.8 10.5 8.5 8.7];
elev = -55:5:65;

sppol = polyfit(elev,sp,5);
elev2 = -55:1:65;
spval = polyval(sppol,elev2);

cline = line([0 0],[5 35]); 
hline1 = line([60 60], [5 35]);
hline2 = line([40 40], [5 35]);
low50 = line([-55 -55], [5 35]);
high50 = line([75 75], [5 35]);
lowDrop = line([-40 -40], [5 35]);

figure(1)
hold on
plot(elev,sp,'+')
plot(elev2,spval)
grid on

xlabel('View angle [deg] from eye normal along "long-axis"')
ylabel('Estimated minimum viewing angle [deg]')
title('Spatial resolution as a function of view angle')