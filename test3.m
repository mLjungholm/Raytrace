% plotting sine

tw = 20;
sw = 20/0.6667*2;
f = sw/90*pi;
shift = 0;
steps = 1000;
stepsize = sw*3/2/steps;
x = (0:stepsize:sw*3/2)./90.*pi;
psi = 0;
y = sin(2*pi*f*x + psi);

plot(x,y)