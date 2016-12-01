%%
[X, Y] = meshgrid (-1:0.01:1, -1:0.01:1);
Z = X.^2 + Y.^2 + X.*Y;
% mesh(X,Y,Z)
contour(X,Y,Z)
hold on
[px, py] = gradient(Z,0.1,0.1);
quiver(px,py);
xlabel('x')
ylabel('y')

%% test gradient
x = 1;
y = 1;
z = x+y;
pxy = gradient(z,x,y) 

%% test quiver
x = 1;
y = 1;
px = 1;
py = 1;
quiver(x,y,px,py)

%% check convex
x = linspace(-2, 2, 100);
y = linspace(-2, 2, 100);
[X, Y] = meshgrid(x, y);
Z = X.^2 - 4 * X .* Y + Y.^2;
mesh(X, Y, Z)