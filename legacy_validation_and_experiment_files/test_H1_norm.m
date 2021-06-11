clear; close all;
f = @(x,y) x.^2 + y.^2;
% f = @(x,y) x.*sin(x) + cos(y);

[x,y] = meshgrid(linspace(0,1));
h_x = x(1,2)-x(1,1);
h_y = y(2,1)-y(1,1);

z = f(x,y);

[gzx,gzy] = gradient(z,h_x,h_y);
gz = cat(3,gzx,gzy);

surf(x,y,z)

compute_H1_norm(z,h_x,h_y)

