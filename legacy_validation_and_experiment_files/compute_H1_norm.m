function out = compute_H1_norm(function_values,h_x,h_y)
    [Fx,Fy] = gradient(function_values,h_x,h_y);
    [M,N] = size(Fx);
    out = sqrt(sum(Fx.^2,'all') + sum(Fy.^2,'all'))*1/M*1/N;
end