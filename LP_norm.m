function out = LP_norm(X,Y,Z,p)
    % Only works for scalar valued functions Z defined on uniform grids X,Y
    % with spacings h_x,h_y
    h_x = X(1,2) - X(1,1); h_y = Y(2,1) - Y(1,1);
    out = (sum(Z.^p,'all')*h_x*h_y)^(1/p);
end