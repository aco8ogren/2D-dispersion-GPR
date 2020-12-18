function out = H1_norm(X,Y,Z)
    h_x = X(1,2) - X(1,1); h_y = Y(2,1) - Y(1,1);    
    [Z_x,Z_y] = gradient(Z,h_x,h_y);
    
    out = sqrt(LP_norm(X,Y,Z_x,2)^2 + LP_norm(X,Y,Z_y,2)^2);
end