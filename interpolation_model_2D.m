function out = interpolation_model_2D(fr,wv,N_sample,N_evaluate,interpolation_method)
       
    a = 1;
    
    N_wv = size(wv,1);
    % N_k = get_N_k(N_wv); % tri
    N_k = sqrt(N_wv); % rect
    
    X = reshape(wv(:,1),N_k,N_k)';
    Y = reshape(wv(:,2),N_k,N_k)';
    Z = reshape(fr,N_k,N_k)';
    
    original_domain_X = X(1,:);
    original_domain_Y = Y(:,1)';
    
    [X_s,Y_s] = meshgrid(linspace(-pi/a,pi/a,N_sample),linspace(0,pi/a,ceil(N_sample/2)));
    [X_e,Y_e] = meshgrid(linspace(-pi/a,pi/a,N_evaluate),linspace(0,pi/a,ceil(N_evaluate/2)));
    
    h_x = X_e(1,2) - X_e(1,1); h_y = Y_e(2,1) - Y_e(1,1);
    
    Z_s = interp2(X,Y,Z,X_s,Y_s);
    Z_e = interp2(X,Y,Z,X_e,Y_e);
    
    wv_s = [reshape(X_s,1,[]); reshape(Y_s,1,[])];
    
    Z_pred = interp2(X_s,Y_s,Z_s,X_e,Y_e,interpolation_method);

    Z_err = Z_pred - Z_e; % in matrix format
%     [Z_err_x,Z_err_y] = gradient(Z_err,h_x,h_y);
%     derrdgamma = cat(3,Z_err_x,Z_err_y);
    
    e_L2 = LP_norm(X_e,Y_e,Z_err,2);
    e_H1 = H1_norm(X_e,Y_e,Z_err);
    
    out.e_L2 = e_L2;
    out.e_H1 = e_H1;
    out.Z_err = Z_err;
    out.Z_pred = Z_pred;
    out.X_e = X_e;
    out.Y_e = Y_e;
    out.Z_e = Z_e;
    out.X_s = X_s;
    out.Y_s = Y_s;
    out.Z_s = Z_s;
    out.wv_s = wv_s;
    out.original_domain_X = original_domain_X;
    out.original_domain_Y = original_domain_Y;
%     out.model = model;
end