function out = interpolation_model_2D(fr,wv,N_wv,N_sample,N_evaluate,interpolation_method)
       
    a = 1;

    X = reshape(wv(:,1),N_wv(1),N_wv(2))'; % Does the transpose fix the fact that I do the N_wv components out of order?
    Y = reshape(wv(:,2),N_wv(1),N_wv(2))';
    Z = reshape(fr,N_wv(1),N_wv(2))';
    
    original_domain_X = X(1,:);
    original_domain_Y = Y(:,1)';
    
    out.N_points = prod(N_sample) - 2*N_sample(1) - N_sample(2) + 2; % This comes prior to updating N_sample because it's meant to track how many points were paid for in dispersion computation.

%     N_sample = N_sample + [1 0]; % Allow interpolation model to sample an extra column of points on the right side in the X direction. Because if we assume periodicity, these points can be obtained freely by copying the points on the left side.
    [X_s,Y_s] = get_wavevectors(N_sample,a,struct('isTrimRightBoundary',false,'format','grid'));
    [X_e,Y_e] = get_wavevectors(N_evaluate,a,struct('isTrimRightBoundary',false,'format','grid'));

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