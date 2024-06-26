function out = interp_model_2D(fr,wv,model_options)
    % Applies and analyzed error of a 2D interpolation model (including
    % GPR)
    
    N_sample = model_options.N_sample;
       
    a = 1;
    
    N_evaluate = [numel(unique(wv(:,1))) numel(unique(wv(:,2)))];
    
    X = reshape(wv(:,1),N_evaluate(1),N_evaluate(2))'; % Does the transpose fix the fact that I do the N_wv components out of order?
    Y = reshape(wv(:,2),N_evaluate(1),N_evaluate(2))';
    Z = reshape(fr,N_evaluate(1),N_evaluate(2))';
    
    original_domain_X = X(1,:);
    original_domain_Y = Y(:,1)';

    out.N_points = prod(N_sample);

    if strcmp(model_options.model_name,'GPR')
        isTrimRightBoundary = true;
    else
        isTrimRightBoundary = false;
    end
    
%     [X_s,Y_s] = get_wavevectors(N_sample,a,struct('isTrimRightBoundary',isTrimRightBoundary,'format','grid'));
    [X_e,Y_e] = get_wavevectors(N_evaluate,a,struct('isTrimRightBoundary',false,'format','grid'));
        
%     wv_s = get_wavevectors(N_sample,a,struct('isTrimRightBoundary',isTrimRightBoundary,'format','list'));
    wv_s = model_options.sample_points(1:model_options.N_sample,:);
    wv_e = get_wavevectors(N_evaluate,a,struct('isTrimRightBoundary',false,'format','list'));
    
    h_x = X_e(1,2) - X_e(1,1); h_y = Y_e(2,1) - Y_e(1,1);
    
    Z_s = interp2(X,Y,Z,wv_s(:,1),wv_s(:,2));
    Z_e = interp2(X,Y,Z,X_e,Y_e);
    if Z_e ~= Z
        warning('Evaluation points are being interpolated') % This is actually undesirable, and doesn't really need to be done if I have high resolution 'ground truth' datasets. Maybe interpolation of ground truth should never be done since it requires an interpolation model itself and could be inaccurate.
    end
    
    fr_s = reshape(Z_s,[],1);
    fr_e = reshape(Z_e,[],1);
      
    x_train = wv_s;
    y_train = fr_s;
    
    if strcmp(model_options.model_name,'GPR')
        model = create_GPR_model3(x_train,y_train,model_options.sigma,model_options.kfcn,model_options.kfcn_grad,'scattered');
        fr_pred = model.pred(wv_e,'scattered')';
        Z_pred = reshape(fr_pred,[N_evaluate(2) N_evaluate(1)]);
    else
        Z_pred = interp2(X_s,Y_s,Z_s,X_e,Y_e,model_options.model_name);
        model = {};
    end
    
    Z_err = Z_pred - Z_e; % in matrix format
    [grad_x,grad_y] = gradient(Z_err,h_x,h_y);
    derrdgamma = cat(3,grad_x,grad_y);
    
    e_L2 = LP_norm(X_e,Y_e,Z_err,2);
    e_H1 = H1_norm(X_e,Y_e,Z_err);
    
    out.e_L2 = e_L2;
    out.e_H1 = e_H1;
    out.Z_err = Z_err;
    out.Z_pred = Z_pred;
    out.X_e = X_e;
    out.Y_e = Y_e;
    out.Z_e = Z_e;
%     out.X_s = X_s;
%     out.Y_s = Y_s;
    out.Z_s = Z_s;
    out.wv_s = wv_s;
    out.original_domain_X = original_domain_X;
    out.original_domain_Y = original_domain_Y;
%     out.covariance = covariance;
    out.model = model;
end