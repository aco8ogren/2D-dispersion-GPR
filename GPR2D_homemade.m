function out = GPR2D_homemade(fr,wv,N_wv,kfcn,N_sample,N_evaluate,options)
    
    if ~exist('options','var')
        options = struct();
        options.isMakePlots = false;
        options.isUseEmpiricalCovariance = true;
        options.sigma_GPR = 1e-3;
    end
       
    a = 1;
    
    X = reshape(wv(:,1),N_wv(1),N_wv(2))';
    Y = reshape(wv(:,2),N_wv(1),N_wv(2))';
    Z = reshape(fr,N_wv(1),N_wv(2))';
    
    original_domain_X = X(1,:);
    original_domain_Y = Y(:,1)';
    
    [X_s,Y_s] = meshgrid(linspace(-pi/a,pi/a,N_sample),linspace(0,pi/a,ceil(N_sample/2)));
    [X_e,Y_e] = meshgrid(linspace(-pi/a,pi/a,N_evaluate),linspace(0,pi/a,ceil(N_evaluate/2)));
        
    h_x = X_e(1,2) - X_e(1,1); h_y = Y_e(2,1) - Y_e(1,1);
    
    Z_s = interp2(X,Y,Z,X_s,Y_s);
    Z_e = interp2(X,Y,Z,X_e,Y_e);
    
    wv_s = [reshape(X_s,1,[]); reshape(Y_s,1,[])];
    fr_s = reshape(Z_s,1,[]);
    wv_e = [reshape(X_e,1,[]); reshape(Y_e,1,[])];
    fr_e = reshape(Z_e,1,[]);
    
%     original_covariance = covariance;
    
%     sigma = 1e-3;
    x_train = wv_s';
    y_train = fr_s';
    
    if options.isUseEmpiricalCovariance
        sigma = options.sigma_GPR;
        model = create_GPR_model3(x_train,y_train,sigma,kfcn,[],'gridded');
        fr_pred = model.pred(wv_e','gridded')';
    else
        % Define a strict squared exponential so that GPR doesn't try to
        % optimize the fit with kernel parameters
        sigma = options.sigma_GPR;
        phi = [mean(std(wv_s'));std(fr_s')/sqrt(2)];        
        kfcn = @(XN,XM,query_format) (phi(2)^2)*exp(-(pdist2(XN,XM).^2)/(phi(1)^2));
        model = create_GPR_model3(x_train,y_train,sigma,kfcn,[],'gridded');
        fr_pred = model.pred(wv_e','gridded')';
    end
    
    Z_pred = reshape(fr_pred,ceil(N_evaluate/2),N_evaluate);
    
    % err = fr_pred - fr_e; % in vector format
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
    out.X_s = X_s;
    out.Y_s = Y_s;
    out.Z_s = Z_s;
    out.wv_s = wv_s;
    out.original_domain_X = original_domain_X;
    out.original_domain_Y = original_domain_Y;
%     out.covariance = covariance;
    out.model = model;
end