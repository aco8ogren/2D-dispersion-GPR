function out = GPR2D(fr,wv,covariance,N_sample,N_evaluate,options)
    
    if ~exist('options','var')
        options = struct();
        options.isMakePlots = false;
        options.isUseEmpiricalCovariance = true;
    end
       
    
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
    fr_s = reshape(Z_s,1,[]);
    wv_e = [reshape(X_e,1,[]); reshape(Y_e,1,[])];
    fr_e = reshape(Z_e,1,[]);
    
    original_covariance = covariance;
    

    
    if options.isUseEmpiricalCovariance
        kfcn = @(wv_i,wv_j) covariance_function(wv_i,wv_j,original_domain_X,original_domain_Y,original_covariance);
        x_train = wv_s';
        y_train = fr_s';
        sigma = 1e-4;
        model = @(x_pred) create_model(x_pred,x_train,y_train,sigma,kfcn);
        fr_pred = model(wv_e')';
    else
        % Define a strict squared exponential so that GPR doesn't try to
        % optimize the fit with kernel parameters
        phi = [mean(std(wv_s'));std(fr_s')/sqrt(2)];        
        kfcn = @(XN,XM,theta) (phi(2)^2)*exp(-(pdist2(XN,XM).^2)/(phi(1)^2));
        
        model = fitrgp(wv_s',fr_s',...
            'Sigma',1e-14,...
            'ConstantSigma',true,...
            'KernelParameters',0,...
            'KernelFunction',kfcn,...
            'BasisFunction','none');
        sigma_L = phi(1);
        sigma_F = phi(2);
        covariance = [];        
        for i = 1:size(wv_s,2)
            for j = 1:size(wv_s,2)
                covariance(i,j) = sigma_F^2*exp(-.5*((wv_s(:,i) - wv_s(:,j))'*(wv_s(:,i) - wv_s(:,j)))/sigma_L^2);
            end
        end
        
        fr_pred = predict(model,wv_e')';
    end
    
    Z_pred = reshape(fr_pred,ceil(N_evaluate/2),N_evaluate);
    
    % err = fr_pred - fr_e; % in vector format
    Z_err = Z_pred - Z_e; % in matrix format
    [grad_x,grad_y] = gradient(Z_err,h_x,h_y);
    derrdgamma = cat(3,grad_x,grad_y);
    
    e_L2_old = sqrt(sum(Z_err.^2,'all')*h_x*h_y);
    e_L2 = LP_norm(X_e,Y_e,Z_err,2);
    e_H1_old = sqrt(sum(sqrt(derrdgamma(:,:,1).^2 + derrdgamma(:,:,2).^2),'all')*h_x*h_y);
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
    out.covariance = covariance;
    out.model = model;
end

function y_pred = create_model(x_pred,x_train,y_train,sigma,kfcn)
    % This function is intended to facilitate making function handles. I.e.
    % call like model = @(x_pred)
    % create_model(x_pred,x_train,y_train,sigma);
    y_pred = kfcn(x_pred,x_train)*...
        inv(kfcn(x_train,x_train) + sigma^2*eye(size(x_train,1),size(x_train,1)))*...
        y_train;
end