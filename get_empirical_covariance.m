function [Cs,C_grads,kfcns,kfcn_grads,Vs,V_grads,vfcns,vfcn_grads,X_grid_vec,Y_grid_vec,covariance_info] = get_empirical_covariance(WAVEVECTOR_DATA,EIGENVALUE_DATA,covariance_options)
    variance_options.isComputeVarianceGradient = false; % HARD CODED
    
    [~,~,N_struct] = size(EIGENVALUE_DATA);

    X_grid_vec = unique(sort(WAVEVECTOR_DATA(:,1,1)));
    Y_grid_vec = unique(sort(WAVEVECTOR_DATA(:,2,1)));
    
    h_x = X_grid_vec(2) - X_grid_vec(1);
    h_y = Y_grid_vec(2) - Y_grid_vec(1);
    
    N_wv(1) = length(X_grid_vec);
    N_wv(2) = length(Y_grid_vec);
    
    original_C_struct.X_grid_vec = X_grid_vec;
    original_C_struct.Y_grid_vec = Y_grid_vec;
    original_C_struct.C = zeros([N_wv N_wv]);
    original_C_struct = repmat(original_C_struct,2,1);
    
    original_V_struct.X_grid_vec = X_grid_vec;
    original_V_struct.Y_grid_vec = Y_grid_vec;
    original_V_struct.V = zeros(N_wv);
    
    covariance_info.N_wv = N_wv;
    covariance_info.N_struct = N_struct;
    covariance_info.eig_idxs = covariance_options.eig_idxs;
    
    for i = 1:length(covariance_options.eig_idxs)
        eig_idx = covariance_options.eig_idxs(i);
        temp = cov(squeeze(EIGENVALUE_DATA(:,eig_idx,:))');
        Cs{i} = reshape(temp,N_wv(1),N_wv(2),N_wv(1),N_wv(2));
        Vs{i} = reshape(diag(temp),N_wv(1),N_wv(2));
        original_C_struct(1).C = Cs{i};
        original_V_struct(1).V = Vs{i};
        if covariance_options.isAllowGPU
%             original_C_struct.C_gpu = gpuArray(reshape(Cs{eig_idx},N_k,N_k,N_k,N_k));
        end
        kfcns{i} = @(wv_i,wv_j,query_format) covariance_function(wv_i,wv_j,query_format,original_C_struct(1),covariance_options);
        vfcns{i} = @(wv,query_format) variance_function(wv,query_format,original_V_struct(1),variance_options);
        if covariance_options.isComputeCovarianceGradient
            [C_wv_i2,C_wv_i1,~,~] = gradient(Cs{i},h_y,h_x,h_x,h_y); % [C_wv_i2,C_wv_i1,C_wv_j1,C_wv_j2] = gradient(Cs{i},h_y,h_x,h_x,h_y);
            C_grads{i} = cat(5,C_wv_i1,C_wv_i2);
            original_C_struct(1).C = C_wv_i1;
            original_C_struct(2).C = C_wv_i2;
            kfcn_grads{i} = @(wv_i,wv_j,grad_comp,query_format) covariance_function(wv_i,wv_j,query_format,original_C_grad_struct(grad_comp),covariance_options);
        end
    end
    
    if ~covariance_options.isComputeCovarianceGradient
        C_grads = {};
        kfcn_grads = {};
    end
    if ~variance_options.isComputeVarianceGradient
       V_grads = {};
       vfcn_grads = {};
    end
end
