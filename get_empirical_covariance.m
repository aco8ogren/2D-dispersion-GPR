function [Cs,C_grads,kfcns,kfcn_grads,X_grid_vec,Y_grid_vec] = get_empirical_covariance(WAVEVECTOR_DATA,EIGENVALUE_DATA,covariance_options)
    [N_wv,N_eig,N_struct] = size(EIGENVALUE_DATA);
    N_k = sqrt(N_wv);
%     N_wv = covariance_options.N_wv;
    
%     X_grid_vec = linspace(-pi,pi,N_k);
%     Y_grid_vec = linspace(0,pi,N_k);

%     X_grid_vec = linspace(-pi,pi,N_wv(1));
%     Y_grid_vec = linspace(0,pi,N_wv(2));

    X_grid_vec = unique(sort(WAVEVECTOR_DATA(:,1,1)));
    Y_grid_vec = unique(sort(WAVEVECTOR_DATA(:,2,1)));
    
    N_wv(1) = length(X_grid_vec);
    N_wv(2) = length(Y_grid_vec);
    
    for eig_idx = covariance_options.eig_idxs
        temp = cov(squeeze(EIGENVALUE_DATA(:,eig_idx,:))');
%         temp = cov(gpuArray(squeeze(EIGENVALUE_DATA(:,eig_idx,:))'));
%         Cs{eig_idx} = reshape(temp,N_k,N_k,N_k,N_k);
        Cs{eig_idx} = reshape(temp,N_wv(1),N_wv(2),N_wv(1),N_wv(2)); % is this correct?
%         kfcns{eig_idx_iter} = @(wv_i,wv_j) covariance_function(wv_i,wv_j,X_grid_vec,Y_grid_vec,Cs{eig_idx_iter});
        original_C_struct.X_grid_vec = X_grid_vec;
        original_C_struct.Y_grid_vec = Y_grid_vec;
        original_C_struct.C = Cs{eig_idx};
        if covariance_options.isAllowGPU
            original_C_struct.C_gpu = gpuArray(reshape(Cs{eig_idx},N_k,N_k,N_k,N_k));
        end
        kfcns{eig_idx} = @(wv_i,wv_j,query_format) covariance_function(wv_i,wv_j,query_format,original_C_struct,covariance_options);
    end
    
    %     min_wv_x = min(WAVEVECTOR_DATA(:,1,1));
    %     max_wv_x = max(WAVEVECTOR_DATA(:,1,1));
    %     min_wv_y = min(WAVEVECTOR_DATA(:,2,1));
    %     max_wv_y = max(WAVEVECTOR_DATA(:,2,1));
    
    if covariance_options.isComputeCovarianceGradient
    h_x = X_grid_vec(2) - X_grid_vec(1);
    h_y = Y_grid_vec(2) - Y_grid_vec(1);
    for eig_idx = 1:N_eig
        [C_wv_i2,C_wv_i1,C_wv_j1,C_wv_j2] = gradient(Cs{eig_idx},h_y,h_x,h_x,h_y);
        C_grads{eig_idx} = cat(5,C_wv_i1,C_wv_i2);
        kfcn_grads{eig_idx} = @(wv_i,wv_j,grad_comp) covariance_function(wv_i,wv_j,X_grid_vec,Y_grid_vec,C_grads{eig_idx}(:,:,:,:,grad_comp));
    end
    else
        C_grads = {};
        kfcn_grads = {};
    end
end
