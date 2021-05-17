function [Cs,C_grads,kfcns,kfcn_grads,X_grid_vec,Y_grid_vec] = get_empirical_covariance(WAVEVECTOR_DATA,EIGENVALUE_DATA)
    [N_wv,N_eig,N_struct] = size(EIGENVALUE_DATA);
    N_k = sqrt(N_wv);
    
    X_grid_vec = linspace(-pi,pi,N_k);
    Y_grid_vec = linspace(0,pi,N_k);
    
    for eig_idx_iter = 1:N_eig
        temp = cov(squeeze(EIGENVALUE_DATA(:,eig_idx_iter,:))');
        %         disp(['The rank of the 2D covariance matrix for the ' num2str(eig_idx_iter) ' branch is ' num2str(rank(temp)) ...
        %                 '. Compare to size of matrix which is ' num2str(size(temp))])
        Cs{eig_idx_iter} = reshape(temp,N_k,N_k,N_k,N_k);
        kfcns{eig_idx_iter} = @(wv_i,wv_j) covariance_function(wv_i,wv_j,X_grid_vec,Y_grid_vec,Cs{eig_idx_iter});
    end
    
    %     min_wv_x = min(WAVEVECTOR_DATA(:,1,1));
    %     max_wv_x = max(WAVEVECTOR_DATA(:,1,1));
    %     min_wv_y = min(WAVEVECTOR_DATA(:,2,1));
    %     max_wv_y = max(WAVEVECTOR_DATA(:,2,1));
    
    h_x = X_grid_vec(2) - X_grid_vec(1);
    h_y = Y_grid_vec(2) - Y_grid_vec(1);
    for eig_idx_iter = 1:N_eig
        [C_wv_i2,C_wv_i1,C_wv_j1,C_wv_j2] = gradient(Cs{eig_idx_iter},h_y,h_x,h_x,h_y);
        C_grads{eig_idx_iter} = cat(5,C_wv_i1,C_wv_i2);
        kfcn_grads{eig_idx_iter} = @(wv_i,wv_j,grad_comp) covariance_function(wv_i,wv_j,X_grid_vec,Y_grid_vec,C_grads{eig_idx_iter}(:,:,:,:,grad_comp));
    end
end
