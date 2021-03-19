function [Cs,X_grid_vec,Y_grid_vec] = get_empirical_covariance(WAVEVECTOR_DATA,EIGENVALUE_DATA)    
    [N_wv,N_eig,N_struct] = size(EIGENVALUE_DATA);
    N_k = sqrt(N_wv);
    for eig_idx_iter = 1:N_eig
        temp = cov(squeeze(EIGENVALUE_DATA(:,eig_idx_iter,:))');
%         disp(['The rank of the 2D covariance matrix for the ' num2str(eig_idx_iter) ' branch is ' num2str(rank(temp)) ...
%                 '. Compare to size of matrix which is ' num2str(size(temp))])
        Cs{eig_idx_iter} = reshape(temp,N_k,N_k,N_k,N_k);
    end
    
%     min_wv_x = min(WAVEVECTOR_DATA(:,1,1));
%     max_wv_x = max(WAVEVECTOR_DATA(:,1,1));
%     min_wv_y = min(WAVEVECTOR_DATA(:,2,1));
%     max_wv_y = max(WAVEVECTOR_DATA(:,2,1));
    X_grid_vec = linspace(-pi,pi,N_k);
    Y_grid_vec = linspace(0,pi,N_k);
end
