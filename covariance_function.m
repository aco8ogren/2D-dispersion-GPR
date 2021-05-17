function C_interp = covariance_function(wv_i,wv_j,X_grid_vec,Y_grid_vec,original_covariance)
    M = length(X_grid_vec);
    N = length(Y_grid_vec);
    if numel(original_covariance)~=M*N*M*N
       warning('error in covariance_function') 
    end
    C4D = reshape(original_covariance,M,N,M,N);
    
    wv = combvec(wv_i',wv_j');
    C_interp_4D = interpn(X_grid_vec,Y_grid_vec,X_grid_vec,Y_grid_vec,C4D,wv(1,:),wv(2,:),wv(3,:),wv(4,:));
    
    [N_h,~] = size(wv_i);
    [M_h,~] = size(wv_j);
    
    C_interp = reshape(C_interp_4D,N_h,M_h);
end