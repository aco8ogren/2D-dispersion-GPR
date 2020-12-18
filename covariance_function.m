function C_interp = covariance_function(wv_i,wv_j,original_domain_X,original_domain_Y,original_covariance)
    M = length(original_domain_X);
    N = length(original_domain_Y);
    C4D = reshape(original_covariance,M,N,M,N);
    
    wv = combvec(wv_i',wv_j');
    C_interp_4D = interpn(original_domain_X,original_domain_Y,original_domain_X,original_domain_Y,C4D,wv(1,:),wv(2,:),wv(3,:),wv(4,:));
    
    [N_h,~] = size(wv_i);
    [M_h,~] = size(wv_j);
    
    C_interp = reshape(C_interp_4D,N_h,M_h);
end