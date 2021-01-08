function Cs = cov4D(WAVEVECTOR_DATA,EIGENVALUE_DATA)
    
    [N_wv,N_eig,~] = size(EIGENVALUE_DATA);
    N_k = sqrt(N_wv); % Contains assumption that input data is square grid (even if that means h_x ~= h_y, yet array of wv is still square)
    
    for eig_idx = 1:N_eig
        temp = cov(squeeze(EIGENVALUE_DATA(:,eig_idx,:))');
        Cs{eig_idx} = reshape(temp,N_k,N_k,N_k,N_k); % Contains assumption that input data is square grid (even if that means h_x ~= h_y, yet array of wv is still square)
    end    
end