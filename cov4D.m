function Cs = cov4D(WAVEVECTOR_DATA,EIGENVALUE_DATA)
    
    [N_wv,N_eig,~] = size(EIGENVALUE_DATA);
    N_k = sqrt(N_wv);
    
    for eig_idx = 1:N_eig
        temp = cov(squeeze(EIGENVALUE_DATA(:,eig_idx,:))');
        Cs{eig_idx} = reshape(temp,51,51,51,51);
    end    
end