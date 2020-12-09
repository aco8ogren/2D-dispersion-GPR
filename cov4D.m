function Cs = cov4D(WAVEVECTOR_DATA,EIGENVALUE_DATA)
    [N_wv,N_eig,N_struct] = size(EIGENVALUE_DATA);
    N_k = sqrt(N_wv);
    original_domain_X = linspace(-pi,pi,N_k);
    original_domain_Y = linspace(0,pi,N_k);
    for eig_idx = 1:N_eig
        fr_stack = squeeze(EIGENVALUE_DATA(:,eig_idx,:));
        wv = WAVEVECTOR_DATA(:,:,1);
        for X_idx1 = 1:length(original_domain_X)
            for Y_idx1 = 1:length(original_domain_Y)
                for X_idx2 = 1:length(original_domain_X)
                    for Y_idx2 = 1:length(original_domain_Y)
                        idx1 = find(all(wv == [original_domain_X(X_idx1) original_domain_Y(Y_idx1)],2));
                        idx2 = find(all(wv == [original_domain_X(X_idx2) original_domain_Y(Y_idx2)],2));
                        C4D(X_idx1,Y_idx1,X_idx2,Y_idx2) = sum((fr_stack(idx1,:) - mean(fr_stack(idx1,:))).*(fr_stack(idx2,:) -  mean(fr_stack(idx2,:))))/(N_struct - 1);
                    end
                end
            end
        end
        Cs{eig_idx} = C4D;
    end
    
end