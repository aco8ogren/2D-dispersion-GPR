clear; close all;

data_path = ['C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\'...
    '2D-dispersion\OUTPUT\ground_truth output 20-May-2021 17-11-07\DATA N_struct128 N_k201 RNG_offset0 20-May-2021 17-11-07.mat'];

[WAVEVECTOR_DATA,EIGENVALUE_DATA] = load_dispersion_dataset(data_path);

[N_wv,N_eig,N_struct] = size(EIGENVALUE_DATA);
N_k = sqrt(N_wv);

symmetry_loss = zeros(N_eig,N_struct);
boundary_loss = zeros(N_eig,N_struct);

for struct_idx = 1:N_struct
    for eig_idx = 1:N_eig
        X = reshape(WAVEVECTOR_DATA(:,1,struct_idx),N_k,N_k);
        Y = reshape(WAVEVECTOR_DATA(:,2,struct_idx),N_k,N_k);
        Z = reshape(EIGENVALUE_DATA(:,eig_idx,struct_idx),N_k,N_k);
        Z_flip = flipud(Z);
        symmetry_loss(eig_idx,struct_idx) = LP_norm(X',Y',(Z-Z_flip)',2);
        boundary_loss(eig_idx,struct_idx) = max(abs(Z(1,:) - Z(end,:)),[],'all');
        boundary_loss2(eig_idx,struct_idx) = max(abs(Z(2,:) - Z(end - 1,:)),[],'all');
%         figure
%         imagesc(Z)
%         figure

%         surf(X,Y,Z)
%         view(2)
%         pause
    end
end
        
figure
boxchart(symmetry_loss')

figure
boxchart(boundary_loss')

figure
boxchart(boundary_loss2')
