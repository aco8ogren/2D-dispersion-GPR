clear; close all;

data_path = ['C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\'...
    '2D-dispersion\OUTPUT\ground_truth output 20-May-2021 17-11-07\DATA N_struct128 N_k201 RNG_offset0 20-May-2021 17-11-07.mat'];

[WAVEVECTOR_DATA,EIGENVALUE_DATA] = load_dispersion_dataset(data_path);

[N_wv,N_eig,N_struct] = size(EIGENVALUE_DATA);
N_k = sqrt(N_wv);

symmetry_loss = zeros(N_eig,N_struct);
boundary_loss = zeros(N_eig,N_struct);
mean_abs_symmetry_loss_surface = zeros(N_k,N_k,N_eig);
        
% This grid is fixed, so can just use the first struct_idx
X = reshape(WAVEVECTOR_DATA(:,1,1),N_k,N_k);
Y = reshape(WAVEVECTOR_DATA(:,2,1),N_k,N_k);

% Reorient so that imagesc(Z) plots are oriented correctly
X = X';
Y = Y';
x = sort(unique(X));
y = sort(unique(Y));

for struct_idx = 1:N_struct
    for eig_idx = 1:N_eig
        Z = reshape(EIGENVALUE_DATA(:,eig_idx,struct_idx),N_k,N_k);
        
        % Reorient so that imagesc(Z) plots are oriented correctly
        Z = Z';
        
        Z_flip = fliplr(Z);
        symmetry_loss(eig_idx,struct_idx) = LP_norm(X',Y',(Z-Z_flip)',2);
        symmetry_loss_normalized(eig_idx,struct_idx) = LP_norm(X',Y',((Z-Z_flip)/mean(Z,'all'))',2);
        mean_abs_symmetry_loss_surface(:,:,eig_idx) = mean_abs_symmetry_loss_surface(:,:,eig_idx) + abs(Z-Z_flip);
%         if symmetry_loss_normalized(eig_idx,struct_idx) 
%             imagesc(flipud(Z'))
%             daspect([1 2 1])
%             colorbar
%             title(['dispersion relation, eig\_idx = ' num2str(eig_idx)])
%             xlabel('x1')
%             ylabel('x2')
%             
%             title('symmetry loss surface')
%             imagesc(flipud((Z-Z_flip)'))
%             colorbar
%             daspect([1 2 1])
%             twelve = 12;
%         end
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
mean_abs_symmetry_loss_surface = mean_abs_symmetry_loss_surface/N_struct;
        
figure
boxchart(symmetry_loss')
title('symmetry loss')

figure
boxchart(symmetry_loss_normalized')
title('symmetry loss normalized')

figure
boxchart(boundary_loss')
title('boundary loss')

figure
boxchart(boundary_loss2')
title('boundary loss2 (left side minus one row in from right side)')

figure
tiledlayout('flow')
for eig_idx = 1:N_eig
    nexttile
    imagesc(x,y,mean_abs_symmetry_loss_surface(:,:,eig_idx))
    title(['mean abs symmetry loss surface, eig\_idx = ' num2str(eig_idx)])
    daspect([1 1 1])
    xlabel('x1')
    ylabel('x2')
    set(gca,'YDir','normal')
    colorbar
end

figure
tiledlayout('flow')
tol = 5;
for eig_idx = 1:N_eig
    nexttile
    spy(flipud(mean_abs_symmetry_loss_surface(:,:,eig_idx)) < tol) % flipud because Z is upside-down
    title(['mean abs symmetry loss surface < tol = 1e' num2str(log10(tol)) ', eig\_idx = ' num2str(eig_idx)])
    daspect([1 2 1])
    xlabel('x1')
    ylabel('x2')
end


