clear; close all;

addpath('../2D-dispersion/')

data_path = "C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\2D-dispersion\OUTPUT\output 10-Jan-2022 14-40-30\DATA N_pix16x16 N_ele1x1 N_wv101x51 N_disp2000 N_eig10 offset0 10-Jan-2022 14-40-30.mat";
normalized_symmetry_loss_threshold = 1;

data = load(data_path);
EIGENVALUE_DATA = data.EIGENVALUE_DATA;
WAVEVECTOR_DATA = data.WAVEVECTOR_DATA;
N_wv = data.const.N_wv;

[~,N_eig,N_struct] = size(EIGENVALUE_DATA);

symmetry_loss = zeros(N_eig,N_struct);
boundary_loss = zeros(N_eig,N_struct);
mean_abs_symmetry_loss_surface = zeros([flip(N_wv) N_eig]);
        
% This grid is fixed, so can just use the first struct_idx
X = reshape(WAVEVECTOR_DATA(:,1,1),flip(N_wv)); 
Y = reshape(WAVEVECTOR_DATA(:,2,1),flip(N_wv));

% Reorient so that imagesc(Z) plots are oriented correctly
% X = X';
% Y = Y';
x = sort(unique(X));
y = sort(unique(Y));

list = [];
for struct_idx = 1:N_struct
    for eig_idx = 1:N_eig
        Z = reshape(EIGENVALUE_DATA(:,eig_idx,struct_idx),flip(N_wv));
        
        % Reorient so that imagesc(Z) plots are oriented correctly
%         Z = Z';
        
        Z_flip = fliplr(Z);
        symmetry_loss(eig_idx,struct_idx) = LP_norm(X,Y,(Z-Z_flip),2);
        symmetry_loss_normalized(eig_idx,struct_idx) = LP_norm(X,Y,((Z-Z_flip)/mean(Z,'all')),2);
        mean_abs_symmetry_loss_surface(:,:,eig_idx) = mean_abs_symmetry_loss_surface(:,:,eig_idx) + abs(Z-Z_flip);

        if symmetry_loss_normalized(eig_idx,struct_idx) > normalized_symmetry_loss_threshold
            list(end+1,:) = [struct_idx eig_idx];
        end
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


