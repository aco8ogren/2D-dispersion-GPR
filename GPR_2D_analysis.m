clear; % close all;

% data = load("C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\2D-dispersion\OUTPUT\IBZ surface 1 output 09-Sep-2020 18-42-52\DATA N_struct500 09-Sep-2020 18-42-52.mat")
% data = load("C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\2D-dispersion\OUTPUT\IBZ surface 1 output 09-Sep-2020 18-42-52\DATA N_struct500 09-Sep-2020 18-42-52.mat",'EIGENVALUE_DATA','WAVEVECTOR_DATA','CONSTITUTIVE_DATA');
data = load("C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\2D-dispersion\OUTPUT\IBZ surface 3 output 11-Nov-2020 21-11-32\DATA N_struct100 RNG_offset0 11-Nov-2020 21-11-32.mat",'EIGENVALUE_DATA','WAVEVECTOR_DATA','CONSTITUTIVE_DATA');

EIGENVALUE_DATA = data.EIGENVALUE_DATA;
WAVEVECTOR_DATA = data.WAVEVECTOR_DATA;
CONSTITUTIVE_DATA = data.CONSTITUTIVE_DATA;

IBZ_shape = 'rectangle';
isUsePolar = true;

N_h = 20; % number of sample points sampled in the long direction of the rectangle for GPR
a = 1; % lattice parameter

[N_struct,N_eig,N_wv] = size(EIGENVALUE_DATA);
% N_k = get_N_k(N_wv); % tri
N_k = sqrt(N_wv); % rect
N_e = 3*N_k; % number of points to evaluate error on

e_L2 = zeros(N_struct,N_eig);
e_H1 =  zeros(N_struct,N_eig);

parfor struct_idx = 1:N_struct
    disp(['computing for structure ' num2str(struct_idx)])
    for eig_idx = 1:N_eig        
        wv = squeeze(WAVEVECTOR_DATA(struct_idx,:,:));
        fr = squeeze(EIGENVALUE_DATA(struct_idx,eig_idx,:))';
        
%         if strcmp(IBZ_shape,'tri')
%             Z(triu(true(N_k))) = squeeze(fr);
%             X(triu(true(N_k))) = squeeze(wv(1,:));
%             Y(triu(true(N_k))) = squeeze(wv(2,:));
%         elseif strcmp(IBZ_shape,'rectangle')
            Z = reshape(squeeze(fr),N_k,N_k);
            X = reshape(squeeze(wv(1,:)),N_k,N_k);
            Y = reshape(squeeze(wv(2,:)),N_k,N_k);
%         end
        
        [X_h,Y_h] = meshgrid(linspace(-pi/a,pi/a,N_h),linspace(0,pi/a,ceil(N_h/2)));
        [X_e,Y_e] = meshgrid(linspace(-pi/a,pi/a,N_e),linspace(0,pi/a,ceil(N_e/2)));
        
        Z_h = interp2(X,Y,Z,X_h,Y_h);
        Z_e = interp2(X,Y,Z,X_e,Y_e);
        
        wv_h = [X_h(true(size(X_h))),Y_h(true(size(Y_h)))]'; % these could be replaced by reshape operations
        fr_h = Z_h(true(size(Z_h))); % these could be replaced by reshape operations
        wv_e = [X_e(true(size(X_e))),Y_e(true(size(Y_e)))]';
        fr_e = Z_e(true(size(Z_e)))';
        
        if isUsePolar
            [wv_h_r,wv_h_theta] = cart2pol(wv_h(1,:),wv_h(2,:));
            wv_h_polar = [wv_h_r; wv_h_theta];
            model = fitrgp(wv_h_polar',fr_h','Sigma',1e-14,'ConstantSigma',true);
            [wv_e_r,wv_e_theta] = cart2pol(wv_e(1,:),wv_e(2,:));
            wv_e_polar = [wv_e_r; wv_e_theta];
            fr_pred = predict(model,wv_e_polar');
        else
            model = fitrgp(wv_h',fr_h','Sigma',1e-14,'ConstantSigma',true);
            fr_pred = predict(model,wv_e')';
        end
        
        
        Z_pred = reshape(fr_pred,ceil(N_e/2),N_e);
        
        h_x = X_e(1,2) - X_e(1,1); h_y = Y_e(2,1) - Y_e(1,1);
%         [dZdgamma_pred(1,:,:),dZdgamma_pred(2,:,:)] = gradient(Z_pred,h_x,h_y);
%         [dZdgamma_e(1,:,:),dZdgamma_e(2,:,:)] = gradient(Z_e,h_x,h_y);
        
        err = fr_pred - fr_e; % in vector format
        Z_err = Z_pred - Z_e; % in matrix format
        [grad_x,grad_y] = gradient(Z_err,h_x,h_y);
        derrdgamma = cat(1,grad_x,grad_y);
        
        e_L2(struct_idx,eig_idx) = sqrt(sum(Z_err.^2,'all')*h_x*h_y);
        e_H1(struct_idx,eig_idx) = sqrt(sum(sqrt(derrdgamma(1,:,:).^2 + derrdgamma(2,:,:).^2),'all')*h_x*h_y);
    end
end

figure
imagesc(e_L2)
colorbar
set(gca,'ColorScale','log')
title('L2 error')
xlabel('branch index')
ylabel('structure index')

figure
imagesc(e_H1)
colorbar
set(gca,'ColorScale','log')
title('H1 error')
xlabel('branch index')
ylabel('structure index')