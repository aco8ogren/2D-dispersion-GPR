clear; close all;

data = load("C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\2D-dispersion\OUTPUT\IBZ surface 3 output 11-Nov-2020 21-11-32\DATA N_struct100 RNG_offset0 11-Nov-2020 21-11-32.mat",'EIGENVALUE_DATA','WAVEVECTOR_DATA','CONSTITUTIVE_DATA');

EIGENVALUE_DATA = data.EIGENVALUE_DATA;
WAVEVECTOR_DATA = data.WAVEVECTOR_DATA;
CONSTITUTIVE_DATA = data.CONSTITUTIVE_DATA;

IBZ_shape = 'rectangle';
isUsePolar = true;

struct_idx = 1;
eig_idx = 1;
N_h = 7; % number of sample points sampled in the long direction of the rectangle for GPR
a = 1; % lattice parameter

[N_struct,N_eig,N_wv] = size(EIGENVALUE_DATA);
N_k = sqrt(N_wv); % rect
N_e = N_k; % number of points to evaluate error on

e_L2 = zeros(N_struct,N_eig);
e_H1 =  zeros(N_struct,N_eig);

wv = squeeze(WAVEVECTOR_DATA(struct_idx,:,:));
fr = squeeze(EIGENVALUE_DATA(struct_idx,eig_idx,:))';

fig = figure();
ax = axes(fig);
hold on
for eig_idx_iter = 1:8
    wv_iter = squeeze(WAVEVECTOR_DATA(struct_idx,:,:));
    fr_iter = squeeze(EIGENVALUE_DATA(struct_idx,eig_idx_iter,:))';
    plot_dispersion_surface(wv_iter,fr_iter,IBZ_shape,N_k,N_k,ax);
end
view(3)

Z = reshape(squeeze(fr),N_k,N_k);
X = reshape(squeeze(wv(1,:)),N_k,N_k);
Y = reshape(squeeze(wv(2,:)),N_k,N_k);
%         end

[X_h,Y_h] = meshgrid(linspace(-pi/a,pi/a,N_h),linspace(0,pi/a,ceil(N_h/2)));
[X_e,Y_e] = meshgrid(linspace(-pi/a,pi/a,N_e),linspace(0,pi/a,ceil(N_e/2)));

h_x = X_e(1,2) - X_e(1,1); h_y = Y_e(2,1) - Y_e(1,1);

Z_h = interp2(X,Y,Z,X_h,Y_h);
Z_e = interp2(X,Y,Z,X_e,Y_e);

wv_h = [X_h(true(size(X_h))),Y_h(true(size(Y_h)))]'; % these could be replaced by reshape operations
fr_h = reshape(Z_h,1,numel(Z_h)); %Z_h(true(size(Z_h))); % these could be replaced by reshape operations
wv_e = [X_e(true(size(X_e))),Y_e(true(size(Y_e)))]';
fr_e = Z_e(true(size(Z_e)))';

if isUsePolar
    [wv_h_theta,wv_h_r] = cart2pol(wv_h(1,:),wv_h(2,:));
    wv_h_polar = [wv_h_theta; wv_h_r];
    additional_wv_h_polar = [linspace(0,pi,N_h); zeros(1,N_h)]; additional_wv_h_polar = additional_wv_h_polar(:,2:end);
    wv_h_polar_train = [wv_h_polar additional_wv_h_polar];
    additional_fr_h = interp2(X,Y,Z,0,0)*ones(1,N_h-1); % Match whatever the sampled value is at the location of the singularity
    fr_h_train = [fr_h additional_fr_h];
    model = fitrgp(wv_h_polar_train',fr_h_train','Sigma',1e-14,'ConstantSigma',true);
    [wv_e_theta,wv_e_r] = cart2pol(wv_e(1,:),wv_e(2,:));
    wv_e_polar = [wv_e_theta; wv_e_r];
    fr_pred = predict(model,wv_e_polar')';
    make_plots(wv_e_polar,wv_h_polar,fr_e,fr_h,fr_pred,IBZ_shape,N_e,struct_idx,eig_idx,h_x,h_y,Z_e,'polar',additional_wv_h_polar,additional_fr_h)
else
    model = fitrgp(wv_h',fr_h','Sigma',1e-14,'ConstantSigma',true);
    fr_pred = predict(model,wv_e')';
end

make_plots(wv_e,wv_h,fr_e,fr_h,fr_pred,IBZ_shape,N_e,struct_idx,eig_idx,h_x,h_y,Z_e,'cartesian',[],[])

function make_plots(wv_e,wv_h,fr_e,fr_h,fr_pred,IBZ_shape,N_e,struct_idx,eig_idx,h_x,h_y,Z_e,plot_type,additional_wv_h_polar,additional_fr_h)
    % plot true dispersion relation
    fig = figure('position',[-1693 514 560 420]);
    ax = axes(fig);
    plot_dispersion_surface(wv_e,fr_e,IBZ_shape,N_e,ceil(N_e/2),ax);
    hold on
    scatter3(ax,wv_h(1,:),wv_h(2,:),fr_h,6,'r','MarkerFaceColor','r') % overlay the sample points
    title(['true dispersion relation' newline 'eig idx = ' num2str(eig_idx) ' struct idx = ' num2str(struct_idx)])
    if strcmp(plot_type,'polar')
        ylabel('$$\sqrt{\gamma_x^2 + \gamma_y^2}$$','interpreter','latex')
        xlabel('$$arctan(\frac{\gamma_y}{\gamma_x})$$','interpreter','latex')       
        scatter3(ax,additional_wv_h_polar(1,:),additional_wv_h_polar(2,:),additional_fr_h,6,'MarkerFaceColor','m','MarkerEdgeColor','m') % overlay additional sample points (polar needs supplemented points at r = 0)
    end
    
    % plot predicted dispersion relation
    fig = figure('position',[-1130 514 560 420]);
    ax = axes(fig);
    plot_dispersion_surface(wv_e,fr_pred,IBZ_shape,N_e,ceil(N_e/2),ax);
    hold on
    scatter3(ax,wv_h(1,:),wv_h(2,:),fr_h,6,'r','MarkerFaceColor','r') % overlay the sample points
    title(['predicted dispersion relation' newline 'eig idx = ' num2str(eig_idx) ' struct idx = ' num2str(struct_idx)])
    if strcmp(plot_type,'polar')
        ylabel('$$\sqrt{\gamma_x^2 + \gamma_y^2}$$','interpreter','latex')
        xlabel('$$arctan(\frac{\gamma_y}{\gamma_x})$$','interpreter','latex')
        scatter3(ax,additional_wv_h_polar(1,:),additional_wv_h_polar(2,:),additional_fr_h,6,'MarkerFaceColor','m','MarkerEdgeColor','m') % overlay additional sample points (polar needs supplemented points at r = 0)
    end
    
    % plot error surface
    fig = figure('position',[-567 514 560 420]);
    ax = axes(fig);
    plot_dispersion_surface(wv_e,fr_pred-fr_e,IBZ_shape,N_e,ceil(N_e/2),ax);
    hold on
    [max_abs_err, idx_abs] = max(abs(fr_pred-fr_e));
%     max_err_perc = abs(fr_pred(idx)-fr_e(idx))/fr_e(idx)*100;
    [max_err_perc,idx_perc] = max(100*abs(fr_pred-fr_e)./fr_e);
    s1 = scatter3(ax,wv_h(1,:),wv_h(2,:),max_abs_err*ones(size(wv_h(1,:))),6,'r','MarkerFaceColor','r'); s1.DisplayName = 'sample points'; % overlay the sample points
    
    [e_L2,e_H1] = get_error(fr_pred,Z_e,N_e,h_x,h_y);
    
    view(2)
    colorbar
    title(['error: L2 error = ' num2str(e_L2,4) ' H1 error = ' num2str(e_H1,4) newline 'max abs error = ' num2str(max_abs_err,4) ' Hz, max perc error = ' num2str(max_err_perc,3) '%' newline 'eig idx = ' num2str(eig_idx) ' struct idx = ' num2str(struct_idx)])
    if strcmp(plot_type,'polar')
        ylabel('$$\sqrt{\gamma_x^2 + \gamma_y^2}$$','interpreter','latex')
        xlabel('$$arctan(\frac{\gamma_y}{\gamma_x})$$','interpreter','latex')
        s4 = scatter3(ax,additional_wv_h_polar(1,:),additional_wv_h_polar(2,:),additional_fr_h,6,'MarkerFaceColor','m','MarkerEdgeColor','m'); s4.DisplayName = 'free degeneracy points'; % overlay additional sample points (polar needs supplemented points at r = 0)
    end
    
    s2 = scatter3(ax,wv_e(1,idx_abs),wv_e(2,idx_abs),fr_pred(idx_abs)-fr_e(idx_abs),50,'r','s'); s2.DisplayName = 'max abs err';
    s3 = scatter3(ax,wv_e(1,idx_perc),wv_e(2,idx_perc),fr_pred(idx_perc)-fr_e(idx_perc),50,'b','s'); s3.DisplayName = 'max perc err';   
    
    l = legend([s1,s2,s3],'Position',[0.7601 0.0222 0.2250 0.1262]);
    
end

function [e_L2,e_H1] = get_error(fr_pred,Z_e,N_e,h_x,h_y)
Z_pred = reshape(fr_pred,ceil(N_e/2),N_e);

% err = fr_pred - fr_e; % in vector format
Z_err = Z_pred - Z_e; % in matrix format
[grad_x,grad_y] = gradient(Z_err,h_x,h_y);
derrdgamma = cat(1,grad_x,grad_y);

e_L2 = sqrt(sum(Z_err.^2,'all')*h_x*h_y);
e_H1 = sqrt(sum(sqrt(derrdgamma(1,:,:).^2 + derrdgamma(2,:,:).^2),'all')*h_x*h_y);

disp(['e_L2 = ' num2str(e_L2)])
disp(['e_H1 = ' num2str(e_H1)])
end