clear; close all;

isPlot = false;

data_path_train = ['C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\2D-dispersion\OUTPUT\'...
    'light_dataset output 01-Jul-2021 17-32-47\DATA N_disp1000 N_wv51x26 01-Jul-2021 17-32-47.mat'];

% data_path_train = ['C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\2D-dispersion\OUTPUT\'...
%     'gold_dataset output 11-Jun-2021 13-24-45\DATA N_wv101x51 N_disp10000 RNG_offset0 11-Jun-2021 13-24-45.mat'];

% data_path_test = ['C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\2D-dispersion\OUTPUT\'...
%     'gold_dataset_small output 08-Jul-2021 12-54-14\DATA N_disp100 N_wv101x51 RNG_offset0 08-Jul-2021 12-54-14.mat'];

N_wv_test = [1001 501];

sample_order_path = ['C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\2D-dispersion-GPR\OUTPUT\sample_orders\'...
    'sample_order_information_LD_1326_1.mat'];

disp('Loading training set...')
[WAVEVECTOR_DATA_TRAIN,EIGENVALUE_DATA_TRAIN] = load_dispersion_dataset(data_path_train);
disp('done.')

% disp('Loading test set...')
% [WAVEVECTOR_DATA_TEST,EIGENVALUE_DATA_TEST] = load_dispersion_dataset(data_path_test);
% disp('done.')

disp('Loading sample order...')
sample_order_data = load(sample_order_path);
disp('done.')

eig_idxs = 1;
covariance_options.isAllowGPU = false;
covariance_options.isComputeCovarianceGradient = false;
covariance_options.eig_idxs = eig_idxs;
[Cs,C_grads,kfcns,kfcn_grads,Vs,V_grads,vfcns,vfcn_grads,X_grid_vec,Y_grid_vec,covariance_info] = ...
    get_empirical_covariance(WAVEVECTOR_DATA_TRAIN,EIGENVALUE_DATA_TRAIN,covariance_options);

sample_count = 690;
eig_idx_idx = 1;
eig_idx = eig_idxs(eig_idx_idx);

wv_s = sample_order_data.sampling_orders(1:sample_count,:,eig_idx);

X_grid_vec_train = sort(unique(WAVEVECTOR_DATA_TRAIN(:,1,1)));
Y_grid_vec_train = sort(unique(WAVEVECTOR_DATA_TRAIN(:,2,1)));

% X_grid_vec_test = sort(unique(WAVEVECTOR_DATA_TEST(:,1,1)));
% Y_grid_vec_test = sort(unique(WAVEVECTOR_DATA_TEST(:,2,1)));

X_grid_vec_test = linspace(-pi,pi,N_wv_test(1));
Y_grid_vec_test = linspace(0,pi,N_wv_test(2))';

[X_temp,Y_temp] = meshgrid(X_grid_vec_test,Y_grid_vec_test);
wv_e = [reshape(X_temp,[],1) reshape(Y_temp,[],1)];

[X_temp,Y_temp] = meshgrid(X_grid_vec_train,Y_grid_vec_train);
wv_o = [reshape(X_temp,[],1) reshape(Y_temp,[],1)];

kfcn = kfcns{eig_idx_idx};

tic
C_interp_comp = kfcn(wv_s,wv_e,'scattered');
toc

if isPlot
    figure
    tiledlayout('flow')
    for sample_idx = 1:sample_count
        nexttile
        plot(C_interp_comp(sample_idx,:))
        title(['Covariance using standard interpolation,' newline 'sample\_idx = ' num2str(sample_idx)])
    end
end

tic
% [val_x,idxs_x] = intersect(X_grid_vec_train,wv_s(:,1));
% [val_y,idxs_y] = intersect(Y_grid_vec_train,wv_s(:,2));

[~,idxs_x] = ismember(wv_s(:,1),X_grid_vec_train);
[~,idxs_y] = ismember(wv_s(:,2),Y_grid_vec_train);
% idxs_x(~idxs_x) = [];
% idxs_y(~idxs_y) = [];

C = Cs{eig_idx_idx};
lin_idxs = sub2ind(size(squeeze(C(1,1,:,:))),idxs_x,idxs_y);
C_sub = C(:,:,lin_idxs);
C_sub = permute(C_sub,[2 1 3]);

% C_sub = interpn(X_grid_vec_train,Y_grid_vec_train,X_grid_vec_train,Y_grid_vec_train,...
%     C,...
%     wv_o(:,1),wv_o(:,2),wv_s(:,1),wv_s(:,2));

    

% [vals,idxs] = intersect([val_x val_y],wv_s,'rows','stable');

C_interp = interp3(X_grid_vec_train,Y_grid_vec_train,1:sample_count,C_sub,X_grid_vec_test,Y_grid_vec_test,1:sample_count);
C_interp = reshape(C_interp,prod(size(C_interp,1,2)),size(C_interp,3));
C_interp = C_interp';
toc

if isPlot
    figure
    tiledlayout('flow')
    for sample_idx = 1:sample_count
        nexttile
        plot(C_interp(sample_idx,:))
        title(['Covariance using pre-sliced interpolation,' newline 'sample\_idx = ' num2str(sample_idx)])
    end
    
    figure
    tiledlayout('flow')
    for sample_idx = 1:sample_count
        nexttile
        plot(abs(C_interp_comp(sample_idx,:) - C_interp(sample_idx,:)))
        title(['Covariance using pre-sliced interpolation,' newline 'sample\_idx = ' num2str(sample_idx)])
    end
end

disp(['Frobenius deviation between interpolations is ' num2str(norm(C_interp - C_interp_comp,'fro'))])

