

clear; close all;

data_path_train = ['C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\'...
    '2D-dispersion\OUTPUT\ground_truth output 20-May-2021 17-11-07\DATA N_struct128 N_k201 RNG_offset0 20-May-2021 17-11-07.mat'];

% data_path_train = ['C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\'...
%     '2D-dispersion-GPR\OUTPUT\Homog w dataset N_k51\DATA N_struct128 N_k51 RNG_offset0 14-Mar-2021 16-46-17.mat'];

% data_path_train = ['C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\'...
%     '2D-dispersion\OUTPUT\concatenated_dataset1\concatenated_dataset.mat'];
% covariance_options.N_wv = [201 201];

% data_path_train = ['C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\'...
%     '2D-dispersion\OUTPUT\covariance_singularity N_wv51x26 N_disp1000 04-Jun-2021 14-37-48\DATA N_struct1000 N_k RNG_offset0 04-Jun-2021 14-37-48.mat'];
% covariance_options.N_wv = [51 26];

% data_path_train = ['C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\'...
%     '2D-dispersion\OUTPUT\covariance_singularity N_wv31x16 N_disp1000 output 04-Jun-2021 15-35-59\DATA N_struct1000 N_k RNG_offset0 04-Jun-2021 15-35-59.mat'];
% covariance_options.N_wv = [31 16];

% data_path_train = ['C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\'...
%     '2D-dispersion\OUTPUT\covariance_singularity N_wv51x26 N_disp5000 output 04-Jun-2021 16-31-13\DATA N_struct5000 N_k RNG_offset0 04-Jun-2021 16-31-13.mat'];
% covariance_options.N_wv = [51 26];

% data_path_train = ['C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\'...
%     '2D-dispersion\OUTPUT\covariance_singularity N_wv51x26 N_disp20000 output 10-Jun-2021 14-58-54\DATA N_struct20000 N_k RNG_offset0 10-Jun-2021 14-58-54.mat'];

% data_path_train = ['C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\'...
%     '2D-dispersion\OUTPUT\covariance_singularity N_wv101x51 N_disp10000 output 11-Jun-2021 13-24-45\DATA N_struct10000 N_k RNG_offset0 11-Jun-2021 13-24-45.mat'];

N_sample = 14;
N_evaluate = [51 NaN]; N_evaluate(2) = ceil(N_evaluate(1)/2);
sigma = 0;

isPause = false;

[WAVEVECTOR_DATA,EIGENVALUE_DATA] = load_dispersion_dataset(data_path_train);

[N_wv,N_eig,N_struct] = size(EIGENVALUE_DATA);

covariance_options.isComputeCovarianceGradient = false;
covariance_options.isAllowGPU = false;
eig_idxs = 1:N_eig;

[X_e,Y_e] = meshgrid(linspace(-pi,pi,N_evaluate(1)),linspace(0,pi,N_evaluate(2)));
wv_e = [reshape(X_e,[],1) reshape(Y_e,[],1)];
[~,idxs] = sort(wv_e(:,1));
wv_e = wv_e(idxs,:);

figure
tiledlayout('flow')

for eig_idx_idx = eig_idxs
    eig_idx = eig_idxs(eig_idx_idx);
    covariance_options.eig_idxs = eig_idx;
    [Cs,C_grads,kfcns,kfcn_grads,X_grid_vec,Y_grid_vec] = get_empirical_covariance(WAVEVECTOR_DATA,EIGENVALUE_DATA,covariance_options);
    kfcn = kfcns{1};
    
    nexttile
    title(['eig\_idx = ' num2str(eig_idx)])
    wv_s = [];
    sample_idx_labels = {};
    %     hold on
    
%     C_evaluate = kfcn(wv_e,wv_e,'gridded');
%     covar = C_evaluate;
%     
%     variances = diag(covar);
%     V = reshape(variances,N_evaluate(2),N_evaluate(1));
    
    variances = zeros(size(wv_e,1),1);
    for wv_idx = 1:size(wv_e,1)
        variances(wv_idx) = kfcn(wv_e(wv_idx,:),wv_e(wv_idx,:),'gridded');
    end
    V = reshape(variances,N_evaluate(2),N_evaluate(1));

    cla
    hold on
    imagesc(wv_e(:,1),wv_e(:,2),V)
    set(gca,'YDir','normal')
    colorbar
    daspect([1 1 1])
    axis tight
    
    if isPause
        pause
    end
    
    for sample_idx = 1:N_sample
        
        [~,idx] = max(variances);
        wv_s = [wv_s; wv_e(idx,:)];
        
        sample_idx_labels{end + 1} = num2str(sample_idx);
        
        scatter(wv_s(:,1),wv_s(:,2),100,'r')
        text(wv_s(:,1),wv_s(:,2),sample_idx_labels,'HorizontalAlignment','center','VerticalAlignment','middle')
        
        if isPause
            pause
        end
        
%         C_sample = kfcn(wv_s,wv_s,'scattered');
%         C_mix1 = kfcn(wv_e,wv_s,'scattered');
%         C_mix2 = kfcn(wv_s,wv_e,'scattered');
%         C_evaluate = kfcn(wv_e,wv_e,'gridded');
%         
%         covar = C_evaluate - C_mix1*((C_sample + sigma^2*eye(size(C_sample)))\C_mix2);
%         
%         variances = diag(covar);
%         V = reshape(variances,N_evaluate(2),N_evaluate(1));
        
        variances = zeros(size(wv_e,1),1);
        for wv_idx = 1:size(wv_e,1)
            variances(wv_idx) = kfcn(wv_e(wv_idx,:),wv_e(wv_idx,:),'gridded') - kfcn(wv_e(wv_idx,:),wv_s,'scattered')*(kfcn(wv_s,wv_s,'scattered')\kfcn(wv_s,wv_e(wv_idx,:),'scattered'));
        end
        V = reshape(variances,N_evaluate(2),N_evaluate(1));

        cla
        hold on
        title(['eig\_idx = ' num2str(eig_idx)]) 
        imagesc(wv_e(:,1),wv_e(:,2),V)
        set(gca,'YDir','normal')
        scatter(wv_s(:,1),wv_s(:,2),100,'r')
        text(wv_s(:,1),wv_s(:,2),sample_idx_labels,'HorizontalAlignment','center','VerticalAlignment','middle')
        colorbar
        daspect([1 1 1])
        axis tight
        
        if isPause
            pause
        end
        
    end
end





