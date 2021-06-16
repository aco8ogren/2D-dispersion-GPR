clear; close all;

% data_path_train = ['C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\'...
%     '2D-dispersion\OUTPUT\ground_truth output 20-May-2021 17-11-07\DATA N_struct128 N_k201 RNG_offset0 20-May-2021 17-11-07.mat'];

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

data_path_train = ['C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\'...
    '2D-dispersion\OUTPUT\covariance_singularity N_wv101x51 N_disp10000 output 11-Jun-2021 13-24-45\DATA N_struct10000 N_k RNG_offset0 11-Jun-2021 13-24-45.mat'];

N_samples = 3:51;
sigmas = [ 0 1e-4 1e-3 1e-2 1e-1 ];
% kernel_type = 'squared exponential';

% row_idxs = [100 101];

[WAVEVECTOR_DATA,EIGENVALUE_DATA] = load_dispersion_dataset(data_path_train);

[N_wv,N_eig,N_struct] = size(EIGENVALUE_DATA);

covariance_options.isComputeCovarianceGradient = false;
covariance_options.isAllowGPU = false;
eig_idxs = 1:N_eig;

ranks = zeros(length(eig_idxs),length(N_samples),length(sigmas));
conds = zeros(length(eig_idxs),length(N_samples),length(sigmas));
rconds = zeros(length(eig_idxs),length(N_samples),length(sigmas));
sizes = zeros(length(eig_idxs),length(N_samples));
memories = zeros(length(eig_idxs),length(N_samples));
norms = zeros(length(eig_idxs),length(N_samples),length(sigmas));

for eig_idx_idx = eig_idxs
    eig_idx = eig_idxs(eig_idx_idx);
    covariance_options.eig_idxs = eig_idx;
    [Cs,C_grads,kfcns,kfcn_grads,X_grid_vec,Y_grid_vec] = get_empirical_covariance(WAVEVECTOR_DATA,EIGENVALUE_DATA,covariance_options);
    kfcn = kfcns{eig_idx};
    for N_sample_idx = 1:length(N_samples)
        N_sample = N_samples(N_sample_idx);
        [X,Y] = meshgrid(linspace(-pi,pi,N_sample),linspace(0,pi,ceil(N_sample/2)));
        X(:,end) = []; % chop off redundant points
        Y(:,end) = []; % chop off redundant points
        wv = [reshape(X,[],1) reshape(Y,[],1)];
        [~,idxs] = sort(wv(:,2));
        wv = wv(idxs,:);
        C = kfcn(wv,wv,'gridded');
%         C = nearestSPD(C);
        sizes(eig_idx_idx,N_sample_idx) = size(C,1);
        temp = whos('C');
        memories(eig_idx_idx,N_sample_idx) = temp.bytes/1e9;
        for sigma_idx = 1:length(sigmas)
            sigma = sigmas(sigma_idx);
            MAT = C + sigma^2*eye(size(C));
%             MAT = nearestSPD(MAT);
            ranks(eig_idx_idx,N_sample_idx,sigma_idx) = rank(MAT);
            rconds(eig_idx_idx,N_sample_idx,sigma_idx) = rcond(MAT);
            norms(eig_idx_idx,N_sample_idx,sigma_idx) = norm(MAT(1,:) - MAT(2,:));
            norms_normalized(eig_idx_idx,N_sample_idx,sigma_idx) = norm((MAT(1,:) - MAT(2,:))/max(MAT(1,:)));
        end
    end
end

fig = figure;
n = 3;
tiledlayout(length(sigmas),n);
for i = 1:length(sigmas)
    for j = 1:n
        nexttile
        ax(i,j) = gca;
    end
end

for sigma_idx = 1:length(sigmas)
    sigma = sigmas(sigma_idx);
    
    %     figure
    axes(ax(sigma_idx,1))
    imagesc(N_samples,eig_idxs,ranks(:,:,sigma_idx))
    daspect([1,1,1])
    ylabel('eig\_idx')
    xlabel('N\_sample')
    title(['Rank of covariance matrix + \sigma^2*I' newline '\sigma = 1e' num2str(log10(sigma))])
    colorbar
    
    %     figure
    axes(ax(sigma_idx,2))
    imagesc(N_samples,eig_idxs,ranks(:,:,sigma_idx)./sizes)
    daspect([1,1,1])
    ylabel('eig\_idx')
    xlabel('N\_sample')
    title(['Rank/size of covariance matrix + \sigma^2*I' newline '\sigma = 1e' num2str(log10(sigma))])
    colorbar
    
    %     figure
    axes(ax(sigma_idx,3))
    imagesc(N_samples,eig_idxs,log10(rconds(:,:,sigma_idx)))
    daspect([1,1,1])
    ylabel('eig\_idx')
    xlabel('N\_sample')
    title(['log10(Reciprocal condition) of covariance matrix + \sigma^2*I' newline '\sigma = 1e' num2str(log10(sigma))])
    % set(gca,'colorscale','log')
    colorbar
    
    %     figure
    %     imagesc(N_samples,eig_idxs,log10(norms(:,:,sigma_idx)))
    %     daspect([1,1,1])
    %     ylabel('eig\_idx')
    %     xlabel('N\_sample')
    %     title(['log10(norm(row1 - row2)) of covariance matrix + \sigma^2*I' newline '\sigma = 1e' num2str(log10(sigma))])
    %     % set(gca,'colorscale','log')
    %     colorbar
    %
    %     figure
    %     imagesc(N_samples,eig_idxs,log10(norms_normalized(:,:,sigma_idx)))
    %     daspect([1,1,1])
    %     ylabel('eig\_idx')
    %     xlabel('N\_sample')
    %     title(['log10(norm((row1 - row2)./max(row1))) of covariance matrix + \sigma^2*I' newline '\sigma = 1e' num2str(log10(sigma))])
    %     % set(gca,'colorscale','log')
    %     colorbar
    
end

fig2 = figure;
copyobj(fig.Children,fig2);

for j = [1 3]
    uniform_colorscale(ax(:,j));
end
set(ax(:,2),'clim',[0 1]);
