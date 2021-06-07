clear; close all;

data_path_train = ['C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\'...
    '2D-dispersion\OUTPUT\ground_truth output 20-May-2021 17-11-07\DATA N_struct128 N_k201 RNG_offset0 20-May-2021 17-11-07.mat'];

% data_path_train = ['C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\'...
%     '2D-dispersion-GPR\OUTPUT\Homog w dataset N_k51\DATA N_struct128 N_k51 RNG_offset0 14-Mar-2021 16-46-17.mat'];

N_samples = 3:51;

% row_idxs = [100 101];

[WAVEVECTOR_DATA,EIGENVALUE_DATA] = load_dispersion_dataset(data_path_train);
[r_idxs,c_idxs] = find(WAVEVECTOR_DATA(:,2,:)~=0 | WAVEVECTOR_DATA(:,1,:) < 0);
WAVEVECTOR_DATA(r_idxs,:,:) = [];
EIGENVALUE_DATA(r_idxs,:,:) = [];

[N_wv,N_eig,N_struct] = size(EIGENVALUE_DATA);

covariance_options.isComputeCovarianceGradient = false;
covariance_options.isAllowGPU = false;
covariance_options.N_wv = [51 26];
eig_idxs = 1:N_eig;

ranks = zeros(length(eig_idxs),length(N_samples));
conds = zeros(length(eig_idxs),length(N_samples));
rconds = zeros(length(eig_idxs),length(N_samples));
sizes = zeros(length(eig_idxs),length(N_samples));
memories = zeros(length(eig_idxs),length(N_samples));
norms = zeros(length(eig_idxs),length(N_samples));

for eig_idx_idx = eig_idxs
    eig_idx = eig_idxs(eig_idx_idx);
%     covariance_options.eig_idxs = eig_idx;
%     [Cs,C_grads,kfcns,kfcn_grads,X_grid_vec,Y_grid_vec] = get_empirical_covariance(WAVEVECTOR_DATA,EIGENVALUE_DATA,covariance_options);
    
%     kfcn = kfcns{eig_idx};
    C = cov(squeeze(EIGENVALUE_DATA(:,eig_idx,:))');
    kfcn = @(xi,xj) kfcn_local(xi,xj,unique(sort(WAVEVECTOR_DATA(:,1,1))),C);
    for N_sample_idx = 1:length(N_samples)
        N_sample = N_samples(N_sample_idx);
        wv = linspace(0,pi,N_sample);
%         [X,Y] = meshgrid(linspace(-pi,pi,N_sample),linspace(0,pi,ceil(N_sample/2)));
%         X(:,end) = []; % chop off redundant points
%         Y(:,end) = []; % chop off redundant points
%         wv = [reshape(X,[],1) reshape(Y,[],1)];
%         [~,idxs] = sort(wv(:,2));
%         wv = wv(idxs,:);
        C = kfcn(wv,wv);
        sizes(eig_idx_idx,N_sample_idx) = size(C,1);
        ranks(eig_idx_idx,N_sample_idx) = rank(C);
        %         conds(eig_idx_idx,N_sample_idx) = cond(C);
        rconds(eig_idx_idx,N_sample_idx) = rcond(C);
        temp = whos('C');
        memories(eig_idx_idx,N_sample_idx) = temp.bytes/1e9;
        norms(eig_idx_idx,N_sample_idx) = norm(C(1,:) - C(2,:));
        norms_normalized(eig_idx_idx,N_sample_idx) = norm((C(1,:) - C(2,:))/max(C(1,:)));
    end
end

% figure
% imagesc(N_samples,eig_idxs,ranks)
% daspect([1,1,1])
% ylabel('eig\_idx')
% xlabel('N\_sample')
% title('Rank of covariance matrix')

figure
imagesc(N_samples,eig_idxs,ranks./sizes)
daspect([1,1,1])
ylabel('eig\_idx')
xlabel('N\_sample')
title('Rank/size of covariance matrix')
colorbar

figure
imagesc(N_samples,eig_idxs,log10(rconds))
daspect([1,1,1])
ylabel('eig\_idx')
xlabel('N\_sample')
title('log10(Reciprocal condition) of covariance matrix')
% set(gca,'colorscale','log')
colorbar

figure
imagesc(N_samples,eig_idxs,log10(norms))
daspect([1,1,1])
ylabel('eig\_idx')
xlabel('N\_sample')
title('log10(norm(row1 - row2)) of covariance matrix')
% set(gca,'colorscale','log')
colorbar

figure
imagesc(N_samples,eig_idxs,log10(norms_normalized))
daspect([1,1,1])
ylabel('eig\_idx')
xlabel('N\_sample')
title('log10(norm((row1 - row2)./max(row1))) of covariance matrix')
% set(gca,'colorscale','log')
colorbar

function C_interp = kfcn_local(Xi,Xj,original_domain,original_covariance)
[X,Y] = meshgrid(original_domain,original_domain);
V = original_covariance;
[Xq,Yq] = meshgrid(Xi,Xj);
Vq = interp2(X,Y,V,Xq,Yq);
C_interp = Vq';
end