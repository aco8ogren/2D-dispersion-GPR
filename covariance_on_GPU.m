clear; close all;

% dispersion_dataset_filename = "C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\2D-dispersion\OUTPUT\N_pix4x4 N_ele2x2 N_wv101x51 N_disp10000 N_eig20 offset0 output 11-Jun-2021 13-24-45\DATA N_wv101x51 N_disp10000 RNG_offset0 11-Jun-2021 13-24-45.mat";
%
% [WAVEVECTOR_DATA,EIGENVALUE_DATA] = load_dispersion_dataset(dispersion_dataset_filename);
%
% [N_disp,N_eig,pN_wv] = size(EIGENVALUE_DATA);
% eig_idx = 1;
% temp = cov(squeeze(EIGENVALUE_DATA(:,eig_idx,:))');

lb = 10;
ub = 10000;
N = 5;
N_var_iter = round(logspace(log10(lb),log10(ub),N));
lb = 10;
ub = 10000;
M = 5;
N_sample_iter = round(logspace(log10(lb),log10(ub),M));

for i = 1:length(N_sample_iter)
    N_sample = N_sample_iter(i);
    for j = 1:length(N_var_iter)
        N_var = N_var_iter(j);
        dataset = rand(N_sample,N_var);
        dataset_info = whos('dataset');
        dataset_memory(i,j) = dataset_info.bytes;
        t0 = tic;
        C = cov(dataset);
        t = toc(t0);
        cpu_times(i,j) = t;
        F = @() cov(dataset);
        cpu_times_2(i,j) = timeit(F);
        C_info = whos('C');
        C_memory(i,j) = C_info.bytes;
        dataset_gpu = gpuArray(dataset);
        t0 = tic;
        C_gpu = cov(dataset_gpu);
        t = toc(t0);
        gpu_times(i,j) = t;
        F = @() cov(dataset_gpu);
        gpu_times_2(i,j) = gputimeit(F);
    end
end

f = figure;
m = 2; n = 3;
t = tiledlayout(m,n);
for i = 1:m
    for j = 1:n
        ax(i,j) = nexttile;
        %         daspect([1 1 1])
        %         set(ax(i,j),'Visible','off')
    end
end

i = 1; j = 1;
axes(ax(i,j))
im = imagesc([min(N_var_iter) max(N_var_iter)],[min(N_sample_iter) max(N_sample_iter)],dataset_memory/1e9);
apply_common_settings()
title('memory required to store dataset [GB]')
set(gca,'XScale','log')
set(gca,'XTick',N_var_iter)
lv = log10(N_var_iter);
dlv = diff(lv);
mp = mean(dlv);
mpb2 = mp/2;
lb = min(N_var_iter)/(10^mpb2);
ub = max(N_var_iter)*(10^mpb2);
axis([lb ub 0 1000])

i = 2; j = 1;
axes(ax(i,j))
imagesc(N_var_iter,N_sample_iter,C_memory/1e9)
apply_common_settings()
title('memory required to store covariance [GB]')

i = 1; j = 2;
axes(ax(i,j))
imagesc(N_var_iter,N_sample_iter,cpu_times)
apply_common_settings()
title('time taken to compute covariance [s] (CPU)')

i = 2; j = 2;
axes(ax(i,j))
imagesc(N_var_iter,N_sample_iter,gpu_times)
apply_common_settings()
title('time taken to compute covariance [s] (GPU)')

i = 1; j = 3;
axes(ax(i,j))
imagesc(N_var_iter,N_sample_iter,cpu_times_2)
apply_common_settings()
title('time taken to compute covariance [s] (CPU) (timeit)')

i = 2; j = 3;
axes(ax(i,j))
imagesc(N_var_iter,N_sample_iter,gpu_times_2)
apply_common_settings()
title('time taken to compute covariance [s] (GPU) (gputimeit)')

function apply_common_settings()
%     daspect([1 1 1])
    xlabel('number of variables')
    ylabel('number of observations')
    set(gca,'YDir','normal')
    colorbar
    set(gca,'colorscale','log')
%     set(gca,'XScale','log')
end