clear; close all;

err_names = {'L2','H1'};

Q = containers.Map;
sample_resolutions = containers.Map;

% datas{1} = load(['C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\2D-dispersion-GPR\'...
%     'OUTPUT\error_analysis_data\test\error_analysis_data1.mat']);

% datas{1} = load(['C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\2D-dispersion-GPR\'...
%     'OUTPUT\error_analysis_data\error_analysis_data_up_to_128_out_of_sample.mat']);

% datas{1} = load(['C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\2D-dispersion-GPR\'...
%     'OUTPUT\error_analysis_data\error_analysis_data_up_to_1326_out_of_sample.mat']);

% datas{1} = load(['C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\2D-dispersion-GPR\'...
%     'OUTPUT\error_analysis_data\error_analysis_data_up_to_1326_orig_test_set.mat']);

% datas{1} = load(['C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\2D-dispersion-GPR\'...
%     'OUTPUT\error_analysis_data\error_analysis_data_up_to_128_in_sample.mat']);

% datas{1} = load(['C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\2D-dispersion-GPR\'...
%     'OUTPUT\error_analysis_data\error_analysis_data_up_to_1326_in_sample.mat']);

% datas{1} = load(['C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\2D-dispersion-GPR\'...
%     'OUTPUT\error_analysis_data\error_analysis_data_newest.mat']);

% datas{1} = load(['C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\2D-dispersion-GPR\'...
%     'OUTPUT\error_analysis_data\error_analysis_data_newest_sigma0.mat']);

% datas{1} = load(['C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\'...
%     '2D-dispersion-GPR\OUTPUT\error_analysis_data\error_analysis_data_LD_1326_1.mat']);

% datas{1} = load(['C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\'...
%     '2D-dispersion-GPR\OUTPUT\error_analysis_data\error_analysis_data_gold1.mat']);

% datas{1} = load(['C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\'...
%     '2D-dispersion-GPR\error_analysis_data.mat']);

% datas{2} = load(['C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\2D-dispersion-GPR\'...
%     'OUTPUT\error_analysis_data\test\error_analysis_data2.mat']);

% datas{1} = load(['C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\'...
%     '2D-dispersion-GPR\OUTPUT\error_analysis_data\error_analysis_data_gold4x4_ground_truth3.mat']);
% % 
% datas{2} = load(['C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\'...
%     '2D-dispersion-GPR\OUTPUT\error_analysis_data\error_analysis_data_gold8x8_ground_truth3_partial.mat']);
% % 
% datas{3} = load(['C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\'...
%     '2D-dispersion-GPR\OUTPUT\error_analysis_data\error_analysis_data_linear_ground_truth3.mat']);
% 
% datas{4} = load(['C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\2D-dispersion-GPR\OUTPUT\error_analysis_data\'...
%     'error_analysis_data_gold4x4_ground_truth5_partial.mat']);
% 
% datas{1} = load(['C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\2D-dispersion-GPR\OUTPUT\error_analysis_data\'...
%     'error_analysis_data_light_dataset_1326_1e-16.mat']);
% 
% datas{2} = load(['C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\2D-dispersion-GPR\OUTPUT\error_analysis_data\'...
%     'error_analysis_data_linear_light_dataset2.mat']);

% datas{1} = load(['C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\2D-dispersion-GPR\OUTPUT\error_analysis_data\'...
%     'error_analysis_data_ultralight_dataset_test_is_train.mat']);

% datas{1} = load(['C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\2D-dispersion-GPR\OUTPUT\error_analysis_data\'...
%     'error_analysis_data_superultralight_dataset_test_is_train.mat']);

% datas{1} = load(['C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\2D-dispersion-GPR\OUTPUT\error_analysis_data\'...
%     'error_analysis_data_superultralight_dataset_out_of_sample.mat']);

% datas{1} = load(['C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\2D-dispersion-GPR\OUTPUT\error_analysis_data\'...
%     'error_analysis_data_sul_sulhr_linear.mat']);

% datas{1} = load(['C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\2D-dispersion-GPR\OUTPUT\error_analysis_data\'...
%     'error_analysis_data_sul_sulhr_cubic.mat']);

% datas{4} = load(['C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\2D-dispersion-GPR\OUTPUT\error_analysis_data\'...
%     'error_analysis_data_gold8x8_ground_truth3_cubic_covariance_interp.mat']);

% datas{1} = load(['C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\'...
%     '2D-dispersion-GPR\error_analysis_data.mat']);

% datas{1} = load(['C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\2D-dispersion-GPR\OUTPUT\error_analysis_data\'...
%     'error_analysis_data_train4x4_test4x4_N_wv31x16.mat']);

% datas{2} = load(['C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\2D-dispersion-GPR\OUTPUT\error_analysis_data\'...
%     'error_analysis_data_N_pix4x4_N_wv31x16_linear.mat']);

% datas{1} = load(['C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\2D-dispersion-GPR\OUTPUT\'...
%     'error_analysis_data\error_analysis_data_N_pix4x4_N_wv101x51.mat']);

% datas{2} = load(['C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\2D-dispersion-GPR\OUTPUT\'...
%     'error_analysis_data\error_analysis_data_N_pix4x4_N_wv101x51_linear.mat']);

% datas{1} = load(['C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\2D-dispersion-GPR\OUTPUT\'...
%     'error_analysis_data\error_analysis_data_N_pix4x4_N_wv151x76.mat']);
% 
% datas{2} = load(['C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\2D-dispersion-GPR\OUTPUT\'...
%     'error_analysis_data\error_analysis_data_N_pix4x4_N_wv151x76_linear.mat']);

datas{1} = load("C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\2D-dispersion-GPR\OUTPUT\error_analysis_data\error_analysis_data_new_ID_ss.mat");

% In group format
% sigmas_of_interest = {1e-2, NaN};
% model_names = {'GPR','linear'};
colormap_names = {'cool'}; % cool, autumn
data_idxs = [1 1 1 1 1 1 1 1 1]; % indicates which dataset to draw from for the group_idx group. data_idx(group_idx) will be called later.
inner_group_idxs = [1 2 3 4 5 6 7 8 9];
quantiles_of_interest = 0.95*ones(length(data_idxs),1);
isAverageOverEigenvalue = true;
isXScaleLog = true;

N_outer_group = length(data_idxs); % Should assert that this is the length of the sig of int and quant of int also and data idxs

group_names = {'50','100','250','500','1000','2000','3000','4000','5000'};
% group_names = {'GPR (train 4x4 test 8x8)','GPR (train 8x8 test 8x8)','linear',['GPR (train 8x8 test 8x8),' newline 'cubic covariance interpolation']};
% group_names = {'GPR, train = SUL dataset, test = SUL 2 dataset'};
% group_names = {'GPR (train 4x4 test 4x4 same set)','linear'};

for group_idx = 1:N_outer_group
    legend_name_functions{group_idx} = @(eig_idx,quantile_of_interest,model_name,group_name,sample_resolution,sigma) group_name;
end

% legend_name_functions = {...
%     @(eig_idx,quantile_of_interest,model_name,sample_resolution,sigma) ...
%     [model_name ' w/ \sigma = 1e' num2str(log10(sigma)) ' eig\_idx = ' num2str(eig_idx)],...
%     @(eig_idx,quantile_of_interest,model_name,sample_resolution,sigma) ...
%     [model_name ' eig\_idx = ' num2str(eig_idx)]...
%     };

eig_idx = []; % should these empty initializations be in the for outer_group_idx loop?
quantile_of_interest = [];
model_name = [];
sample_resolution = [];
sigma_of_interest = [];

colors = cell(0,0);
for outer_group_idx = 1:N_outer_group
    data_idx = data_idxs(outer_group_idx);
    inner_group_idx = inner_group_idxs(outer_group_idx);
    if ~isAverageOverEigenvalue
        colors{outer_group_idx} = get_N_colors(colormap_names{outer_group_idx},length(datas{data_idx}.eig_idxs_iter{inner_group_idx}));
    elseif isAverageOverEigenvalue
        colors = get_N_colors(colormap_names{1},N_outer_group);
    end
end

markers = {'none','none'};

legend_location = 'northeastoutside';

% make plots of original quantile data

figs(1) = figure;
hold on
figs(2) = figure;
hold on


for err_idx = 1:2
    err_name = err_names{err_idx};
    figure(figs(err_idx));
    plot_counter = 0;
    clear p
    for outer_group_idx = 1:N_outer_group
        group_name = group_names{outer_group_idx};
        
        data_idx = data_idxs(outer_group_idx);
        inner_group_idx = inner_group_idxs(outer_group_idx);
        quantile_of_interest = quantiles_of_interest(outer_group_idx);

        data = datas{data_idx};
        
        inner_group_idx = inner_group_idxs(outer_group_idx);
        
        quantile_idx = find(data.q == quantile_of_interest);
        
        sigma = data.sigmas{inner_group_idx};
        
        assert(~isempty(quantile_idx))
        
        Q('L2') = squeeze(data.Q_L2{inner_group_idx}(:,quantile_idx,:));
        Q('H1') = squeeze(data.Q_H1{inner_group_idx}(:,quantile_idx,:));
        
        model_name = data.model_names{inner_group_idx};
        eig_idxs = data.eig_idxs_iter{inner_group_idx};
        if isAverageOverEigenvalue
            Q_temp = Q(err_name);
%             sample_resolutions_temp = data.sample_resolutions_iter{inner_group_idx};
            sample_counts = data.sample_counts_iter{inner_group_idx};
%             colors_temp = colors{outer_group_idx};
            color = colors(outer_group_idx,:);
%             marker = markers{outer_group_idx};
            marker = 'none';
            Q_temp_average = mean(Q_temp,1);
%             p_temp = plot(prod(sample_resolutions_temp,2),Q_temp_average,'color',color,'marker',marker);
            p_temp = plot(sample_counts,Q_temp_average,'color',color,'marker',marker);
            plot_counter = plot_counter + 1;
            p(plot_counter) = p_temp(1);
            legend_name_function = legend_name_functions{outer_group_idx};
            eig_idx = 'average';
            legend_label = legend_name_function(eig_idx,quantile_of_interest,model_name,group_name,sample_resolution,sigma);
            p(plot_counter).DisplayName = legend_label;
        else
            for eig_idx_idx = 1:length(eig_idxs)
                eig_idx = eig_idxs(eig_idx_idx);
                Q_temp = Q(err_name);
%                 sample_resolutions_temp = data.sample_resolutions_iter{inner_group_idx};
                sample_counts = data.sample_counts_iter{inner_group_idx};
                colors_temp = colors{outer_group_idx};
                color = colors_temp(eig_idx_idx,:);
                marker = markers{outer_group_idx};
%                 p_temp = plot(prod(sample_resolutions_temp,2),Q_temp(eig_idx_idx,:),'color',color,'marker',marker);
                p_temp = plot(sample_counts,Q_temp(eig_idx_idx,:),'color',color,'marker',marker);
                plot_counter = plot_counter + 1;
                p(plot_counter) = p_temp(1);
                legend_name_function = legend_name_functions{outer_group_idx};
                legend_label = legend_name_function(eig_idx,quantile_of_interest,model_name,sample_resolution,sigma);
                p(plot_counter).DisplayName = legend_label;
            end
        end
        set(gca,'yscale','log')
        if isXScaleLog
            set(gca,'xscale','log')
        end
        axis tight
        legend(p,'location',legend_location)
        title([err_names{err_idx} ' error quantiles ' num2str(quantile_of_interest)])
        xlabel('number of sample points')
        ylabel([err_names{err_idx} ' error'])
        grid on
        grid minor
    end
end

