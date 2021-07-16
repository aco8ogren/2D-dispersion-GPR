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
%     '2D-dispersion-GPR\OUTPUT\error_analysis_data\error_analysis_data_gold_ground_truth3_unfinished.mat']);

datas{1} = load(['C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\'...
    '2D-dispersion-GPR\OUTPUT\error_analysis_data\error_analysis_data_gold_ground_truth3.mat']);

% datas{1} = load(['C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\'...
%     '2D-dispersion-GPR\error_analysis_data.mat']);

% datas{2} = load(['C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\2D-dispersion-GPR\'...
%     'OUTPUT\error_analysis_data\test\error_analysis_data2.mat']);

% In group format
colormap_names = {'cool','autumn'}; % cool, autumn
data_idxs = [1 1]; % indicates which dataset to draw from for the group_idx group. data_idx(group_idx) will be called later.
inner_group_idxs = [1 2];
quantiles_of_interest = 0.95*ones(length(data_idxs),1);
isAverageOverEigenvalue = true; accuracies{1} = [40 10 2]; accuracies{2} = [ 70 50 35 ];
isFindMatchingAccuracy = true;

N_outer_group = length(data_idxs); % Should assert that this is the length of the sig of int and quant of int also and data idxs

legend_name_functions = {...
    @(eig_idx,quantile_of_interest,model_name,sample_resolution,sigma) ...
    [model_name ' w/ \sigma = 1e' num2str(log10(sigma)) ' eig\_idx = ' num2str(eig_idx)],...
    @(eig_idx,quantile_of_interest,model_name,sample_resolution,sigma) ...
    [model_name ' eig\_idx = ' num2str(eig_idx)]...
    };
eig_idx = []; % should these empty initializations be in the for outer_group_idx loop?
quantile_of_interest = [];
model_name = [];
sample_resolution = [];
sigma_of_interest = [];

colors = cell(0,0);
for outer_group_idx = 1:N_outer_group
    data_idx = data_idxs(outer_group_idx);
    inner_group_idx = inner_group_idxs(outer_group_idx);
    colors{outer_group_idx} = get_N_colors(colormap_names{outer_group_idx},length(datas{data_idx}.eig_idxs_iter{inner_group_idx}));
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
            colors_temp = colors{outer_group_idx};
            color = colors_temp(1,:);
            marker = markers{outer_group_idx};
            Q_temp_average = mean(Q_temp,1);
%             p_temp = plot(prod(sample_resolutions_temp,2),Q_temp_average,'color',color,'marker',marker);
            p_temp = plot(sample_counts,Q_temp_average,'color',color,'marker',marker);
            plot_counter = plot_counter + 1;
            p(plot_counter) = p_temp(1);
            legend_name_function = legend_name_functions{outer_group_idx};
            eig_idx = 'average';
            legend_label = legend_name_function(eig_idx,quantile_of_interest,model_name,sample_resolution,sigma);
            p(plot_counter).DisplayName = legend_label;
            if isFindMatchingAccuracy
                log_sample_counts_interp{outer_group_idx,err_idx} = interp1(log10(Q_temp_average),log10(sample_counts),log10(accuracies{err_idx}));
                scatter(10.^(log_sample_counts_interp{outer_group_idx,err_idx}),accuracies{err_idx},'ko')
                for acc_idx = 1:length(accuracies{err_idx})
                    text(10.^(log_sample_counts_interp{outer_group_idx,err_idx}(acc_idx))*1.05,accuracies{err_idx}(acc_idx)*1.05,num2str(10.^(log_sample_counts_interp{outer_group_idx,err_idx}(acc_idx))))
                    yline(accuracies{err_idx}(acc_idx));
                end
            end
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
        set(gca,'xscale','log')
        axis tight
        legend(p,'location',legend_location)
        title([err_names{err_idx} ' error quantiles ' num2str(quantile_of_interest)])
        xlabel('number of sample points')
        ylabel([err_names{err_idx} ' error'])
        grid on
        grid minor
    end
end

