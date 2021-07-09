clear; close all;

err_names = {'L2','H1'};

Q = containers.Map;
sample_resolutions = containers.Map;
eig_idxs = containers.Map;


% load GPR quantile data
% data = load('C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\2D-dispersion-GPR\OUTPUT\ground_truth1 quantile data\quantile_data_GPR.mat');
% data = load('C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\2D-dispersion-GPR\OUTPUT\ground_truth3 quantile data\quantile_data_GPR_confusing.mat');
% data = load('C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\2D-dispersion-GPR\OUTPUT\ground_truth3 quantile data\quantile_data_GPR_ttt.mat');

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

datas{1} = load(['C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\'...
    '2D-dispersion-GPR\OUTPUT\error_analysis_data\error_analysis_data_gold1.mat']);

% datas{1} = load(['C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\'...
%     '2D-dispersion-GPR\error_analysis_data.mat']);

% datas{2} = load(['C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\2D-dispersion-GPR\'...
%     'OUTPUT\error_analysis_data\test\error_analysis_data2.mat']);

% In group format
% sigmas_of_interest = {1e-2, NaN};
% model_names = {'GPR','linear'};
colormap_names = {'cool','autumn'}; % cool, autumn
data_idxs = [1 1]; % indicates which dataset to draw from for the group_idx group. data_idx(group_idx) will be called later.
inner_group_idxs = [1 2];
quantiles_of_interest = 0.95*ones(length(data_idxs),1);

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
%     sigma_of_interest = sigmas_of_interest{outer_group_idx};
    
    data = datas{data_idx};
    
    inner_group_idx = inner_group_idxs(outer_group_idx);
    
    quantile_idx = find(data.q == quantile_of_interest);
    
    sigma = data.sigmas{inner_group_idx};
    
    % sigma_idx no longer exists.
%     if ~isnan(sigma_of_interest)
%     sigma_idx = find([data.sigmas{:}] == sigma_of_interest);
%     else
%         sigma_idx = 1;
%     end
%     model_idx = find(strcmp(data.model_names,model_names{outer_group_idx}));
    assert(~isempty(quantile_idx))
    
    % sigma_idx and model_idx are no longer a thing
%     assert(~isempty(sigma_idx))
%     assert(~isempty(model_idx))
    
    Q('L2') = squeeze(data.Q_L2{inner_group_idx}(:,quantile_idx,:));
    Q('H1') = squeeze(data.Q_H1{inner_group_idx}(:,quantile_idx,:));
    

        
%         for model_idx = 1:length(model_names)
%             model_name = model_names{model_idx};
        model_name = data.model_names{inner_group_idx};
%             eig_idxs_temp = data.eig_idxs;
        eig_idxs_temp = data.eig_idxs_iter{inner_group_idx};
            for eig_idx_idx = 1:length(eig_idxs_temp)
                eig_idx = eig_idxs_temp(eig_idx_idx);
%                 Q_temp = Q([model_name ' ' err_name]);
                Q_temp = Q(err_name);
%                 sample_resolutions_temp = data.sample_resolutions;
                sample_resolutions_temp = data.sample_resolutions_iter{inner_group_idx};
                colors_temp = colors{outer_group_idx};
                color = colors_temp(eig_idx_idx,:);
                marker = markers{outer_group_idx};
                p_temp = plot(prod(sample_resolutions_temp,2),Q_temp(eig_idx_idx,:),'color',color,'marker',marker);
                plot_counter = plot_counter + 1;
                p(plot_counter) = p_temp(1);
                legend_name_function = legend_name_functions{outer_group_idx};
                legend_label = legend_name_function(eig_idx,quantile_of_interest,model_name,sample_resolution,sigma);
                p(plot_counter).DisplayName = legend_label;   
            end
%         end
        set(gca,'yscale','log')
        axis tight
        legend(p,'location',legend_location)
        title([err_names{err_idx} ' error quantiles ' num2str(quantile_of_interest)])
        xlabel('number of sample points')
        ylabel([err_names{err_idx} ' error'])
        grid on
        grid minor
    end
end

% % make plots of smoothed quantile data
% for key_cell = Q.keys
%     key = char(key_cell);
%     Q([key ' smooth']) = smoothdata(Q(key),2,'movmedian',10);
% end
% 
% for err_idx = 1:2
%     figure
%     hold on
%     for model_idx = 1:2
%         for eig_idx = 1:min(cell2mat(values(eig_idxs)))
%             model_name = model_names{model_idx};
%             Q_temp = Q([model_name ' ' err_names{err_idx} ' smooth']);
%             colors_temp = colors(model_name);
%             p_temp = plot(sample_resolutions(model_name),Q_temp(eig_idx,:),'color',colors_temp(eig_idx,:),'marker',markers{model_idx});
%             p(eig_idx,model_idx) = p_temp(1);
%             p(eig_idx,model_idx).DisplayName = [model_name ' eig\_idx = ' num2str(eig_idx)];
%         end
%     end
%     set(gca,'yscale','log')
%     axis tight
%     legend(reshape(p,[],1),'location',legend_location)
%     title([err_names{err_idx} ' error quantiles ' num2str(quantile_of_interest) newline 'smoothed'])
%     xlabel('N\_sample')
%     ylabel([err_names{err_idx} ' error'])
%     grid on
%     grid minor
% end
% 
% % make plots of distances to GPR errors at specific GPR sample resolutions
% N_samples_of_interest = [7 13 22];
% 
% for err_idx = 1:2
%     for N_samples_of_interest_idx = 1:length(N_samples_of_interest)
%         clear p
%         figure
%         hold on
%         for eig_idx_idx = 1:min(cell2mat(values(eig_idxs)))
%             eig_idx = eig_idxs(eig_idx_idx);
%             N_sample_GPR = N_samples_of_interest(N_samples_of_interest_idx);
%             Q_temp_GPR = Q(['GPR ' err_names{err_idx} ' smooth']);
%             idx = find(sample_resolutions('GPR') == N_sample_GPR);
%             Q_temp_GPR = Q_temp_GPR(:,idx);
%             Q_temp_linear = Q(['linear ' err_names{err_idx} ' smooth']);
%             distances = abs(Q_temp_linear(eig_idx_idx,:) - Q_temp_GPR(eig_idx_idx));
%             N_samp_lin = sample_resolutions('linear');
%             colors_temp = colors('linear');
%             p(eig_idx_idx) = plot(N_samp_lin,distances,'color',colors_temp(eig_idx_idx,:),'marker',markers{find(strcmp(model_names,'linear'))});
%             p(eig_idx_idx).DisplayName = ['eig\_idx = ' num2str(eig_idx)];
%             [pks,locs] = findpeaks(-log10(distances),'SortStr','descend');
%             
%             N_samp_lin_match(eig_idx_idx) = N_samp_lin(locs(1));
%             scatter(N_samp_lin(locs(1)),10^(-pks(1)),'ko');
%         end
%         set(gca,'yscale','log')
%         axis tight
%         legend(reshape(p,[],1),'location',legend_location)
%         title(['Distance of linear ' err_names{err_idx} ' error quantile ' num2str(quantile_of_interest) ' to GPR error quantile ' num2str(quantile_of_interest)...
%             newline 'GPR N\_sample = ' num2str(N_sample_GPR)...
%             newline 'Matching linear N\_sample = ' regexprep(num2str(N_samp_lin_match),' +',' ')])
%         xlabel('N\_sample for linear model')
%         ylabel(['Distance from linear ' err_names{err_idx} ' error quantile' newline 'to fixed GPR error quantile'])
%         grid on
%         grid minor
%     end
% end