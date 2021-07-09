clear; close all;
% =========================================================================
% MUST BE RUN IN MATLAB R2020a (or maybe a later version would be okay too)
% =========================================================================

warning('off','MATLAB:MKDIR:DirectoryExists')

addpath('../2D-dispersion-optimization-GPR')
addpath('../2D-dispersion')

isPlot = true;
isSaveData = true;

% Iteratable variables (entered in group format)
eig_idxs_iter = {1:5}; % can be 'all'
struct_idxs_iter = {1:3}; % can be 'all'
model_names = {'GPR'}; % can be 'GPR','linear','nearest','cubic','makima','spline'
sample_methods = {'scattered'}; % can be 'gridded', 'scattered'. Must be 'gridded' for any model except 'GPR'.
% sample_resolutions_iter = {[(6:5:51)'; 52], [(6:5:51)'; 52]}; for i = 1:length(sample_resolutions_iter); sample_resolutions_iter{i}(:,2) = ceil(sample_resolutions_iter{i}(:,1)/2); end;
sample_resolutions_iter = {[51]}; for i = 1:length(sample_resolutions_iter); sample_resolutions_iter{i}(:,2) = ceil(sample_resolutions_iter{i}(:,1)/2); end;
for i = 1:length(sample_resolutions_iter); sample_counts_iter{i} = prod(sample_resolutions_iter{i},2) - 2*floor(sample_resolutions_iter{i}(:,1)/2) - sample_resolutions_iter{i}(:,2) + 2; end; % subtraction to account for trimming, trimming to account for symmetry
sigmas = {0};

temp = [length(eig_idxs_iter) length(struct_idxs_iter) length(model_names) length(sample_methods) length(sample_resolutions_iter) length(sample_counts_iter) length(sigmas)];
assert(all(temp(1) == temp));
N_group = temp(1);
clear temp;

% Define training sets ====================================================

for group_idx = 1:N_group
    data_path_train_iter{group_idx} = ['C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\2D-dispersion\OUTPUT\'...
        'gold_dataset output 11-Jun-2021 13-24-45\DATA N_wv101x51 N_disp10000 RNG_offset0 11-Jun-2021 13-24-45.mat'];
%     data_path_train_iter{group_idx} = ['C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\2D-dispersion\OUTPUT\'...
%         'gold_dataset_small output 08-Jul-2021 12-54-14\DATA N_disp100 N_wv101x51 RNG_offset0 08-Jul-2021 12-54-14.mat'];
end

% Define test sets ========================================================

for group_idx = 1:N_group
%     data_path_test_iter{group_idx} = ['C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\2D-dispersion\OUTPUT\'...
%         'ground_truth3 output 28-May-2021 16-41-37\DATA N_struct100 N_k RNG_offset0 28-May-2021 16-41-37.mat'];
%     data_path_test_iter{group_idx} = ['C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\2D-dispersion\OUTPUT\'...
%         'gold_dataset_small output 08-Jul-2021 12-54-14\DATA N_disp100 N_wv101x51 RNG_offset0 08-Jul-2021 12-54-14.mat'];
    data_path_test_iter{group_idx} = ['C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\2D-dispersion\OUTPUT\'...
        'ground_truth4 output 08-Jul-2021 16-18-23\DATA N_disp100 N_wv501x251 RNG_offset0 08-Jul-2021 16-18-23.mat'];
end

% Define sample order =====================================================

for group_idx = 1:N_group
    sample_order_path_iter{group_idx} = ['C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\2D-dispersion-GPR\OUTPUT\sample_orders\'...
        'sample_order_data_gold1.mat'];
end

% =========================================================================

% Datasets ================================================================

% data_path = ['C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\2D-dispersion\OUTPUT\'...
%     'light_dataset output 01-Jul-2021 17-32-47\DATA N_disp1000 N_wv51x26 01-Jul-2021 17-32-47.mat'];
%
% data_path = ['C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\2D-dispersion\OUTPUT\'...
%     'ground_truth3 output 28-May-2021 16-41-37\DATA N_struct100 N_k RNG_offset0 28-May-2021 16-41-37.mat'];
%
% data_path = ['C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\2D-dispersion\OUTPUT\'...
%     'ground_truth output 20-May-2021 17-11-07\DATA N_struct128 N_k201 RNG_offset0 20-May-2021 17-11-07.mat'];
%
% data_path = ['C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\2D-dispersion\OUTPUT\'...
%     'covariance_singularity N_wv101x51 N_disp10000 output 11-Jun-2021 13-24-45\DATA N_struct10000 N_k RNG_offset0 11-Jun-2021 13-24-45.mat'];
%
% data_path = ['C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\2D-dispersion\OUTPUT\'...
%     'covariance_singularity N_wv51x26 N_disp5000 output 04-Jun-2021 16-31-13\DATA N_struct5000 N_k RNG_offset0 04-Jun-2021 16-31-13.mat'];
%
% data_path = ['C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\2D-dispersion\OUTPUT\'...
%     'covariance_singularity N_wv51x26 N_disp20000 output 10-Jun-2021 14-58-54\DATA N_struct20000 N_k RNG_offset0 10-Jun-2021 14-58-54.mat'];

% data_path__iter{group_idx} = ['C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\2D-dispersion\OUTPUT\'...
%     'covariance_singularity N_wv101x51 N_disp10000 output 11-Jun-2021 13-24-45\DATA N_struct10000 N_k RNG_offset0 11-Jun-2021 13-24-45.mat'];

% Training sets ===========================================================

% data_path_train = ['C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\'...
%     '2D-dispersion\OUTPUT\ground_truth output 20-May-2021 17-11-07\DATA N_struct128 N_k201 RNG_offset0 20-May-2021 17-11-07.mat'];

% data_path_train = ['C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\'...
%     '2D-dispersion\OUTPUT\covariance_singularity N_wv51x26 N_disp5000 output 04-Jun-2021 16-31-13\DATA N_struct5000 N_k RNG_offset0 04-Jun-2021 16-31-13.mat'];

% data_path_train = ['C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\'...
%     '2D-dispersion\OUTPUT\covariance_singularity N_wv51x26 N_disp20000 output 10-Jun-2021 14-58-54\DATA N_struct20000 N_k RNG_offset0 10-Jun-2021 14-58-54.mat'];

% data_path_train = ['C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\'...
%     '2D-dispersion\OUTPUT\covariance_singularity N_wv101x51 N_disp10000 output 11-Jun-2021 13-24-45\DATA N_struct10000 N_k RNG_offset0 11-Jun-2021 13-24-45.mat'];

% data_path_train = ['C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\'...
%     '2D-dispersion\OUTPUT\light_dataset output 01-Jul-2021 17-32-47\DATA N_disp1000 N_wv51x26 01-Jul-2021 17-32-47.mat'];

% test sets ===============================================================

% data_path_test = ['C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\'...
%     '2D-dispersion\OUTPUT\ground_truth3 output 28-May-2021 16-41-37\DATA N_struct100 N_k RNG_offset0 28-May-2021 16-41-37.mat'];

% data_path_test = ['C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\'...
%     '2D-dispersion\OUTPUT\ground_truth output 20-May-2021 17-11-07\DATA N_struct128 N_k201 RNG_offset0 20-May-2021 17-11-07.mat'];

% data_path_test = ['C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\'...
%     '2D-dispersion-GPR\OUTPUT\Homog w dataset N_k51\DATA N_struct128 N_k51 RNG_offset0 14-Mar-2021 16-46-17.mat'];

% data_path_test = ['C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\'...
%     '2D-dispersion\OUTPUT\covariance_singularity N_wv101x51 N_disp10000 output 11-Jun-2021 13-24-45\DATA N_struct10000 N_k RNG_offset0 11-Jun-2021 13-24-45.mat'];

% data_path_test = ['C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\'...
%     '2D-dispersion\OUTPUT\covariance_singularity N_wv51x26 N_disp5000 output 04-Jun-2021 16-31-13\DATA N_struct5000 N_k RNG_offset0 04-Jun-2021 16-31-13.mat'];

% data_path_test = ['C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\'...
%     '2D-dispersion\OUTPUT\light_dataset output 01-Jul-2021 17-32-47\DATA N_disp1000 N_wv51x26 01-Jul-2021 17-32-47.mat'];

% =========================================================================

% Sample orders ===========================================================

% sample_order_path_iter{group_idx} = ['C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\2D-dispersion-GPR\OUTPUT\sample_orders\'...
%     'sampling_order_information_LD_1326_1.mat'];

% sample_order_path_iter{group_idx} = ['C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\2D-dispersion-GPR\OUTPUT\sample_orders\'...
%     'sample_order_data_gold1.mat'];

% =========================================================================

for group_idx = 1:N_group
    disp(repmat('=',1,40))
    disp(['Processing group ' num2str(group_idx)])
    disp(repmat('-',1,40))
    sample_order_path = sample_order_path_iter{group_idx};
    data_path_train = data_path_train_iter{group_idx};
    data_path_test = data_path_test_iter{group_idx};
    
    sigma = sigmas{group_idx};
    sample_method = sample_methods{group_idx};
    model_name = model_names{group_idx};
    eig_idxs = eig_idxs_iter{group_idx};
    struct_idxs = struct_idxs_iter{group_idx};
    sample_resolutions = sample_resolutions_iter{group_idx};
    sample_counts = sample_counts_iter{group_idx};
    
    disp('Loading sampling_order_data...')
    sample_order_data = load(sample_order_path);
    disp('done.')
    
    disp('Loading training set...')
    [WAVEVECTOR_DATA_TRAIN,EIGENVALUE_DATA_TRAIN] = load_dispersion_dataset(data_path_train);
    disp('done.')
    
    disp('Loading test set...')
    [WAVEVECTOR_DATA_TEST,EIGENVALUE_DATA_TEST] = load_dispersion_dataset(data_path_test);
    disp('done.')
    
    disp('Applying models and evaluating errors...')
    
    regexp_idx = regexp(data_path_train,'\');
    data_dir = data_path_train(1:(regexp_idx(end)));
    script_start_time = replace(char(datetime),':','-');
    
    N_evaluate = [numel(unique(WAVEVECTOR_DATA_TEST(:,1,1))) numel(unique(WAVEVECTOR_DATA_TEST(:,2,1)))];
    
    [N_wv_test,N_eig_test,N_struct_test] = size(EIGENVALUE_DATA_TEST);
    
    if strcmp(struct_idxs,'all')
        struct_idxs = 1:N_struct_test;
    end
    
    if strcmp(eig_idxs,'all')
        eig_idxs = 1:N_eig_test;
        eig_idxs_iter{group_idx} = eig_idxs;
    end
    
    err_L2{group_idx} = zeros(length(eig_idxs_iter{group_idx}),length(struct_idxs_iter{group_idx}),size(sample_resolutions_iter{group_idx},1));
    err_H1{group_idx} = zeros(length(eig_idxs_iter{group_idx}),length(struct_idxs_iter{group_idx}),size(sample_resolutions_iter{group_idx},1));
    
    covariance_options.isAllowGPU = false;
    covariance_options.isComputeCovarianceGradient = false;
    
    for eig_idx_idx = 1:length(eig_idxs)
        eig_idx = eig_idxs(eig_idx_idx);
        model_options = struct();
        if ismember('GPR',model_names)
            original_wv_x = unique(sort(WAVEVECTOR_DATA_TRAIN(:,1,1)));
            original_wv_y = unique(sort(WAVEVECTOR_DATA_TRAIN(:,2,1)));
            
            covariance_options.eig_idxs = eig_idx;
            [Cs,C_grads,kfcns,kfcn_grads] = get_empirical_covariance(WAVEVECTOR_DATA_TRAIN,EIGENVALUE_DATA_TRAIN,covariance_options); %#ok<ASGLU>
            model_options.kfcn = kfcns{1};
            model_options.kfcn_grad = {};
        end
        
        wb_counter = 0;
        wb = waitbar(0,['Processing eig\_idx = ' num2str(eig_idx)]);
        for sample_idx = 1:size(sample_resolutions,1)
            N_sample = sample_resolutions(sample_idx,:);
            for struct_idx_idx = 1:length(struct_idxs)
                struct_idx = struct_idxs(struct_idx_idx);
                fr = squeeze(EIGENVALUE_DATA_TEST(:,eig_idx,struct_idx));
                wv = squeeze(WAVEVECTOR_DATA_TEST(:,:,struct_idx));
                
                %                 for model_idx = 1:length(model_names)
                %                     model_name = model_names{model_idx};
                model_options.model_name = model_name;
                a = 1;
                if strcmp(model_name,'GPR')
                    model_options.kfcn = kfcns{1};
                    model_options.sample_interpolation_format = 'scattered'; % specify this because these scattered points cannot be interpolated as a grid
                    model_options.train_format = 'scattered'; % only applies to GPR. must be 'scattered' even for 'gridded' sampling since the trim operation no longer leaves a gridded domain
                    model_options.predict_format = 'scattered'; % it has to be 'scattered' if train_format is 'scattered'.
                    if strcmp(sample_method,'scattered')
                        model_options.sample_points = sample_order_data.sample_orders(1:sample_counts(sample_idx),:,eig_idx); % get sample points from loaded file
                    elseif strcmp(sample_method,'gridded')
                        model_options.sample_points = get_wavevectors(N_sample,a,struct('isTrimRightBoundary',true,'format','list'));
                    else
                        error('sampling_method not recognized')
                    end
                else
                    
                    model_options.sample_points = get_wavevectors(N_sample,a,struct('isTrimRightBoundary',false,'format','list')); % get them as a list for uniformity.
                    model_options.sample_interpolation_format = 'gridded'; % specify this for speed
                end
                
                %                     for sigma_idx = 1:length(sigmas)
                %                         sigma = sigmas(sigma_idx);
                model_options.sigma = sigma;
                out = interp_model_2D(fr,wv,model_options);
                err_L2{group_idx}(eig_idx_idx,struct_idx_idx,sample_idx) = out.e_L2;
                err_H1{group_idx}(eig_idx_idx,struct_idx_idx,sample_idx) = out.e_H1;
                %                     end
            end
            %             end
            wb_counter = wb_counter + 1;
            waitbar(wb_counter/size(sample_resolutions,1),wb)
        end
        close(wb)
    end
    disp('done.')
    disp(repmat('=',1,40))
end

if isSaveData
    disp('Computing quantiles and saving data...')
    q = linspace(0,1,101);
    for group_idx = 1:N_group
    Q_L2{group_idx} = quantile(err_L2{group_idx},q,2);
    Q_H1{group_idx} = quantile(err_H1{group_idx},q,2);
    end
    save('error_analysis_data',...
        'q','Q_L2','Q_H1','err_L2','err_H1',...
        'eig_idxs_iter','struct_idxs_iter','model_names','sample_resolutions_iter','sample_counts_iter','sigmas',...
        'data_path_train_iter','data_path_test_iter')
    disp('done.')
end

error('the end of this script doesn''t work, but don''t worry - the data was saved')

% Plot the qth error quantile for each model & eigensurface
q = .95;
Q_L2_temp = quantile(err_L2,q,2);
Q_H1_temp = quantile(err_H1,q,2);
cm = lines(length(eig_idxs));
mm = {'s','d','^','x','*','o'};
if isPlot
    clear p
    for err_idx = 1:2
        if err_idx == 1
            Q = Q_L2_temp;
            err_name = 'L2';
        elseif err_idx == 2
            Q = Q_H1_temp;
            err_name = 'H1';
        end
        for sigma_idx = 1:length(sigmas)
            sigma = sigmas(sigma_idx);
            
            for model_idx = 1:length(model_names)
                figure
                hold on
                for eig_idx_idx = 1:length(eig_idxs)
                    eig_idx = eig_idxs(eig_idx_idx);
                    p_temp = plot(prod(sample_resolutions,2),squeeze(Q(eig_idx_idx,:,model_idx,:,sigma_idx)),'color',cm(eig_idx_idx,:),'marker',mm{model_idx});
                    p(eig_idx_idx) = p_temp(1);
                    p(eig_idx_idx).DisplayName = [model_names{model_idx} ' eig\_idx = ' num2str(eig_idx)];
                end
                set(gca,'yscale','log')
                legend(reshape(p,[],1),'location','northeast')
                title([err_name ' error quantiles ' regexprep(num2str(q),' +',',') newline model_names{model_idx} ' \sigma = 1e' num2str(log10(sigma))])
                xlabel('prod(N\_sample)')
                ylabel([err_name ' error'])
                grid on
                grid minor
            end
        end
    end
end

