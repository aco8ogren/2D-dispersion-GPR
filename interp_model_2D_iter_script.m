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
eig_idxs_iter = repmat({'all'},9,1); % can be 'all'
struct_idxs_iter = repmat({'all'},9,1); % can be 'all'
model_names = repmat({'GPR'},9,1); % can be 'GPR','linear','nearest','cubic','makima','spline'
sample_methods = repmat({'scattered'},9,1); % can be 'gridded', 'scattered'. Must be 'gridded' for any model except 'GPR'.
% sample_resolutions_iter = {[6:5:151]'}; for i = 1:length(sample_resolutions_iter); sample_resolutions_iter{i}(:,2) = ceil(sample_resolutions_iter{i}(:,1)/2); end;
sample_resolutions_iter = repmat({[3:2:9]'},9,1); for i = 1:length(sample_resolutions_iter); sample_resolutions_iter{i}(:,2) = ceil(sample_resolutions_iter{i}(:,1)/2); end;
for i = 1:length(sample_resolutions_iter); sample_counts_iter{i} = prod(sample_resolutions_iter{i},2) - 2*floor(sample_resolutions_iter{i}(:,1)/2) - sample_resolutions_iter{i}(:,2) + 2; end; % subtraction to account for trimming, trimming to account for symmetry
% sample_counts_iter = {[4 34 79 164 254 394 529 724 904 1154 1181]'}; sample_resolutions_iter = {[]};
sigmas = {0,0,0,0,0,0,0,0,0};
N_disps = {50,100,250,500,1000,2000,3000,4000,5000};

temp = [length(eig_idxs_iter) length(struct_idxs_iter) length(model_names) length(sample_methods) length(sample_resolutions_iter) length(sample_counts_iter) length(sigmas)];
assert(all(temp(1) == temp));
N_group = temp(1);
clear temp;

% Define training sets ====================================================

for group_idx = 1:N_group
    %     data_path_train_iter{group_idx} = ['C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\2D-dispersion\OUTPUT\'...
    %         'gold4x4_dataset output 11-Jun-2021 13-24-45\DATA N_wv101x51 N_disp10000 RNG_offset0 11-Jun-2021 13-24-45.mat'];
    %     data_path_train_iter{group_idx} = ['C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\2D-dispersion\OUTPUT\'...
    %         'gold4x4_dataset output 11-Jun-2021 13-24-45\DATA N_wv101x51 N_disp9900 RNG_offset100 11-Jun-2021 13-24-45.mat'];
    %     data_path_train_iter{group_idx} = ['C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\2D-dispersion\OUTPUT\'...
    %         'gold4x4_dataset_small output 08-Jul-2021 12-54-14\DATA N_disp100 N_wv101x51 RNG_offset0 08-Jul-2021 12-54-14.mat'];
    %     data_path_train_iter{group_idx} = ['C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\2D-dispersion\OUTPUT\'...
    %         'gold8x8_dataset output 11-Jul-2021 19-02-27\DATA N_disp9900 N_wv101x51 RNG_offset100 11-Jul-2021 19-02-27.mat'];
    %     data_path_train_iter{group_idx} = ['C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\'...
    %         '2D-dispersion\OUTPUT\light_dataset output 01-Jul-2021 17-32-47\DATA N_disp1000 N_wv51x26 01-Jul-2021 17-32-47.mat'];
    %     data_path_train_iter{group_idx} = ['C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\2D-dispersion\OUTPUT\'...
    %         'ultralight_dataset output 16-Jul-2021 10-39-56\DATA N_disp100 N_wv51x26 RNG_offset0 16-Jul-2021 10-39-56.mat'];
    %     data_path_train_iter{group_idx} = ['C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\2D-dispersion\OUTPUT\'...
    %         'superultralight_dataset1 output 16-Jul-2021 13-46-13\DATA N_disp100 N_wv13x7 RNG_offset0 16-Jul-2021 13-46-13.mat'];
    %     data_path_train_iter{group_idx} = ['C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\2D-dispersion\OUTPUT\'...
    %         'N_pix4x4 N_ele2x2 N_wv31x16 N_disp1000 N_eig10 offset0 output 02-Oct-2021 01-07-49\DATA 02-Oct-2021 01-07-49.mat'];
%     data_path_train_iter{group_idx} = ['C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\2D-dispersion\OUTPUT\'...
%         'N_pix4x4 N_ele2x2 N_wv101x51 N_disp10000 N_eig20 offset0 output 11-Jun-2021 13-24-45\DATA N_wv101x51 N_disp9900 RNG_offset100 11-Jun-2021 13-24-45.mat'];
%     data_path_train_iter{group_idx} = ['C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\2D-dispersion\OUTPUT\'...
%         'N_pix4x4 N_ele2x2 N_wv151x76 N_disp2000 offset100 output 04-Oct-2021 11-27-08\DATA 04-Oct-2021 11-27-08.mat'];
    data_path_train_iter{group_idx} = "C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\2D-dispersion\OUTPUT\new train output 17-Mar-2022 14-21-36\DATA N_pix16x16 N_ele1x1 N_wv31x16 N_disp5000 N_eig5 offset0 17-Mar-2022 14-21-36.mat";
end

% Define test sets ========================================================

for group_idx = 1:N_group
    %     data_path_test_iter{group_idx} = ['C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\2D-dispersion\OUTPUT\'...
    %         'ground_truth3 output 28-May-2021 16-41-37\DATA N_struct100 N_k RNG_offset0 28-May-2021 16-41-37.mat'];
    %     data_path_test_iter{group_idx} = ['C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\2D-dispersion\OUTPUT\'...
    %         'gold_dataset_small output 08-Jul-2021 12-54-14\DATA N_disp100 N_wv101x51 RNG_offset0 08-Jul-2021 12-54-14.mat'];
    %     data_path_test_iter{group_idx} = ['C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\2D-dispersion\OUTPUT\'...
    %         'ground_truth4 output 08-Jul-2021 16-18-23\DATA N_disp100 N_wv501x251 RNG_offset0 08-Jul-2021 16-18-23.mat'];
    %     data_path_test_iter{group_idx} = ['C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\2D-dispersion\OUTPUT\'...
    %         'ground_truth5 4x4 output 14-Jul-2021 21-46-14\DATA N_disp100 N_wv1001x501 RNG_offset0 14-Jul-2021 21-46-14.mat'];
    %     data_path_test_iter{group_idx} = ['C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\'...
    %         '2D-dispersion\OUTPUT\light_dataset output 01-Jul-2021 17-32-47\DATA N_disp1000 N_wv51x26 01-Jul-2021 17-32-47.mat'];
    %     data_path_test_iter{group_idx} = ['C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\2D-dispersion\OUTPUT\'...
    %         'ultralight_dataset output 16-Jul-2021 10-39-56\DATA N_disp100 N_wv51x26 RNG_offset0 16-Jul-2021 10-39-56.mat'];
    %     data_path_test_iter{group_idx} = ['C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\2D-dispersion\OUTPUT\'...
    %         'superultralight_dataset output 16-Jul-2021 13-46-13\DATA N_disp100 N_wv13x7 RNG_offset0 16-Jul-2021 13-46-13.mat'];
    %     data_path_test_iter{group_idx} = ['C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\2D-dispersion\OUTPUT\'...
    %         'superultralight_dataset2 output 16-Jul-2021 16-07-08\DATA N_disp100 N_wv13x7 RNG_offset100 16-Jul-2021 16-07-08.mat'];
    %     data_path_test_iter{group_idx} = ['C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\2D-dispersion\OUTPUT\'...
    %         'superultralight_highres_dataset output 19-Jul-2021 13-38-42\DATA N_disp100 N_wv51x26 RNG_offset200 19-Jul-2021 13-38-42.mat'];
    %     data_path_test_iter{group_idx} = ['C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\2D-dispersion\OUTPUT\'...
    %         'N_pix4x4 N_ele2x2 N_wv31x16 N_disp100 N_eig10 offset1000 output 02-Oct-2021 01-16-33\DATA 02-Oct-2021 01-16-33.mat'];
    %     data_path_test_iter{group_idx} = ['C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\2D-dispersion\OUTPUT\'...
    %         'N_pix4x4 N_ele2x2 N_wv31x16 N_disp1000 N_eig10 offset0 output 02-Oct-2021 01-07-49\DATA 02-Oct-2021 01-07-49.mat'];
%     data_path_test_iter{group_idx} = ['C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\2D-dispersion\OUTPUT\'...
%         'N_pix4x4 N_ele2x2 N_wv101x51 N_disp10000 N_eig20 offset0 output 11-Jun-2021 13-24-45\DATA N_wv101x51 N_disp100 offset0 11-Jun-2021 13-24-45.mat'];
%     data_path_test_iter{group_idx} = ['C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\2D-dispersion\OUTPUT\'...
%         'N_pix4x4 N_ele2x2 N_wv151x76 N_disp100 N_eig20 offset0 output 20-Jul-2021 11-41-48\DATA N_disp100 N_wv151x76 RNG_offset0 20-Jul-2021 11-41-48.mat'];
    data_path_test_iter{group_idx} = "C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\2D-dispersion\OUTPUT\new test in distr LR output 17-Mar-2022 15-11-13\DATA N_pix16x16 N_ele1x1 N_wv31x16 N_disp100 N_eig5 offset5000 17-Mar-2022 15-11-13.mat";
%     data_path_test_iter{group_idx} = "C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\2D-dispersion\OUTPUT\new test out of distr LR output 17-Mar-2022 15-22-27\DATA N_pix16x16 N_ele1x1 N_wv31x16 N_disp100 N_eig5 offset5000 17-Mar-2022 15-22-27.mat";
end

% Define sample order =====================================================

for group_idx = 1:N_group
    %     sample_order_path_iter{group_idx} = ['C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\2D-dispersion-GPR\OUTPUT\sample_orders\'...
    %         'sample_order_data_gold4x4.mat'];
    %     sample_order_path_iter{group_idx} = ['C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\2D-dispersion-GPR\OUTPUT\sample_orders\'...
    %         'sample_order_data_gold8x8_offset100_1326.mat'];
    %     sample_order_path_iter{group_idx} = ['C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\2D-dispersion-GPR\OUTPUT\sample_orders\'...
    %         'sample_order_data_light_dataset_1326_1e-16.mat'];
    %     sample_order_path_iter{group_idx} = ['C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\2D-dispersion-GPR\OUTPUT\sample_orders\'...
    %         'sample_order_data_ultralight_dataset.mat'];
    %     sample_order_path_iter{group_idx} = ['C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\2D-dispersion-GPR\OUTPUT\sample_orders\'...
    %         'sample_order_data_superultralight_dataset1.mat'];
    %     sample_order_path_iter{group_idx} = ['C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\2D-dispersion-GPR\OUTPUT\sample_orders\'...
    %         'sample_order_data_gold8x8_offset100_5000.mat'];
    %     sample_order_path_iter{group_idx} = ['C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\2D-dispersion\OUTPUT\'...
    %         'N_pix4x4 N_ele2x2 N_wv31x16 N_disp1000 N_eig10 offset0 output 02-Oct-2021 01-07-49\sample_order_data.mat'];
%     sample_order_path_iter{group_idx} = ['C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\2D-dispersion\OUTPUT\'...
%         'N_pix4x4 N_ele2x2 N_wv101x51 N_disp10000 N_eig20 offset0 output 11-Jun-2021 13-24-45\sample_order_data.mat'];
%     sample_order_path_iter{group_idx} = ['C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\2D-dispersion\OUTPUT\'...
%         'N_pix4x4 N_ele2x2 N_wv151x76 N_disp2000 offset100 output 04-Oct-2021 11-27-08\sample_order_data.mat'];
    sample_order_path_iter{group_idx} = "C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\2D-dispersion\OUTPUT\new train output 17-Mar-2022 14-21-36\sample_order_data_50_sig0.mat";
end

% =========================================================================

for group_idx = 1:N_group
    disp(repmat('=',1,40))
    disp(['Processing group ' num2str(group_idx)])
    disp(repmat('-',1,40))
    sample_order_path = char(sample_order_path_iter{group_idx});
    data_path_train = char(data_path_train_iter{group_idx});
    data_path_test = char(data_path_test_iter{group_idx});
    
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

    disp_idxs = 1:N_disps{group_idx};
    
    if strcmp(model_names{group_idx},'GPR')
        disp('Loading training set...')
        [WAVEVECTOR_DATA_TRAIN,EIGENVALUE_DATA_TRAIN] = load_dispersion_dataset(data_path_train);
        WAVEVECTOR_DATA_TRAIN = WAVEVECTOR_DATA_TRAIN(:,:,disp_idxs);
        EIGENVALUE_DATA_TRAIN = EIGENVALUE_DATA_TRAIN(:,:,disp_idxs);
        disp('done.')
    end
    
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
        struct_idxs_iter{group_idx} = struct_idxs;
    end
    
    if strcmp(eig_idxs,'all')
        eig_idxs = 1:N_eig_test;
        eig_idxs_iter{group_idx} = eig_idxs;
    end
    
    if strcmp(model_names{group_idx},'GPR')
        assert(all(ismember(eig_idxs_iter{group_idx},sample_order_data.eig_idxs)))
    end
    
    err_L2{group_idx} = zeros(length(eig_idxs_iter{group_idx}),length(struct_idxs_iter{group_idx}),size(sample_counts_iter{group_idx},1));
    err_H1{group_idx} = zeros(length(eig_idxs_iter{group_idx}),length(struct_idxs_iter{group_idx}),size(sample_counts_iter{group_idx},1));
    
    covariance_options.isAllowGPU = false;
    covariance_options.isComputeCovarianceGradient = false;
    
    for eig_idx_idx = 1:length(eig_idxs)
        eig_idx = eig_idxs(eig_idx_idx);
        model_options = struct();
        if strcmp(model_names{group_idx},'GPR')
            original_wv_x = unique(sort(WAVEVECTOR_DATA_TRAIN(:,1,1)));
            original_wv_y = unique(sort(WAVEVECTOR_DATA_TRAIN(:,2,1)));
            
            covariance_options.eig_idxs = eig_idx;
            [Cs,C_grads,kfcns,kfcn_grads] = get_empirical_covariance(WAVEVECTOR_DATA_TRAIN,EIGENVALUE_DATA_TRAIN,covariance_options); %#ok<ASGLU>
            model_options.kfcn = kfcns{1};
            model_options.kfcn_grad = {};
        end
        
        wb_counter = 0;
        wb = waitbar(0,['Processing eig\_idx = ' num2str(eig_idx)]);
        
        %         if strcmp(model_names{group_idx},'GPR')
        %             sample_counts = [sample_counts;
        %             sample_counts_iter{group_idx} = sample_counts;
        %         end
        
        for sample_idx = 1:size(sample_counts,1)
            if ~strcmp(model_names{group_idx},'GPR')
                N_sample = sample_resolutions(sample_idx,:);
            end
            
            % to be used with model_options.predict_format = 'scattered -
            % precomputed combvec'
            if strcmp(model_names{group_idx},'GPR')
                wv_e = squeeze(WAVEVECTOR_DATA_TEST(:,:,1)); % redundant with wv = squeeze(WAVEVECTOR_DATA_TEST(:,:,struct_idx));
                [~,idxs] = sort(wv_e(:,1));
                wv_e = wv_e(idxs,:);
                wv_s = sample_order_data.sample_orders(1:sample_counts(sample_idx),:,eig_idx);
                model_options.precomputed_combvec = combvec2(wv_e',wv_s');
            end
            
            for struct_idx_idx = 1:length(struct_idxs)
                struct_idx = struct_idxs(struct_idx_idx);
                fr = squeeze(EIGENVALUE_DATA_TEST(:,eig_idx,struct_idx));
                wv = squeeze(WAVEVECTOR_DATA_TEST(:,:,struct_idx));
                model_options.model_name = model_name;
                a = 1;
                if strcmp(model_name,'GPR')
                    model_options.kfcn = kfcns{1};
                    model_options.sample_interpolation_format = 'scattered'; % specify this because these scattered points cannot be interpolated as a grid
                    model_options.train_format = 'scattered'; % only applies to GPR. must be 'scattered' even for 'gridded' sampling since the trim operation no longer leaves a gridded domain
                    model_options.predict_format = 'scattered - precomputed combvec'; % it has to be 'scattered' if train_format is 'scattered'.
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
                model_options.sigma = sigma;
                out = interp_model_2D(fr,wv,model_options);
                err_L2{group_idx}(eig_idx_idx,struct_idx_idx,sample_idx) = out.e_L2;
                err_H1{group_idx}(eig_idx_idx,struct_idx_idx,sample_idx) = out.e_H1;
            end
            wb_counter = wb_counter + 1;
            waitbar(wb_counter/size(sample_counts,1),wb)
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