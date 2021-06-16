clear; close all;
% =========================================================================
% MUST BE RUN IN MATLAB R2020a (or maybe a later version would be okay too)
% =========================================================================

warning('off','MATLAB:MKDIR:DirectoryExists')
% warning('off','MATLAB:handle_graphics:Layout:NoPositionSetInTiledChartLayout')

addpath('../2D-dispersion-optimization-GPR')
addpath('../2D-dispersion')

% isSavePlots = false;
% isUseHomemade = true;
% isHighlightFirstStructure = true;
% isMakeBoxPlots = false;
% isMakeQuantilePlots = true;
isPlot = true;
isSaveData = true;

% Iteratable variables
eig_idxs = 'all';
struct_idxs = 'all';
model_names = {'GPR'}; % {'GPR','linear','nearest','cubic','makima','spline'};
sample_resolutions(:,1) = (2:50)'; sample_resolutions(:,2) = ceil(sample_resolutions(:,1)/2) + 1;
sigmas = [1e-2 1e-1]; % If not a scalar value or NaN, any non-GPR models will be wastefully run multiple times.


model_idxs = 1:length(model_names);

% data_path_train = ['C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\'...
%     '2D-dispersion\OUTPUT\ground_truth output 20-May-2021 17-11-07\DATA N_struct128 N_k201 RNG_offset0 20-May-2021 17-11-07.mat'];

% data_path_train = ['C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\'...
%     '2D-dispersion-GPR\OUTPUT\Homog w dataset N_k51\DATA N_struct128 N_k51 RNG_offset0 14-Mar-2021 16-46-17.mat'];

% data_path_train = ['C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\'...
%     '2D-dispersion\OUTPUT\covariance_singularity N_wv51x26 N_disp5000 output 04-Jun-2021 16-31-13\DATA N_struct5000 N_k RNG_offset0 04-Jun-2021 16-31-13.mat'];
% covariance_options.N_wv = [51 26];

data_path_train = ['C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\'...
    '2D-dispersion\OUTPUT\covariance_singularity N_wv51x26 N_disp20000 output 10-Jun-2021 14-58-54\DATA N_struct20000 N_k RNG_offset0 10-Jun-2021 14-58-54.mat'];

regexp_idx = regexp(data_path_train,'\');
data_dir = data_path_train(1:(regexp_idx(end)));
script_start_time = replace(char(datetime),':','-');
% plot_folder = replace([data_dir 'plots/covariance_analysis ' save_appendage ' N_sample_' num2str(min(N_samples)) 'to' num2str(max(N_samples)) ' N_evaluate_' num2str(N_evaluate) ' ' script_start_time '/'],'\','/');
% plot_folder = replace([data_dir 'plots/covar_analys N_samp_' num2str(min(N_samples)) 'to' num2str(max(N_samples)) ' N_eval_' num2str(N_evaluate) ' ' script_start_time '/'],'\','/');
% if isSavePlots
%     mkdir(plot_folder)
%     copyfile([mfilename('fullpath') '.m'],[plot_folder '/' mfilename '.m']);
% end

[WAVEVECTOR_DATA_TRAIN,EIGENVALUE_DATA_TRAIN] = load_dispersion_dataset(data_path_train);

[N_wv_train,N_eig_train,N_struct_train] = size(EIGENVALUE_DATA_TRAIN);

% data_path_test = ['C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\'...
%     '2D-dispersion\OUTPUT\ground_truth3 output 28-May-2021 16-41-37\DATA N_struct100 N_k RNG_offset0 28-May-2021 16-41-37.mat'];

data_path_test = ['C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\'...
    '2D-dispersion\OUTPUT\ground_truth output 20-May-2021 17-11-07\DATA N_struct128 N_k201 RNG_offset0 20-May-2021 17-11-07.mat'];

% data_path_test = ['C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\'...
%     '2D-dispersion-GPR\OUTPUT\Homog w dataset N_k51\DATA N_struct128 N_k51 RNG_offset0 14-Mar-2021 16-46-17.mat'];

[WAVEVECTOR_DATA_TEST,EIGENVALUE_DATA_TEST] = load_dispersion_dataset(data_path_test);

N_evaluate = [numel(unique(WAVEVECTOR_DATA_TEST(:,1,1))) numel(unique(WAVEVECTOR_DATA_TEST(:,2,1)))];

[N_wv_test,N_eig_test,N_struct_test] = size(EIGENVALUE_DATA_TEST);

if strcmp(struct_idxs,'all')
    struct_idxs = 1:N_struct_test;
end

if strcmp(eig_idxs,'all')
    eig_idxs = 1:N_eig_test;
end

covariance_options.isAllowGPU = false;
covariance_options.isComputeCovarianceGradient = false;

err_L2 = zeros(length(eig_idxs),length(struct_idxs),length(model_names),size(sample_resolutions,1),length(sigmas));
err_H1 = zeros(length(eig_idxs),length(struct_idxs),length(model_names),size(sample_resolutions,1),length(sigmas));

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
    for N_sample_idx = 1:size(sample_resolutions,1)
        N_sample = sample_resolutions(N_sample_idx,:);
        model_options.N_sample = N_sample;
        for struct_idx_idx = 1:length(struct_idxs)
            struct_idx = struct_idxs(struct_idx_idx);
            fr = squeeze(EIGENVALUE_DATA_TEST(:,eig_idx,struct_idx));
            wv = squeeze(WAVEVECTOR_DATA_TEST(:,:,struct_idx));
            
            for model_idx = model_idxs
                model_name = model_names{model_idx};
                model_options.model_name = model_name;
                if strcmp(model_name,'GPR')
                   model_options.kfcn = kfcns{1}; 
                end
                
                for sigma_idx = 1:length(sigmas)
                    sigma = sigmas(sigma_idx);
                    model_options.sigma = sigma;
                    out = interp_model_2D(fr,wv,model_options);
                    err_L2(eig_idx_idx,struct_idx_idx,model_idx,N_sample_idx,sigma_idx) = out.e_L2;
                    err_H1(eig_idx_idx,struct_idx_idx,model_idx,N_sample_idx,sigma_idx) = out.e_H1;
                end
            end
        end
        wb_counter = wb_counter + 1;
        waitbar(wb_counter/size(sample_resolutions,1),wb)
    end
    close(wb)
end

if isSaveData
    q = linspace(0,1,101);
    Q_L2 = quantile(err_L2,q,2);
    Q_H1 = quantile(err_H1,q,2);
    save('error_analysis_data',...
        'q','Q_L2','Q_H1','err_L2','err_H1',...
        'eig_idxs','struct_idxs','model_names','sample_resolutions','sigmas',...
        'data_path_train','data_path_test')
end

% %% New section
% 
% if isMakeBoxPlots
%     temp = categorical(model_names,model_names);
%     temp = permute(temp,[1 3 2 4]);
%     model_categories = repmat(temp,1,N_struct_test,1,1);
%     
%     fig = [figure() figure()];
%     
%     for err_idx = 1:2
%         figure(fig(err_idx))
%         tiledlayout(N_eig,size(N_samples,1))
%         for eig_idx = 1:N_eig
%             for N_sample_idx = 1:size(N_samples,1)
%                 nexttile;
%                 ax(eig_idx,N_sample_idx) = gca();
%             end
%         end
%         
%         if err_idx == 1
%             err = err_L2;
%             err_name = 'L2';
%         elseif err_idx == 2
%             err = err_H1;
%             err_name = 'H1';
%         end
%         
%         for eig_idx = eig_idxs
%             for N_sample_idx = 1:size(N_samples,1)
%                 row_idx = eig_idx;
%                 col_idx = N_sample_idx;
%                 b(row_idx,col_idx) = boxchart(ax(row_idx,col_idx),reshape(model_categories,[],1),reshape(err(eig_idx,:,:,N_sample_idx),[],1));
%                 
%                 title_appendage = ['N_{sample} = ' num2str(N_samples(N_sample_idx))];
%                 title(ax(row_idx,col_idx),[err_name ' ' 'branch ' num2str(eig_idx) ' ' title_appendage])
%                 
%                 if isHighlightFirstStructure
%                     hold(ax(row_idx,col_idx),'on');
%                     scatter(ax(row_idx,col_idx),...
%                         categorical(model_names,model_names),...
%                         err(eig_idx,1,:,N_sample_idx),'r.');
%                 end
%                 set(ax(row_idx,col_idx),'YLim',[1 10])
%                 set(ax(row_idx,col_idx),'YScale','log')
%                 set(ax(row_idx,col_idx),'YLimMode','auto')
%             end
%         end
%         set(ax,'YLim',[min(err,[],'all') max(err,[],'all')])
%         %     uniform_yscale(ax);
%         ylims = ax(1,1).YLim;
%         low_log_y = ceil(log10(ylims(1)));
%         high_log_y = floor(log10(ylims(2)));
%         newticks = 10.^(low_log_y:high_log_y);
%         yticks(ax,newticks)
%         for i = 1:numel(ax)
%             grid(ax(i),'on')
%             grid(ax(i),'minor')
%         end
%         
%         if isSavePlots
%             set(fig(err_idx),'WindowState','Maximized');
%             pause(plot_pause)
%             fig(err_idx) = fix_pdf_border(fig(err_idx));
%             save_in_all_formats(fig(err_idx),['GPR_compared_to_interpolation_' err_name '_err'],plot_folder,true);
%         end
%     end
% end
% 
% % q = linspace(0,1,3);
% q = .95;
% Q_L2 = quantile(err_L2,q,2);
% Q_H1 = quantile(err_H1,q,2);
% cm = lines(length(eig_idxs));
% if isMakeQuantilePlots
%     for err_idx = 1:2
%         if err_idx == 1
% %             err = err_L2;
%             Q = Q_L2;
%             err_name = 'L2';
%         elseif err_idx == 2
% %             err = err_H1;
%             Q = Q_H1;
%             err_name = 'H1';
%         end
%         for model_idx = 1:N_model
%             figure
%             hold on
%             for idx = 1:length(eig_idxs)
%                 eig_idx = eig_idxs(idx);
%                 p_temp = plot(prod(N_samples,2),squeeze(Q(eig_idx,:,model_idx,:)),'color',cm(idx,:));
%                 p(idx) = p_temp(1);
%                 p(idx).DisplayName = ['eig\_idx = ' num2str(eig_idx)];
%             end
%             set(gca,'yscale','log')
%             legend(p,'location','northeast')
%             title([err_name ' error quantiles ' regexprep(num2str(q),' +',',') newline model_names{model_idx}])
%             xlabel('prod(N\_sample)')
%             ylabel([err_name ' error'])
%             grid on
%             grid minor
%         end
%     end
% end
% 

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
            
            for model_idx = 1:N_model
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

