clear; close all;
% =========================================================================
% MUST BE RUN IN MATLAB R2020a (or maybe a later version would be okay too)
% =========================================================================

warning('off','MATLAB:MKDIR:DirectoryExists')
% warning('off','MATLAB:handle_graphics:Layout:NoPositionSetInTiledChartLayout')

addpath('../2D-dispersion-optimization-GPR')

isSavePlots = false;
isUseHomemade = true;
isHighlightFirstStructure = true;
isMakeBoxPlots = false;
isMakeQuantilePlots = true;
isMakeQuantilePlots2 = true;
isSaveQuantiles = true;
eig_idxs = 1:4;
covariance_options.isAllowGPU = false;
covariance_options.isComputeCovarianceGradient = false;
sigma_GPR = 1e-4;

% N_sample = 9; % number of sample points sampled in the long direction of the rectangle for GPR
% N_samples = 3:2:21;
% N_samples = [7 9 11 13 21 31 51];
N_samples = 3:51;
N_samp_res = length(N_samples);
N_evaluate = 1001; % number of points to evaluate error on (scalar N maps to [N ceil(N/2)])

% model_names = {'GPR','linear','nearest','cubic','makima','spline'};
model_names = {'GPR'};
isUseGPR = ismember('GPR',model_names);
model_idxs = 1:length(model_names); N_model = length(model_names);

plot_pause = length(N_samples); % Give plots time to resize before trying to fix their border and save them

save_appendage = '';

% data_path_train = ['C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\'...
%     '2D-dispersion\OUTPUT\ground_truth output 20-May-2021 17-11-07\DATA N_struct128 N_k201 RNG_offset0 20-May-2021 17-11-07.mat'];

% data_path_train = ['C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\'...
%     '2D-dispersion-GPR\OUTPUT\Homog w dataset N_k51\DATA N_struct128 N_k51 RNG_offset0 14-Mar-2021 16-46-17.mat'];

data_path_train = ['C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\'...
    '2D-dispersion\OUTPUT\covariance_singularity N_wv51x26 N_disp5000 output 04-Jun-2021 16-31-13\DATA N_struct5000 N_k RNG_offset0 04-Jun-2021 16-31-13.mat'];
covariance_options.N_wv = [51 26];

regexp_idx = regexp(data_path_train,'\');
data_dir = data_path_train(1:(regexp_idx(end)));
script_start_time = replace(char(datetime),':','-');
plot_folder = replace([data_dir 'plots/covariance_analysis ' save_appendage ' N_sample_' num2str(min(N_samples)) 'to' num2str(max(N_samples)) ' N_evaluate_' num2str(N_evaluate) ' ' script_start_time '/'],'\','/');
% plot_folder = replace([data_dir 'plots/covar_analys N_samp_' num2str(min(N_samples)) 'to' num2str(max(N_samples)) ' N_eval_' num2str(N_evaluate) ' ' script_start_time '/'],'\','/');
if isSavePlots
    mkdir(plot_folder)
    copyfile([mfilename('fullpath') '.m'],[plot_folder '/' mfilename '.m']);
end

[WAVEVECTOR_DATA_TRAIN,EIGENVALUE_DATA_TRAIN] = load_dispersion_dataset(data_path_train);

[N_wv_comp,N_eig_train,N_struct] = size(EIGENVALUE_DATA_TRAIN);

% if isUseGPR
%     original_wv_x = unique(sort(WAVEVECTOR_DATA_TRAIN(:,1,1)));
%     original_wv_y = unique(sort(WAVEVECTOR_DATA_TRAIN(:,2,1)));
%     
%     [Cs,C_grads,kfcns,kfcn_grads] = get_empirical_covariance(WAVEVECTOR_DATA_TRAIN,EIGENVALUE_DATA_TRAIN,covariance_options); %#ok<ASGLU>
% end

data_path_test = ['C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\'...
    '2D-dispersion\OUTPUT\ground_truth3 output 28-May-2021 16-41-37\DATA N_struct100 N_k RNG_offset0 28-May-2021 16-41-37.mat'];

% data_path_test = ['C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\'...
%     '2D-dispersion\OUTPUT\ground_truth output 20-May-2021 17-11-07\DATA N_struct128 N_k201 RNG_offset0 20-May-2021 17-11-07.mat'];

% data_path_test = ['C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\'...
%     '2D-dispersion-GPR\OUTPUT\Homog w dataset N_k51\DATA N_struct128 N_k51 RNG_offset0 14-Mar-2021 16-46-17.mat'];

N_wv = [1001 501]; % resolution of the test data

[WAVEVECTOR_DATA,EIGENVALUE_DATA] = load_dispersion_dataset(data_path_test);

[N_wv_comp,N_eig_test,N_struct] = size(EIGENVALUE_DATA);
assert(prod(N_wv) == N_wv_comp)

err_L2 = zeros(length(eig_idxs),N_struct,length(model_names),N_samp_res);
err_H1 = zeros(length(eig_idxs),N_struct,length(model_names),N_samp_res);

for eig_idx = eig_idxs
    %     C_struct = struct();
    %     C_struct.C = Cs{eig_idx}; %#ok<PFBNS>
    %     C_struct.C_gpu = gpuArray(reshape(Cs{eig_idx},length(original_wv_x),length(original_wv_y),length(original_wv_x),length(original_wv_y)));
    %     C_struct.X_grid_vec = original_wv_x;
    %     C_struct.Y_grid_vec = original_wv_y;
    %     kfcn = @(wv_i,wv_j,query_format) covariance_function(wv_i,wv_j,query_format,C_struct,covariance_options);
    
    if isUseGPR
        original_wv_x = unique(sort(WAVEVECTOR_DATA_TRAIN(:,1,1)));
        original_wv_y = unique(sort(WAVEVECTOR_DATA_TRAIN(:,2,1)));
        
        covariance_options.eig_idxs = eig_idx;
        [Cs,C_grads,kfcns,kfcn_grads] = get_empirical_covariance(WAVEVECTOR_DATA_TRAIN,EIGENVALUE_DATA_TRAIN,covariance_options); %#ok<ASGLU>
        kfcn = kfcns{eig_idx};
    end
    
%     if isUseGPR
%         kfcn = kfcns{eig_idx};
%     end
    
    wb_counter = 0;
    wb = waitbar(0,['Processing eig\_idx = ' num2str(eig_idx)]);
    for N_sample_idx = 1:length(N_samples)
        N_sample = N_samples(N_sample_idx);
        for struct_idx = 1:N_struct
            fr = squeeze(EIGENVALUE_DATA(:,eig_idx,struct_idx));
            wv = squeeze(WAVEVECTOR_DATA(:,:,struct_idx));
            
            options = struct();
            options.isMakePlots = false;
            options.isUseEmpiricalCovariance = true;
            options.sigma_GPR = 1e-3;
            
            for model_idx = model_idxs
                if strcmp(model_names{model_idx},'GPR')
                    out = GPR2D_homemade(fr,wv,N_wv,kfcn,N_sample,N_evaluate,options);
                else
                    out = interpolation_model_2D(fr,wv,N_wv,N_sample,N_evaluate,model_names{model_idx});
                end
                err_L2(eig_idx,struct_idx,model_idx,N_sample_idx) = out.e_L2;
                err_H1(eig_idx,struct_idx,model_idx,N_sample_idx) = out.e_H1;
            end
        end
        wb_counter = wb_counter + 1;
        waitbar(wb_counter/length(N_samples),wb)
    end
    close(wb)
end

%% New section

if isMakeBoxPlots
    temp = categorical(model_names,model_names);
    temp = permute(temp,[1 3 2 4]);
    model_categories = repmat(temp,1,N_struct,1,1);
    
    fig = [figure() figure()];
    
    for err_idx = 1:2
        figure(fig(err_idx))
        tiledlayout(N_eig,N_samp_res)
        for eig_idx = 1:N_eig
            for N_sample_idx = 1:N_samp_res
                nexttile;
                ax(eig_idx,N_sample_idx) = gca();
            end
        end
        
        if err_idx == 1
            err = err_L2;
            err_name = 'L2';
        elseif err_idx == 2
            err = err_H1;
            err_name = 'H1';
        end
        
        for eig_idx = eig_idxs
            for N_sample_idx = 1:N_samp_res
                row_idx = eig_idx;
                col_idx = N_sample_idx;
                b(row_idx,col_idx) = boxchart(ax(row_idx,col_idx),reshape(model_categories,[],1),reshape(err(eig_idx,:,:,N_sample_idx),[],1));
                
                title_appendage = ['N_{sample} = ' num2str(N_samples(N_sample_idx))];
                title(ax(row_idx,col_idx),[err_name ' ' 'branch ' num2str(eig_idx) ' ' title_appendage])
                
                if isHighlightFirstStructure
                    hold(ax(row_idx,col_idx),'on');
                    scatter(ax(row_idx,col_idx),...
                        categorical(model_names,model_names),...
                        err(eig_idx,1,:,N_sample_idx),'r.');
                end
                set(ax(row_idx,col_idx),'YLim',[1 10])
                set(ax(row_idx,col_idx),'YScale','log')
                set(ax(row_idx,col_idx),'YLimMode','auto')
            end
        end
        set(ax,'YLim',[min(err,[],'all') max(err,[],'all')])
        %     uniform_yscale(ax);
        ylims = ax(1,1).YLim;
        low_log_y = ceil(log10(ylims(1)));
        high_log_y = floor(log10(ylims(2)));
        newticks = 10.^(low_log_y:high_log_y);
        yticks(ax,newticks)
        for i = 1:numel(ax)
            grid(ax(i),'on')
            grid(ax(i),'minor')
        end
        
        if isSavePlots
            set(fig(err_idx),'WindowState','Maximized');
            pause(plot_pause)
            fig(err_idx) = fix_pdf_border(fig(err_idx));
            save_in_all_formats(fig(err_idx),['GPR_compared_to_interpolation_' err_name '_err'],plot_folder,true);
        end
    end
end

% q = linspace(0,1,3);
q = .95;
Q_L2 = quantile(err_L2,q,2);
Q_H1 = quantile(err_H1,q,2);
cm = lines(length(eig_idxs));
if isMakeQuantilePlots
    for err_idx = 1:2
        if err_idx == 1
%             err = err_L2;
            Q = Q_L2;
            err_name = 'L2';
        elseif err_idx == 2
%             err = err_H1;
            Q = Q_H1;
            err_name = 'H1';
        end
        for model_idx = 1:N_model
            figure
            hold on
            for eig_idx = eig_idxs
                p_temp = plot(N_samples,squeeze(Q(eig_idx,:,model_idx,:)),'color',cm(eig_idx,:));
                p(eig_idx) = p_temp(1);
                p(eig_idx).DisplayName = ['eig\_idx = ' num2str(eig_idx)];
            end
            set(gca,'yscale','log')
            legend(p,'location','northeast')
            title([err_name ' error quantiles ' regexprep(num2str(q),' +',',') newline model_names{model_idx}])
            xlabel('N\_sample')
            ylabel([err_name ' error'])
            grid on
            grid minor
        end
    end
end

q = .95;
Q_L2 = quantile(err_L2,q,2);
Q_H1 = quantile(err_H1,q,2);
cm = lines(length(eig_idxs));
mm = {'s','d','^','x','*','o'};
if isMakeQuantilePlots2
    clear p
    for err_idx = 1:2
        if err_idx == 1
            %             err = err_L2;
            Q = Q_L2;
            err_name = 'L2';
        elseif err_idx == 2
            %             err = err_H1;
            Q = Q_H1;
            err_name = 'H1';
        end
        figure
        hold on
        for model_idx = 1:N_model
            for eig_idx = eig_idxs
                p_temp = plot(N_samples,squeeze(Q(eig_idx,:,model_idx,:)),'color',cm(eig_idx,:),'marker',mm{model_idx});
                p(eig_idx,model_idx) = p_temp(1);
                p(eig_idx,model_idx).DisplayName = [model_names{model_idx} ' eig\_idx = ' num2str(eig_idx)];
            end
            set(gca,'yscale','log')
            legend(reshape(p,[],1),'location','northeast')
            title([err_name ' error quantiles ' regexprep(num2str(q),' +',',') newline model_names{model_idx}])
            xlabel('N\_sample')
            ylabel([err_name ' error'])
            grid on
            grid minor
        end
    end
end

if isSaveQuantiles
    q = linspace(0,1,101);
    Q_L2 = quantile(err_L2,q,2);
    Q_H1 = quantile(err_H1,q,2);
    save('quantile_data','q','Q_L2','Q_H1','N_samples','model_names')
end

