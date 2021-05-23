clear; close all;
% =========================================================================
% MUST BE RUN IN MATLAB R2020a (or maybe a later version would be okay too)
% =========================================================================

warning('off','MATLAB:MKDIR:DirectoryExists')
% warning('off','MATLAB:handle_graphics:Layout:NoPositionSetInTiledChartLayout')

isSavePlots = false;
isUseHomemade = true;
isHighlightFirstStructure = true;

% N_sample = 9; % number of sample points sampled in the long direction of the rectangle for GPR
% N_samples = 3:2:21;
N_samples = [7 13 21]; N_samp_res = length(N_samples);
N_evaluate = 151; % number of points to evaluate error on

% model_names = {'GPR','linear','nearest','cubic','makima','spline'};
model_names = {'GPR','linear'};
model_idxs = 1:length(model_names); N_model = length(model_names);

plot_pause = length(N_samples); % Give plots time to resize before trying to fix their border and save them

save_appendage = '';

% data_path = 'C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\2D-dispersion-GPR\OUTPUT\FOR COVAR EXPER output 07-Dec-2020 15-37-06\DATA N_struct188 RNG_offset0 07-Dec-2020 15-37-06.mat';
% data_path = 'C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\2D-dispersion-GPR\OUTPUT\N_struct1024 output 10-Dec-2020 14-02-57\DATA N_struct1024 RNG_offset0 10-Dec-2020 14-02-57.mat';
data_path = ['C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\'...
    '2D-dispersion-GPR\OUTPUT\Homog w dataset N_k51\DATA N_struct128 N_k51 RNG_offset0 14-Mar-2021 16-46-17.mat'];
data = load(data_path,'EIGENVALUE_DATA','WAVEVECTOR_DATA','CONSTITUTIVE_DATA');
regexp_idx = regexp(data_path,'\');
data_dir = data_path(1:(regexp_idx(end)));
script_start_time = replace(char(datetime),':','-');
plot_folder = replace([data_dir 'plots/covariance_analysis ' save_appendage ' N_sample_' num2str(min(N_samples)) 'to' num2str(max(N_samples)) ' N_evaluate_' num2str(N_evaluate) ' ' script_start_time '/'],'\','/');
% plot_folder = replace([data_dir 'plots/covar_analys N_samp_' num2str(min(N_samples)) 'to' num2str(max(N_samples)) ' N_eval_' num2str(N_evaluate) ' ' script_start_time '/'],'\','/');
if isSavePlots
    mkdir(plot_folder)
    copyfile([mfilename('fullpath') '.m'],[plot_folder '/' mfilename '.m']);
end

EIGENVALUE_DATA = data.EIGENVALUE_DATA;
WAVEVECTOR_DATA = data.WAVEVECTOR_DATA;
wv = WAVEVECTOR_DATA(:,:,1);
[~,idxs] = sort(wv(:,2));
WAVEVECTOR_DATA = WAVEVECTOR_DATA(idxs,:,:);
EIGENVALUE_DATA = EIGENVALUE_DATA(idxs,:,:);

Cs = cov4D(WAVEVECTOR_DATA,EIGENVALUE_DATA);

[N_wv,N_eig,N_struct] = size(EIGENVALUE_DATA);

err_L2 = zeros(N_eig,N_struct,length(model_names),N_samp_res);
err_H1 = zeros(N_eig,N_struct,length(model_names),N_samp_res);

f = figure2();
tl = tiledlayout((length(model_names)*2),length(N_samples));
for model_idx = 1:(length(model_names)*2)
    for N_sample_idx = 1:length(N_samples)
        nexttile;
        ax(model_idx,N_sample_idx) = gca();
    end
end

original_wv_x = unique(sort(WAVEVECTOR_DATA(:,1,1)));
original_wv_y = unique(sort(WAVEVECTOR_DATA(:,2,1)));

for eig_idx = 1:N_eig
    C_struct = struct();
    C_struct.C = Cs{eig_idx}; %#ok<PFBNS>
    C_struct.C_gpu = gpuArray(reshape(Cs{eig_idx},length(original_wv_x),length(original_wv_y),length(original_wv_x),length(original_wv_y)));
    C_struct.X_grid_vec = original_wv_x;
    C_struct.Y_grid_vec = original_wv_y;
    kfcn = @(wv_i,wv_j,query_format) covariance_function(wv_i,wv_j,query_format,C_struct);
    
    for N_sample_idx = 1:length(N_samples)
        N_sample = N_samples(N_sample_idx);
        for struct_idx = 1:N_struct
            fr = squeeze(EIGENVALUE_DATA(:,eig_idx,struct_idx));
            wv = squeeze(WAVEVECTOR_DATA(:,:,struct_idx));
            
            options = struct();
            options.isMakePlots = false;
            options.isUseEmpiricalCovariance = true;
            
            for model_idx = model_idxs
                if strcmp(model_names{model_idx},'GPR')
                    out = GPR2D_homemade(fr,wv,kfcn,N_sample,N_evaluate,options);
                else
                    out = interpolation_model_2D(fr,wv,N_sample,N_evaluate,model_names{model_idx});
                end
                err_L2(eig_idx,struct_idx,model_idx,N_sample_idx) = out.e_L2;
                err_H1(eig_idx,struct_idx,model_idx,N_sample_idx) = out.e_H1;
            end
        end
    end
    
    
    title_appendage = ['N_{sample} = ' num2str(N_sample)];
    
    for model_idx = 1:length(model_names)
        row_idx = model_idx;
        col_idx = N_sample_idx;
        %     plot(ax(model_idx,N_sample_idx),jitter(repmat(1:N_eig,N_struct,1)),err_L2','k.')
        boxchart(ax(row_idx,col_idx),err_L2(:,:,model_idx,N_sample_idx)')
        title(ax(row_idx,col_idx),['L2 ' model_names{model_idx} ' ' title_appendage])
        if isHighlightFirstStructure
            hold(ax(row_idx,col_idx),'on');
            scatter(ax(row_idx,col_idx),categorical(1:N_eig),err_L2(:,1,model_idx,N_sample_idx)','r.');
        end
        
        row_idx = model_idx + max(model_idxs);
        col_idx = N_sample_idx;
        %     plot(ax(model_idx,N_sample_idx),jitter(repmat(1:N_eig,N_struct,1)),err_L2','k.')
        boxchart(ax(row_idx,col_idx),err_H1(:,:,model_idx,N_sample_idx)')
        title(ax(row_idx,col_idx),['H1 ' model_names{model_idx} ' ' title_appendage])
        if isHighlightFirstStructure
            hold(ax(row_idx,col_idx),'on');
            scatter(ax(row_idx,col_idx),categorical(1:N_eig),err_H1(:,1,model_idx,N_sample_idx)','r.');
        end
    end
end

orig_YLims = get_orig_YLims(ax);

for f_idx = 1:7
    switch f_idx
        case 1
            max_dir = 1;
            make_y_axis_uniform(ax,max_dir);
        case 2
            max_dir = 2;
            make_y_axis_uniform(ax(1:2,:),max_dir);
            make_y_axis_uniform(ax(3:4,:),max_dir);
        case 3
            max_dir = 'all';
            make_y_axis_uniform(ax(1:2,:),max_dir);
            make_y_axis_uniform(ax(3:4,:),max_dir);
        case 4
            set(ax,'YLim',[1 10]) % workaround because of matlab's stupid log axis behavior (min(YLim) == 0 causes ax.YLim update issue)
            set(ax,'YLimMode','auto') % workaround because of matlab's stupid log axis behavior (min(YLim) == 0 causes ax.YLim update issue)
            set(ax,'YScale','log')
        case 5
            set(ax,'YLim',[1 10]) % workaround because of matlab's stupid log axis behavior (min(YLim) == 0 causes ax.YLim update issue)
            set(ax,'YLimMode','auto') % workaround because of matlab's stupid log axis behavior (min(YLim) == 0 causes ax.YLim update issue)
            set(ax,'YScale','log')
            max_dir = 1;
            make_y_axis_uniform(ax,max_dir);
        case 6
            set(ax,'YLim',[1 10]) % workaround because of matlab's stupid log axis behavior (min(YLim) == 0 causes ax.YLim update issue)
            set(ax,'YLimMode','auto') % workaround because of matlab's stupid log axis behavior (min(YLim) == 0 causes ax.YLim update issue)
            set(ax,'YScale','log')
            max_dir = 2;
            make_y_axis_uniform(ax(1:2,:),max_dir);
            make_y_axis_uniform(ax(3:4,:),max_dir);
        case 7
            set(ax,'YLim',[1 10]) % workaround because of matlab's stupid log axis behavior (min(YLim) == 0 causes ax.YLim update issue)
            set(ax,'YLimMode','auto') % workaround because of matlab's stupid log axis behavior (min(YLim) == 0 causes ax.YLim update issue)
            set(ax,'YScale','log')
            max_dir = 'all';
            make_y_axis_uniform(ax(1:2,:),max_dir);
            make_y_axis_uniform(ax(3:4,:),max_dir);
    end
    fs(f_idx) = figure2();
    copyobj(f.Children,fs(f_idx));
    set_YLims(ax,orig_YLims);
    set(ax,'YScale','linear')
end

% if isSavePlots
%     set(f,'WindowState','Maximized');
%     for i = 1:3
%         set(fs(i),'WindowState','Maximized')
%     end
%     pause(plot_pause)
%     f = fix_pdf_border(f);
%     for i = 1:3
%         fs(i) = fix_pdf_border(fs(i));
%     end
%     save_in_all_formats(f,'structure_branch_errors_independent_axes',plot_folder,true);
%     save_in_all_formats(fs(1),'structure_branch_errors_compare_N_sample',plot_folder,true);
%     save_in_all_formats(fs(2),'structure_branch_errors_compare_kernel',plot_folder,true);
%     save_in_all_formats(fs(3),'structure_branch_errors_compare_all',plot_folder,true);
% end


%% New section

N_eig = 4;

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
    
    for eig_idx = 1:N_eig
        for N_sample_idx = 1:N_samp_res
            row_idx = eig_idx;
            col_idx = N_sample_idx;
            b = boxchart(ax(row_idx,col_idx),reshape(model_categories,[],1),reshape(err(eig_idx,:,:,N_sample_idx),[],1));
            %             ax(row_idx,col_idx).XAxis.TickValues = model_names;
            
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