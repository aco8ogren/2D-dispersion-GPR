close all;
% =========================================================================
% MUST BE RUN IN MATLAB R2020a (or maybe a later version would be okay too)
% =========================================================================

warning('off','MATLAB:MKDIR:DirectoryExists')

isSavePlots = true;
isUseHomemade = true;

% N_sample = 9; % number of sample points sampled in the long direction of the rectangle for GPR
% N_samples = 3:2:21;
N_samples = 3:4:21;
N_evaluate = 101; % number of points to evaluate error on

plot_pause = length(N_samples); % Give plots time to resize before trying to fix their border and save them

save_appendage = 'Homemade_sig1e4_fwdslash';

data_path = 'C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\2D-dispersion-GPR\OUTPUT\FOR COVAR EXPER output 07-Dec-2020 15-37-06\DATA N_struct188 RNG_offset0 07-Dec-2020 15-37-06.mat';
% data_path = 'C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\2D-dispersion-GPR\OUTPUT\N_struct1024 output 10-Dec-2020 14-02-57\DATA N_struct1024 RNG_offset0 10-Dec-2020 14-02-57.mat';
data = load(data_path,'EIGENVALUE_DATA','WAVEVECTOR_DATA','CONSTITUTIVE_DATA');
regexp_idx = regexp(data_path,'\');
data_dir = data_path(1:(regexp_idx(end)));
script_start_time = replace(char(datetime),':','-');
plot_folder = replace([data_dir 'plots/covariance_analysis ' save_appendage ' N_sample_' num2str(min(N_samples)) 'to' num2str(max(N_samples)) ' N_evaluate_' num2str(N_evaluate) ' ' script_start_time '/'],'\','/');
% plot_folder = replace([data_dir 'plots/covar_analys N_samp_' num2str(min(N_samples)) 'to' num2str(max(N_samples)) ' N_eval_' num2str(N_evaluate) ' ' script_start_time '/'],'\','/');
mkdir(plot_folder)

copyfile([mfilename('fullpath') '.m'],[plot_folder '/' mfilename '.m']);

EIGENVALUE_DATA = data.EIGENVALUE_DATA;
WAVEVECTOR_DATA = data.WAVEVECTOR_DATA;
wv = WAVEVECTOR_DATA(:,:,1);
[~,idxs] = sort(wv(:,2));
WAVEVECTOR_DATA = WAVEVECTOR_DATA(idxs,:,:);
EIGENVALUE_DATA = EIGENVALUE_DATA(idxs,:,:);

Cs = cov4D(WAVEVECTOR_DATA,EIGENVALUE_DATA);

[N_wv,N_eig,N_struct] = size(EIGENVALUE_DATA);

f = figure2();
tl = tiledlayout(4,length(N_samples));
for plot_idx = 1:4
    for N_sample_idx = 1:length(N_samples)
        nexttile;
        ax(plot_idx,N_sample_idx) = gca();
    end
end

for N_sample_idx = 1:length(N_samples)
    N_sample = N_samples(N_sample_idx);
    pfwb = parfor_wait(N_struct,'Waitbar', true);
    parfor struct_idx = 1:N_struct
        for eig_idx = 1:N_eig
            covariance = Cs{eig_idx}; %#ok<PFBNS>
            
            fr = squeeze(EIGENVALUE_DATA(:,eig_idx,struct_idx));
            wv = squeeze(WAVEVECTOR_DATA(:,:,struct_idx));
            
            options = struct();
            options.isMakePlots = false;
            options.isUseEmpiricalCovariance = true;
            
            if isUseHomemade
                out_emp = GPR2D_homemade(fr,wv,covariance,N_sample,N_evaluate,options);
            else
                out_emp = GPR2D(fr,wv,covariance,N_sample,N_evaluate,options);
            end
            
            options.isUseEmpiricalCovariance = false;
            
            if isUseHomemade
                out_sqexp = GPR2D_homemade(fr,wv,covariance,N_sample,N_evaluate,options);
            else
                out_sqexp = GPR2D(fr,wv,covariance,N_sample,N_evaluate,options);
            end
            
            e_L2_emp(eig_idx,struct_idx) = out_emp.e_L2;
            e_H1_emp(eig_idx,struct_idx) = out_emp.e_H1;
            e_L2_sqexp(eig_idx,struct_idx) = out_sqexp.e_L2;
            e_H1_sqexp(eig_idx,struct_idx) = out_sqexp.e_H1;
        end
        pfwb.Send; %#ok<PFBNS>
    end
    pfwb.Destroy;
    
    title_appendage = ['N_{sample} = ' num2str(N_sample)];
    
    plot_idx = 2;
    % fig(plot_idx) = figure2();
    % ax(plot_idx) = axes(fig(plot_idx));
    plot(ax(plot_idx,N_sample_idx),jitter(repmat(1:N_eig,N_struct,1)),e_L2_emp','k.')
    title(ax(plot_idx,N_sample_idx),['L2 emp ' title_appendage])
    
    plot_idx = 4;
    % fig(plot_idx) = figure2();
    % ax(plot_idx) = axes(fig(plot_idx));
    plot(ax(plot_idx,N_sample_idx),jitter(repmat(1:N_eig,N_struct,1)),e_H1_emp','k.')
    title(ax(plot_idx,N_sample_idx),['H1 emp ' title_appendage])
    
    plot_idx = 1;
    % fig(plot_idx) = figure2();
    % ax(plot_idx) = axes(fig(plot_idx));
    plot(ax(plot_idx,N_sample_idx),jitter(repmat(1:N_eig,N_struct,1)),e_L2_sqexp','k.')
    title(ax(plot_idx,N_sample_idx),['L2 sqexp ' title_appendage])
    
    plot_idx = 3;
    % fig(plot_idx) = figure2();
    % ax(plot_idx) = axes(fig(plot_idx));
    plot(ax(plot_idx,N_sample_idx),jitter(repmat(1:N_eig,N_struct,1)),e_H1_sqexp','k.')
    title(ax(plot_idx,N_sample_idx),['H1 sqexp ' title_appendage])
    
    % ax(1).YLim = ax(3).YLim;
    % ax(2).YLim = ax(4).YLim;
end

orig_YLims = get_orig_YLims(ax);

for f_idx = 1:3
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
    end
    fs(f_idx) = figure2();
    copyobj(f.Children,fs(f_idx));
    set_YLims(ax,orig_YLims);
end

if isSavePlots    
    set(f,'WindowState','Maximized');
    for i = 1:3        
        set(fs(i),'WindowState','Maximized')        
    end
    pause(plot_pause)
    f = fix_pdf_border(f);
    for i = 1:3
        fs(i) = fix_pdf_border(fs(i));
    end
    save_in_all_formats(f,'structure_branch_errors_independent_axes',plot_folder,true);
    save_in_all_formats(fs(1),'structure_branch_errors_compare_N_sample',plot_folder,true);
    save_in_all_formats(fs(2),'structure_branch_errors_compare_kernel',plot_folder,true);
    save_in_all_formats(fs(3),'structure_branch_errors_compare_all',plot_folder,true);
end