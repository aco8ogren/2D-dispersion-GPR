clear; close all;

warning('off','MATLAB:MKDIR:DirectoryExists')

isSavePlots = false; % broken
plot_folder = '';
isUseHomemade = true;
isMeasureRank = false;

% struct_idx = 3;
struct_idx = 1;
% eig_idx = 2;
eig_idx = 13;
N_sample = 11; % number of sample points sampled in the long direction of the rectangle for GPR
% N_evaluate = 1001; % number of points to evaluate error on
N_evaluate = 51;
sigma_GPR = 1e-16;

covariance_options.eig_idxs = eig_idx;
covariance_options.isAllowGPU = false;
covariance_options.isComputeCovarianceGradient = false;

save_appendage = '';

data_path_train = ['C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\'...
    '2D-dispersion\OUTPUT\ground_truth output 20-May-2021 17-11-07\DATA N_struct128 N_k201 RNG_offset0 20-May-2021 17-11-07.mat'];

% data_path_train = ['C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\'...
%     '2D-dispersion-GPR\OUTPUT\Homog w dataset N_k51\DATA N_struct128 N_k51 RNG_offset0 14-Mar-2021 16-46-17.mat'];

% data_train = load(data_path_train,'EIGENVALUE_DATA','WAVEVECTOR_DATA','CONSTITUTIVE_DATA');
regexp_idx = regexp(data_path_train,'\');
data_dir = data_path_train(1:(regexp_idx(end)));
script_start_time = replace(char(datetime),':','-');
% plot_folder = replace([data_dir 'plots/covariance_analysis ' save_appendage ' N_sample_' num2str(min(N_samples)) 'to' num2str(max(N_samples)) ' N_evaluate_' num2str(N_evaluate) ' ' script_start_time '/'],'\','/');
% % plot_folder = replace([data_dir 'plots/covar_analys N_samp_' num2str(min(N_samples)) 'to' num2str(max(N_samples)) ' N_eval_' num2str(N_evaluate) ' ' script_start_time '/'],'\','/');
% if isSavePlots
%     mkdir(plot_folder)
%     copyfile([mfilename('fullpath') '.m'],[plot_folder '/' mfilename '.m']);
% end

[WAVEVECTOR_DATA_TRAIN,EIGENVALUE_DATA_TRAIN] = load_dispersion_dataset(data_path_train);

[N_wv_comp,N_eig_train,N_struct] = size(EIGENVALUE_DATA_TRAIN);


original_wv_x = unique(sort(WAVEVECTOR_DATA_TRAIN(:,1,1)));
original_wv_y = unique(sort(WAVEVECTOR_DATA_TRAIN(:,2,1)));

[Cs,C_grads,kfcns,kfcn_grads] = get_empirical_covariance(WAVEVECTOR_DATA_TRAIN,EIGENVALUE_DATA_TRAIN,covariance_options); %#ok<ASGLU>

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

% % data_path = 'C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\2D-dispersion-GPR\OUTPUT\FOR COVAR EXPER output 07-Dec-2020 15-37-06\DATA N_struct188 RNG_offset0 07-Dec-2020 15-37-06.mat';
% % data_path = 'C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\2D-dispersion-GPR\OUTPUT\N_struct1024 output 10-Dec-2020 14-02-57\DATA N_struct1024 RNG_offset0 10-Dec-2020 14-02-57.mat';
% % data_path = 'C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\2D-dispersion-GPR\OUTPUT\Homog w dataset N_k51\DATA N_struct128 N_k51 RNG_offset0 14-Mar-2021 16-46-17.mat';
% regexp_idx = regexp(data_path,'\');
% data_dir = data_path(1:(regexp_idx(end)));
% plot_folder = replace([data_dir 'plots/' '_' save_appendage ' struct_idx_'...
%     num2str(struct_idx) '_eig_idx_' num2str(eig_idx) '_N_samp_'...
%     num2str(N_sample) '_N_eval_' num2str(N_evaluate) '/'],'\','/');
% if isSavePlots
%     mkdir(plot_folder)
% end
% 
% [WAVEVECTOR_DATA,EIGENVALUE_DATA] = load_dispersion_dataset(data_path);
% 
% [Cs,original_domain_X,original_domain_Y] = get_empirical_covariance(WAVEVECTOR_DATA,EIGENVALUE_DATA,covariance_options);
% 
% kfcn = @(wv_i,wv_j) covariance_function(wv_i,wv_j,original_wv_x,original_wv_y,Cs{eig_idx},covariance_options);

fr = EIGENVALUE_DATA(:,eig_idx,struct_idx);
wv = WAVEVECTOR_DATA(:,:,struct_idx);

options.isMakePlots = false;
options.isUseEmpiricalCovariance = true;
options.sigma_GPR = sigma_GPR;

kfcn = kfcns{eig_idx};

if isUseHomemade
    out = GPR2D_homemade(fr,wv,N_wv,kfcn,N_sample,N_evaluate,options);
else
    out = GPR2D(fr,wv,kfcn,N_sample,N_evaluate,options);
end

options.isUseEmpiricalCovariance = false;
if isUseHomemade
    out_sqexp = GPR2D_homemade(fr,wv,N_wv,kfcn,N_sample,N_evaluate,options);
else
    out_sqexp = GPR2D(fr,wv,kfcn,N_sample,N_evaluate,options);
end

disp('Empirical covariance results:')
plot_output(out,false,isSavePlots,save_appendage,plot_folder);

disp('Squared exponential results:')
plot_output(out_sqexp,true,isSavePlots,save_appendage,plot_folder);

function plot_output(out,isUseSqexp,isSavePlots,save_appendage,plot_folder)
    disp(['e_L2 = ' num2str(out.e_L2)])
    disp(['e_H1 = ' num2str(out.e_H1)])
    
    if isUseSqexp
        save_appendage = [save_appendage 'sqexp'];
        sqexp_or_empir = 'sq. exp.';
    else
        save_appendage = [save_appendage 'empirical'];
        sqexp_or_empir = 'empir.';
    end
    
    fig = figure2();
    ax1 = axes(fig);
    hold('on');
    surf(out.X_e,out.Y_e,out.Z_e)
%     scatter3(reshape(out.X_s,1,[]),reshape(out.Y_s,1,[]),reshape(out.Z_s,1,[]),'MarkerFaceColor','r')
    title('True')
    view(2)
%     colorbar;
    fig = fix_pdf_border(fig);
    if isSavePlots
        save_in_all_formats(fig,['true_dispersion_' save_appendage],plot_folder,true)
    end
    
    fig = figure2();
    ax2 = axes(fig);
    hold on
    surf(out.X_e,out.Y_e,out.Z_pred)
    scatter3(reshape(out.X_s,1,[]),reshape(out.Y_s,1,[]),reshape(out.Z_s,1,[]),'MarkerFaceColor','r')
    title(['Predicted - ' sqexp_or_empir newline 'L^2 error: ' num2str(out.e_L2) ' || H^1 error: ' num2str(out.e_H1)])
    view(3)
%     colorbar; ax2.CLim = ax1.CLim;
    ax2.XLim = ax1.XLim; ax2.YLim = ax1.YLim; ax2.ZLim = ax1.ZLim;
    fig = fix_pdf_border(fig);
    if isSavePlots        
        save_in_all_formats(fig,['predicted_dispersion_' save_appendage],plot_folder,false)
    end
end
