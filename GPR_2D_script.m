clear; close all;

% Perform 2D interpolation model (possibly GPR) on one eigensurface of a
% single dispersion relation contained in a dataset.

warning('off','MATLAB:MKDIR:DirectoryExists')

addpath('../2D-dispersion')

isSavePlots = false; % broken
plot_folder = '';

struct_idx = 2;
eig_idx = 10;
N_sample = [40 21]; % number of sample points sampled in the long direction of the rectangle for GPR
% N_evaluate = [1001 501]; % Number of points to evaluate error on
sigma_GPR = 1e-2;
model_name = 'GPR'; % can be GPR or any kind of interpolation method supported by interp2

covariance_options.eig_idxs = eig_idx;
covariance_options.isAllowGPU = false;
covariance_options.isComputeCovarianceGradient = false;

%% Compute covariance function

% data_path_train = ['C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\'...
%     '2D-dispersion\OUTPUT\ground_truth output 20-May-2021 17-11-07\DATA N_struct128 N_k201 RNG_offset0 20-May-2021 17-11-07.mat'];

% data_path_train = ['C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\'...
%     '2D-dispersion-GPR\OUTPUT\Homog w dataset N_k51\DATA N_struct128 N_k51 RNG_offset0 14-Mar-2021 16-46-17.mat'];

data_path_train = ['C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\'...
    '2D-dispersion\OUTPUT\covariance_singularity N_wv51x26 N_disp5000 output 04-Jun-2021 16-31-13\DATA N_struct5000 N_k RNG_offset0 04-Jun-2021 16-31-13.mat'];

% data_path_train = ['C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\'...
%     '2D-dispersion\OUTPUT\covariance_singularity N_wv101x51 N_disp10000 output 11-Jun-2021 13-24-45\DATA N_struct10000 N_k RNG_offset0 11-Jun-2021 13-24-45.mat'];

[WAVEVECTOR_DATA_TRAIN,EIGENVALUE_DATA_TRAIN] = load_dispersion_dataset(data_path_train);

% [~,N_eig_train,N_struct_train] = size(EIGENVALUE_DATA_TRAIN);

% original_wv_x = unique(sort(WAVEVECTOR_DATA_TRAIN(:,1,1)));
% original_wv_y = unique(sort(WAVEVECTOR_DATA_TRAIN(:,2,1)));

[Cs,C_grads,kfcns,kfcn_grads] = get_empirical_covariance(WAVEVECTOR_DATA_TRAIN,EIGENVALUE_DATA_TRAIN,covariance_options); %#ok<ASGLU>

% data_path_test = ['C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\'...
%     '2D-dispersion\OUTPUT\ground_truth3 output 28-May-2021 16-41-37\DATA N_struct100 N_k RNG_offset0 28-May-2021 16-41-37.mat'];

data_path_test = ['C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\'...
    '2D-dispersion\OUTPUT\ground_truth output 20-May-2021 17-11-07\DATA N_struct128 N_k201 RNG_offset0 20-May-2021 17-11-07.mat'];

% data_path_test = ['C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\'...
%     '2D-dispersion-GPR\OUTPUT\Homog w dataset N_k51\DATA N_struct128 N_k51 RNG_offset0 14-Mar-2021 16-46-17.mat'];


[WAVEVECTOR_DATA,EIGENVALUE_DATA] = load_dispersion_dataset(data_path_test);

[N_wv_test,N_eig_test,N_struct_test] = size(EIGENVALUE_DATA);

fr = EIGENVALUE_DATA(:,eig_idx,struct_idx);
wv = WAVEVECTOR_DATA(:,:,struct_idx);

model_options.model_name = model_name;
model_options.N_sample = N_sample;
if strcmp(model_options.model_name,'GPR')
    model_options.sigma = sigma_GPR;
%     model_options.kfcn = kfcns{eig_idx};
    model_options.kfcn = kfcns{1};
    model_options.kfcn_grad = {};
else
    % Do nothing if it's an interpolation method?
end

out = interp_model_2D(fr,wv,model_options);

model_options

disp('Model analysis results:')
disp(['e_L2 = ' num2str(out.e_L2)])
disp(['e_H1 = ' num2str(out.e_H1)])
% plot_output(out,isSavePlots,save_appendage,plot_folder);
% 
% disp('Squared exponential results:')
% plot_output(out_sqexp,isSavePlots,save_appendage,plot_folder);
% 
% function plot_output(out,isSavePlots,save_appendage,plot_folder)
% 
%     
%     if isUseSqexp
%         save_appendage = [save_appendage 'sqexp'];
%         sqexp_or_empir = 'sq. exp.';
%     else
%         save_appendage = [save_appendage 'empirical'];
%         sqexp_or_empir = 'empir.';
%     end
%     
%     fig = figure2();
%     ax1 = axes(fig);
%     hold('on');
%     surf(out.X_e,out.Y_e,out.Z_e)
% %     scatter3(reshape(out.X_s,1,[]),reshape(out.Y_s,1,[]),reshape(out.Z_s,1,[]),'MarkerFaceColor','r')
%     title('True')
%     view(2)
% %     colorbar;
%     fig = fix_pdf_border(fig);
%     if isSavePlots
%         save_in_all_formats(fig,['true_dispersion_' save_appendage],plot_folder,true)
%     end
%     
%     fig = figure2();
%     ax2 = axes(fig);
%     hold on
%     surf(out.X_e,out.Y_e,out.Z_pred)
% %     scatter3(reshape(out.X_s,1,[]),reshape(out.Y_s,1,[]),reshape(out.Z_s,1,[]),'MarkerFaceColor','r')
%     title(['Predicted - ' sqexp_or_empir newline 'L^2 error: ' num2str(out.e_L2) ' || H^1 error: ' num2str(out.e_H1)])
%     view(2)
% %     colorbar; ax2.CLim = ax1.CLim;
%     ax2.XLim = ax1.XLim; ax2.YLim = ax1.YLim; ax2.ZLim = ax1.ZLim;
%     fig = fix_pdf_border(fig);
%     if isSavePlots        
%         save_in_all_formats(fig,['predicted_dispersion_' save_appendage],plot_folder,false)
%     end
% end
