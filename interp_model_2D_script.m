clear; close all;

% Perform 2D interpolation model (possibly GPR) on one eigensurface of a
% single dispersion relation contained in a dataset.

warning('off','MATLAB:MKDIR:DirectoryExists')

addpath('../2D-dispersion')

isSavePlots = false; % broken
plot_folder = '';

struct_idx = 2;
eig_idx = 1;
N_sample = []; % number of sample points sampled in the long direction of the rectangle for GPR
sample_count = 20;
% N_evaluate = [1001 501]; % Number of points to evaluate error on
sigma_GPR = 0;
model_name = 'GPR'; % can be GPR or any kind of interpolation method supported by interp2
sample_format = 'scattered';
train_format = 'scattered';
predict_format = 'scattered';

covariance_options.eig_idxs = eig_idx;
covariance_options.isAllowGPU = false;
covariance_options.isComputeCovarianceGradient = false;

%% Compute covariance function

% data_path_train = ['C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\'...
%     '2D-dispersion\OUTPUT\ground_truth output 20-May-2021 17-11-07\DATA N_struct128 N_k201 RNG_offset0 20-May-2021 17-11-07.mat'];

% data_path_train = ['C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\'...
%     '2D-dispersion-GPR\OUTPUT\Homog w dataset N_k51\DATA N_struct128 N_k51 RNG_offset0 14-Mar-2021 16-46-17.mat'];

% data_path_train = ['C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\'...
%     '2D-dispersion\OUTPUT\covariance_singularity N_wv51x26 N_disp5000 output 04-Jun-2021 16-31-13\DATA N_struct5000 N_k RNG_offset0 04-Jun-2021 16-31-13.mat'];

% data_path_train = ['C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\'...
%     '2D-dispersion\OUTPUT\covariance_singularity N_wv101x51 N_disp10000 output 11-Jun-2021 13-24-45\DATA N_struct10000 N_k RNG_offset0 11-Jun-2021 13-24-45.mat'];

data_path_train = "C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\2D-dispersion\OUTPUT\debug train output 21-Jan-2022 16-05-00\DATA N_pix8x8 N_ele1x1 N_wv31x16 N_disp100 N_eig5 offset0 21-Jan-2022 16-05-00.mat";

if strcmp(model_name,'GPR')
    [WAVEVECTOR_DATA_TRAIN,EIGENVALUE_DATA_TRAIN] = load_dispersion_dataset(data_path_train);
    
    % [~,N_eig_train,N_struct_train] = size(EIGENVALUE_DATA_TRAIN);
    
    % original_wv_x = unique(sort(WAVEVECTOR_DATA_TRAIN(:,1,1)));
    % original_wv_y = unique(sort(WAVEVECTOR_DATA_TRAIN(:,2,1)));
    
    [Cs,C_grads,kfcns,kfcn_grads] = get_empirical_covariance(WAVEVECTOR_DATA_TRAIN,EIGENVALUE_DATA_TRAIN,covariance_options); %#ok<ASGLU>
end

% data_path_test = ['C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\'...
%     '2D-dispersion\OUTPUT\ground_truth3 output 28-May-2021 16-41-37\DATA N_struct100 N_k RNG_offset0 28-May-2021 16-41-37.mat'];

% data_path_test = ['C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\'...
%     '2D-dispersion\OUTPUT\ground_truth output 20-May-2021 17-11-07\DATA N_struct128 N_k201 RNG_offset0 20-May-2021 17-11-07.mat'];

% data_path_test = ['C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\'...
%     '2D-dispersion-GPR\OUTPUT\Homog w dataset N_k51\DATA N_struct128 N_k51 RNG_offset0 14-Mar-2021 16-46-17.mat'];

data_path_test = "C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\2D-dispersion\OUTPUT\debug test output 21-Jan-2022 16-06-49\DATA N_pix16x16 N_ele1x1 N_wv51x26 N_disp20 N_eig5 offset0 21-Jan-2022 16-06-49.mat";

% sample_order_data_path = ['C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\2D-dispersion-GPR\OUTPUT\sample_orders\'...
%     'sample_order_data_light_dataset_1326_1e-16.mat'];

sample_order_data_path = "C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\2D-dispersion\OUTPUT\debug train output 21-Jan-2022 16-05-00\sample_order_data.mat";

sample_order_data = load(sample_order_data_path);

[WAVEVECTOR_DATA,EIGENVALUE_DATA] = load_dispersion_dataset(data_path_test);

[N_wv_test,N_eig_test,N_struct_test] = size(EIGENVALUE_DATA);

fr = EIGENVALUE_DATA(:,eig_idx,struct_idx);
wv = WAVEVECTOR_DATA(:,:,struct_idx);

model_options.model_name = model_name;
model_options.N_sample = N_sample;
if strcmp(model_options.model_name,'GPR')
    model_options.sigma = sigma_GPR;
    model_options.kfcn = kfcns{eig_idx};
    model_options.kfcn_grad = {};
    model_options.sample_interpolation_format = sample_format;
    model_options.train_format = train_format;
    model_options.predict_format = predict_format;
    model_options.sample_points = sample_order_data.sample_orders(1:sample_count,:,eig_idx);
else
    % Do nothing if it's an interpolation method?
end

out = interp_model_2D(fr,wv,model_options);

model_options

disp('Model analysis results:')
disp(['e_L2 = ' num2str(out.e_L2)])
disp(['e_H1 = ' num2str(out.e_H1)])

