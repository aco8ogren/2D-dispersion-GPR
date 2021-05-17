clear; close all;

warning('off','MATLAB:MKDIR:DirectoryExists')

isSavePlots = false;
isUseHomemade = true;
isMeasureRank = false;

struct_idx = 3;
eig_idx = 6;
N_sample = 7; % number of sample points sampled in the long direction of the rectangle for GPR
N_evaluate = 51; % number of points to evaluate error on

save_appendage = '';

% data_path = 'C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\2D-dispersion-GPR\OUTPUT\FOR COVAR EXPER output 07-Dec-2020 15-37-06\DATA N_struct188 RNG_offset0 07-Dec-2020 15-37-06.mat';
data_path = 'C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\2D-dispersion-GPR\OUTPUT\N_struct1024 output 10-Dec-2020 14-02-57\DATA N_struct1024 RNG_offset0 10-Dec-2020 14-02-57.mat';
% data_path = 'C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\2D-dispersion-GPR\OUTPUT\Homog w dataset N_k51\DATA N_struct128 N_k51 RNG_offset0 14-Mar-2021 16-46-17.mat';
regexp_idx = regexp(data_path,'\');
data_dir = data_path(1:(regexp_idx(end)));
plot_folder = replace([data_dir 'plots/' '_' save_appendage ' struct_idx_'...
    num2str(struct_idx) '_eig_idx_' num2str(eig_idx) '_N_samp_'...
    num2str(N_sample) '_N_eval_' num2str(N_evaluate) '/'],'\','/');
if isSavePlots
    mkdir(plot_folder)
end

[WAVEVECTOR_DATA,EIGENVALUE_DATA] = load_dispersion_dataset(data_path);

[Cs,original_domain_X,original_domain_Y] = get_empirical_covariance(WAVEVECTOR_DATA,EIGENVALUE_DATA);

kfcn = @(wv_i,wv_j) covariance_function(wv_i,wv_j,original_domain_X,original_domain_Y,Cs{eig_idx});

fr = EIGENVALUE_DATA(:,eig_idx,struct_idx);
wv = WAVEVECTOR_DATA(:,:,struct_idx);

options.isMakePlots = false;
options.isUseEmpiricalCovariance = true;

if isUseHomemade
    out = GPR2D_homemade(fr,wv,kfcn,N_sample,N_evaluate,options);
else
    out = GPR2D(fr,wv,kfcn,N_sample,N_evaluate,options);
end

options.isUseEmpiricalCovariance = false;
if isUseHomemade
    out_sqexp = GPR2D_homemade(fr,wv,kfcn,N_sample,N_evaluate,options);
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
    scatter3(reshape(out.X_s,1,[]),reshape(out.Y_s,1,[]),reshape(out.Z_s,1,[]),'MarkerFaceColor','r')
    title('True')
    view(3)
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
