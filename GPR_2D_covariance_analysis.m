clear; close all;

warning('off','MATLAB:MKDIR:DirectoryExists')

isSavePlots = false;

N_sample = 9; % number of sample points sampled in the long direction of the rectangle for GPR
N_evaluate = 101; % number of points to evaluate error on

save_appendage = '';

% data_path = 'C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\2D-dispersion-GPR\OUTPUT\FOR COVAR EXPER output 07-Dec-2020 15-37-06\DATA N_struct188 RNG_offset0 07-Dec-2020 15-37-06.mat';
data_path = 'C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\2D-dispersion-GPR\OUTPUT\N_struct1024 output 10-Dec-2020 14-02-57\DATA N_struct1024 RNG_offset0 10-Dec-2020 14-02-57.mat';
data = load(data_path,'EIGENVALUE_DATA','WAVEVECTOR_DATA','CONSTITUTIVE_DATA');

EIGENVALUE_DATA = data.EIGENVALUE_DATA;
WAVEVECTOR_DATA = data.WAVEVECTOR_DATA;
wv = WAVEVECTOR_DATA(:,:,1);
[~,idxs] = sort(wv(:,2));
WAVEVECTOR_DATA = WAVEVECTOR_DATA(idxs,:,:);
EIGENVALUE_DATA = EIGENVALUE_DATA(idxs,:,:);

[N_wv,N_eig,N_struct] = size(EIGENVALUE_DATA);

Cs = cov4D(WAVEVECTOR_DATA,EIGENVALUE_DATA);

[N_wv,N_eig,N_struct] = size(EIGENVALUE_DATA);

pfwb = parfor_wait(N_struct,'Waitbar', true);
parfor struct_idx = 1:N_struct
    for eig_idx = 1:N_eig
        covariance = Cs{eig_idx};
        
        fr = squeeze(EIGENVALUE_DATA(:,eig_idx,struct_idx));
        wv = squeeze(WAVEVECTOR_DATA(:,:,struct_idx));
        
        options = struct();
        options.isMakePlots = false;
        options.isUseEmpiricalCovariance = true;
        
        out_emp = GPR2D(fr,wv,covariance,N_sample,N_evaluate,options);
        
        options.isUseEmpiricalCovariance = false;
        out_sqexp = GPR2D(fr,wv,covariance,N_sample,N_evaluate,options);
        
        e_L2_emp(eig_idx,struct_idx) = out_emp.e_L2;
        e_H1_emp(eig_idx,struct_idx) = out_emp.e_H1;
        e_L2_sqexp(eig_idx,struct_idx) = out_sqexp.e_L2;
        e_H1_sqexp(eig_idx,struct_idx) = out_sqexp.e_H1;
        
        %         if out_sqexp.e_L2 > 1200
        %             alphabeta = 1;
        %             isSavePlots = false; save_appendage = []; plot_folder = [];
        %             plot_output(out_sqexp,true,isSavePlots,save_appendage,plot_folder);
        %         end
    end
    pfwb.Send;
end
pfwb.Destroy;

fig = figure2();
ax(1) = axes(fig);
plot(jitter(repmat(1:N_eig,N_struct,1)),e_L2_emp','k.')
title('L2 emp')

fig = figure2();
ax(2) = axes(fig);
plot(jitter(repmat(1:N_eig,N_struct,1)),e_H1_emp','k.')
title('H1 emp')

fig = figure2();
ax(3) = axes(fig);
plot(jitter(repmat(1:N_eig,N_struct,1)),e_L2_sqexp','k.')
title('L2 sqexp')

fig = figure2();
ax(4) = axes(fig);
plot(jitter(repmat(1:N_eig,N_struct,1)),e_H1_sqexp','k.')
title('H1 sqexp')

ax(1).YLim = ax(3).YLim;
ax(2).YLim = ax(4).YLim;
