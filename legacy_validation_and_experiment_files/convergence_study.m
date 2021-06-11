clear; close all;

struct_idx = 1;
eig_idx = 3;

N_evaluates = [11 21 51 75 101 501 1001 2001 3001 4001 5001];

for idx = 1:length(N_evaluates)
    N_evaluate = N_evaluates(idx);
    N_sample = 3; % number of sample points sampled in the long direction of the rectangle for GPR
    % N_evaluate = 1001; % number of points to evaluate error on
    
    data_path = 'C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\2D-dispersion-GPR\OUTPUT\FOR COVAR EXPER output 07-Dec-2020 15-37-06\DATA N_struct188 RNG_offset0 07-Dec-2020 15-37-06.mat';
    data = load(data_path,'EIGENVALUE_DATA','WAVEVECTOR_DATA','CONSTITUTIVE_DATA');
    
    EIGENVALUE_DATA = data.EIGENVALUE_DATA;
    WAVEVECTOR_DATA = data.WAVEVECTOR_DATA;    
    
    % covariance = cov(squeeze(EIGENVALUE_DATA(:,eig_idx,:))');
    
    Cs = cov4D(WAVEVECTOR_DATA,EIGENVALUE_DATA);
    covariance = Cs{eig_idx};
    
    fr = squeeze(EIGENVALUE_DATA(:,eig_idx,struct_idx));
    wv = squeeze(WAVEVECTOR_DATA(:,:,struct_idx));
    
    options = struct();
    options.isMakePlots = false;
    options.isUseEmpiricalCovariance = true;
    
    out = GPR2D(fr,wv,covariance,N_sample,N_evaluate,options);
    
    options.isUseEmpiricalCovariance = false;
    out_sqexp = GPR2D(fr,wv,covariance,N_sample,N_evaluate,options);
    
    e_L2_emp(idx) = out.e_L2;
    e_H1_emp(idx) = out.e_H1;
    e_L2_sqexp(idx) = out_sqexp.e_L2;
    e_H1_sqexp(idx) = out_sqexp.e_H1;
end

figure2()
plot(N_evaluates,e_L2_emp,'-*')
title('L2 error with empirical covariance')
xlabel('N\_evaluate (~sqrt(number of evaluation points))')
ylabel('error')
set(gca,'xscale','log')
% set(gca,'yscale','log')

figure2()
plot(N_evaluates,e_H1_emp,'-*')
title('H1 error with empirical covariance')
xlabel('N\_evaluate (~sqrt(number of evaluation points))')
ylabel('error')
set(gca,'xscale','log')
% set(gca,'yscale','log')

figure2()
plot(N_evaluates,e_L2_sqexp,'-*')
title('L2 error with squared exponential covariance')
xlabel('N\_evaluate (~sqrt(number of evaluation points))')
ylabel('error')
set(gca,'xscale','log')
% set(gca,'yscale','log')

figure2()
plot(N_evaluates,e_H1_sqexp,'-*')
title('H1 error with squared exponential covariance')
xlabel('N\_evaluate (~sqrt(number of evaluation points))')
ylabel('error')
set(gca,'xscale','log')
% set(gca,'yscale','log')


