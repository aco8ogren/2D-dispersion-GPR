clear; close all;

warning('off','MATLAB:nearlySingularMatrix');

% data_path_train = "C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\2D-dispersion\OUTPUT\N_pix16x16 N_ele1x1 N_wv31x16 N_disp500 N_eig10 offset2100 output 02-Dec-2021 14-13-34\DATA 02-Dec-2021 14-13-34.mat";
% data_path_train = "C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\2D-dispersion\OUTPUT\output 13-Jan-2022 17-52-58\DATA N_pix16x16 N_ele1x1 N_wv101x51 N_disp20000 N_eig10 offset2000 13-Jan-2022 17-52-58.mat";
data_path_train = "C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\2D-dispersion\OUTPUT\debug train output 21-Jan-2022 16-05-00\DATA N_pix8x8 N_ele1x1 N_wv31x16 N_disp100 N_eig5 offset0 21-Jan-2022 16-05-00.mat";

N_sample = 20;
disp_idxs = 'all'; % typically set to 'all';
% N_evaluate = [31 NaN]; N_evaluate(2) = ceil(N_evaluate(1)/2);
sigma = 0;
rcond_tol = 1e-16;

isPause = false;
isPlot = false;

disp('Loading training set...')
% [WAVEVECTOR_DATA,EIGENVALUE_DATA] = load_dispersion_dataset(data_path_train);
data = load(data_path_train);
disp('done.')

if strcmp(disp_idxs,'all')
    disp_idxs = 1:size(data.EIGENVALUE_DATA,3);
end

EIGENVALUE_DATA = data.EIGENVALUE_DATA(:,:,disp_idxs);
WAVEVECTOR_DATA = data.WAVEVECTOR_DATA(:,:,disp_idxs);
[~,idxs] = sort(WAVEVECTOR_DATA(:,2,1));
EIGENVALUE_DATA = EIGENVALUE_DATA(idxs,:,:);
WAVEVECTOR_DATA = WAVEVECTOR_DATA(idxs,:,:);
N_evaluate = data.const.N_wv;

[~,N_eig,N_struct] = size(EIGENVALUE_DATA);

N_wv_train(1) = numel(unique(WAVEVECTOR_DATA(:,1,1)));
N_wv_train(2) = numel(unique(WAVEVECTOR_DATA(:,2,1)));

covariance_options.isComputeCovarianceGradient = false;
covariance_options.isAllowGPU = false;
% eig_idxs = 1:N_eig;
eig_idxs = 1:5;

[X_e,Y_e] = meshgrid(linspace(-pi,pi,N_evaluate(1)),linspace(0,pi,N_evaluate(2)));
wv_e = [reshape(X_e,[],1) reshape(Y_e,[],1)];
[~,idxs] = sort(wv_e(:,1));
wv_e = wv_e(idxs,:);

figure
tiledlayout('flow')

sample_orders = nan(N_sample,2,length(eig_idxs));
ranks = nan(length(eig_idxs),N_sample);
rconds = nan(length(eig_idxs),N_sample);

for eig_idx_idx = eig_idxs
    eig_idx = eig_idxs(eig_idx_idx);
    
    wb = waitbar(0,['Finding sampling order for eig\_idx = ' num2str(eig_idx)]);
    
    isBreak = false;
    
    covariance_options.eig_idxs = eig_idx;
    [Cs,C_grads,kfcns,kfcn_grads,Vs,V_grads,vfcns,vfcn_grads,X_grid_vec,Y_grid_vec] = get_empirical_covariance(WAVEVECTOR_DATA,EIGENVALUE_DATA,covariance_options);
    kfcn = kfcns{1};
    vfcn = vfcns{1};
    
    [V,D] = eig(reshape(Cs{1},[prod(N_wv_train) prod(N_wv_train)]));
    covariance_eigenvalues(eig_idx_idx,:) = sort(diag(D),'descend');
    
    nexttile
    title(['eig\_idx = ' num2str(eig_idx)])
    wv_s = [];
    sample_idx_labels = {};
    %     hold on
    
%     C_evaluate = kfcn(wv_e,wv_e,'gridded');
%     covar = C_evaluate;
%     
%     variances = diag(covar);
%     V = reshape(variances,N_evaluate(2),N_evaluate(1));
    
%     variances = zeros(size(wv_e,1),1);
%     for wv_idx = 1:size(wv_e,1)
%         variances(wv_idx) = kfcn(wv_e(wv_idx,:),wv_e(wv_idx,:),'gridded');
%     end
%     V = reshape(variances,N_evaluate(2),N_evaluate(1));
    
    V = vfcn(wv_e,'gridded')';
    variances = reshape(V,[],1);

    cla
    hold on
    imagesc(wv_e(:,1),wv_e(:,2),V)
    set(gca,'YDir','normal')
    colorbar
    daspect([1 1 1])
    axis tight
    
    if isPause
        pause
    end
    
    for sample_idx = 1:N_sample
        if ceil(sample_idx/10) == sample_idx/10
            waitbar(sample_idx/N_sample,wb,['Finding sampling order for eig\_idx = ' num2str(eig_idx)]);
        end
        if isBreak
            break
        end
        
        [~,idx] = max(variances);
        wv_s = [wv_s; wv_e(idx,:)];
        
        sample_idx_labels{end + 1} = num2str(sample_idx);
        
        if isPlot
            scatter(wv_s(:,1),wv_s(:,2),100,'r')
            text(wv_s(:,1),wv_s(:,2),sample_idx_labels,'HorizontalAlignment','center','VerticalAlignment','middle')
        end
        
        if isPause
            pause
        end
        
%         C_sample = kfcn(wv_s,wv_s,'scattered');
%         C_mix1 = kfcn(wv_e,wv_s,'scattered');
%         C_mix2 = kfcn(wv_s,wv_e,'scattered');
%         C_evaluate = kfcn(wv_e,wv_e,'gridded');
%         
%         covar = C_evaluate - C_mix1*((C_sample + sigma^2*eye(size(C_sample)))\C_mix2);
%         
%         variances = diag(covar);
%         V = reshape(variances,N_evaluate(2),N_evaluate(1));
        
%         variances = zeros(size(wv_e,1),1);
%         clear t1 t2 t3 t4
%         for wv_idx = 1:size(wv_e,1)
%             t1(wv_idx) = kfcn(wv_e(wv_idx,:),wv_e(wv_idx,:),'gridded');
%             t2(:,:,wv_idx) = kfcn(wv_e(wv_idx,:),wv_s,'scattered');
%             t3(:,:,wv_idx) = kfcn(wv_s,wv_s,'scattered');
%             t4(:,:,wv_idx) = kfcn(wv_s,wv_e(wv_idx,:),'scattered');
%             t_subtract(wv_idx) = t2(:,:,wv_idx)*(t3(:,:,wv_idx)\t4(:,:,wv_idx));
%             variances(wv_idx) = t1(wv_idx) - t_subtract(wv_idx);
% %             variances(wv_idx) = kfcn(wv_e(wv_idx,:),wv_e(wv_idx,:),'gridded') - kfcn(wv_e(wv_idx,:),wv_s,'scattered')*(kfcn(wv_s,wv_s,'scattered')\kfcn(wv_s,wv_e(wv_idx,:),'scattered'));
%         end
%         V = reshape(variances,N_evaluate(2),N_evaluate(1));

%         V = vfcn(wv_e,'gridded');
%         t1 = reshape(V,[],1);
%         t2 = kfcn(wv_e,wv_s,'scattered');
%         t3 = kfcn(wv_s,wv_s,'scattered');
% %         t4 = kfcn(wv_s,wv_e,'scattered');
%         t4 = t2';
%         variances = t1 - t2'*(t3\t4)';
%         V = V';

        V_temp = vfcn(wv_e,'gridded');
        t1_comp = reshape(V_temp',[],1);
        t2_comp = kfcn(wv_e,wv_s,'scattered');
        t3_comp = kfcn(wv_s,wv_s,'scattered');
        ranks(eig_idx_idx,sample_idx) = rank(t3_comp);
        rconds(eig_idx_idx,sample_idx) = rcond(t3_comp);
        if rconds(eig_idx_idx,sample_idx) < rcond_tol
            isBreak = true;
            break
        end
        t4_comp = permute(t2_comp,[2 1]);
        t2_comp = permute(t2_comp,[3 2 1]);
        t5_comp = t3_comp\t4_comp;
        t5_comp = permute(t5_comp,[1 3 2]);
        t_subtract_comp = permute(pagemtimes(t2_comp,t5_comp),[3 2 1]);
        variances_comp = t1_comp - t_subtract_comp;
        variances = variances_comp;
        V = reshape(variances,N_evaluate(2),N_evaluate(1));
        
        if isPlot
            cla
            hold on
            title(['eig\_idx = ' num2str(eig_idx)])
            imagesc(wv_e(:,1),wv_e(:,2),V)
            set(gca,'YDir','normal')
            scatter(wv_s(:,1),wv_s(:,2),100,'r')
            text(wv_s(:,1),wv_s(:,2),sample_idx_labels,'HorizontalAlignment','center','VerticalAlignment','middle')
            colorbar
            daspect([1 1 1])
            axis tight
        end
        
        if isPause
            pause
        end
        
    end
    close(wb)
    sample_orders(1:size(wv_s,1),:,eig_idx_idx) = wv_s;
    final_Vs{eig_idx_idx} = V;
end

figure
imagesc(isnan(rconds))
title('isnan(rconds)')
xlabel('sample\_idx')
ylabel('eig\_idx\_idx')

figure
plot(1:size(covariance_eigenvalues,2),covariance_eigenvalues,'.')
set(gca,'yscale','log')
ylabel('eigenvalue of covariance matrix')
xlabel('index of eigenvalue')
for eig_idx_temp = 1:size(covariance_eigenvalues,1)
    legend_labels{eig_idx_temp} = ['eig\_idx = ' num2str(eig_idx_temp)];
end
legend(legend_labels,'location','northeastoutside')

vars_to_save = {'sample_orders','eig_idxs','sigma','rconds','ranks','final_Vs','data_path_train','N_evaluate','disp_idxs'};
save('sample_order_data',vars_to_save{:})



