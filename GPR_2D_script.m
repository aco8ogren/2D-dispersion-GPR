clear; close all;

warning('off','MATLAB:MKDIR:DirectoryExists')

isSavePlots = true;

struct_idx = 1;
eig_idx = 1;

save_appendage = '';

data_path = 'C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\2D-dispersion-GPR\OUTPUT\FOR COVAR EXPER output 07-Dec-2020 15-37-06\DATA N_struct188 RNG_offset0 07-Dec-2020 15-37-06.mat';
plot_folder = replace(['C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\2D-dispersion-GPR\OUTPUT\FOR COVAR EXPER output 07-Dec-2020 15-37-06\plots/' 'struct_idx_' num2str(struct_idx) '_eig_idx_' num2str(eig_idx) '/'],'\','/');
mkdir(plot_folder)
data = load(data_path,'EIGENVALUE_DATA','WAVEVECTOR_DATA','CONSTITUTIVE_DATA');

N_sample = 3; % number of sample points sampled in the long direction of the rectangle for GPR
N_evaluate = 21;

EIGENVALUE_DATA = data.EIGENVALUE_DATA;
WAVEVECTOR_DATA = data.WAVEVECTOR_DATA;

% EIGENVALUE_DATA = reshape(EIGENVALUE_DATA,100,8,60,60);
% EIGENVALUE_DATA = permute(EIGENVALUE_DATA,[1 2 4 3]);
% EIGENVALUE_DATA = reshape(EIGENVALUE_DATA,100,8,3600);

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

plot_output(out,false,isSavePlots,save_appendage,plot_folder);

plot_output(out_sqexp,true,isSavePlots,save_appendage,plot_folder);

function plot_output(out,isUseSqexp,isSavePlots,save_appendage,plot_folder)
    disp(['e_L2 = ' num2str(out.e_L2)])
    disp(['e_H1 = ' num2str(out.e_H1)])
    
    if isUseSqexp
        save_appendage = [save_appendage 'sqexp'];
    else
        save_appendage = [save_appendage 'empirical'];
    end
    
    fig = figure2();
    hold on
    surf(out.X_e,out.Y_e,out.Z_e)
    scatter3(reshape(out.X_s,1,[]),reshape(out.Y_s,1,[]),reshape(out.Z_s,1,[]),'MarkerFaceColor','r')
    title('True')
    view(3)
    if isSavePlots
        fig = fix_pdf_border(fig);
        save_in_all_formats(fig,['true_dispersion_' save_appendage],plot_folder,true)
    end
    
    figure2()
    hold on
    surf(out.X_e,out.Y_e,out.Z_pred)
    scatter3(reshape(out.X_s,1,[]),reshape(out.Y_s,1,[]),reshape(out.Z_s,1,[]),'MarkerFaceColor','r')
    title('Predicted')
    view(3)
    if isSavePlots
        fig = fix_pdf_border(fig);
        save_in_all_formats(fig,['predicted_dispersion_' save_appendage],plot_folder,false)
    end
    
    figure2()
    hold on
    if ~isUseSqexp
        new_cov_s = covariance_function(out.wv_s',out.wv_s',out.original_domain_X,out.original_domain_Y,out.covariance);
        imagesc(new_cov_s)
    else
        imagesc(out.covariance)
    end
    set(gca,'YDir','reverse')
    colorbar('location','west')
    set(gca,'colorscale','log')
    add_top_labels(gca,out)
    if isSavePlots
        fig = fix_pdf_border(fig);
        save_in_all_formats(fig,['covariance_' save_appendage],plot_folder,false)
    end
    
end

function add_top_labels(ax,out)
    % First, store the handle to those axes.
    % Next create a second set of axes,
    % position This on top of the first and make it transparent.
    ax1=ax;
    ax2 = axes('Position', get(ax1, 'Position'),'Color', 'none');
    set(ax2, 'XAxisLocation', 'top','YAxisLocation','Right');
    % set the same Limits and Ticks on ax2 as on ax1;
    set(ax2, 'XLim', get(ax1, 'XLim'),'YLim', get(ax1, 'YLim'));
    set(ax2, 'XTick', get(ax1, 'XTick'), 'YTick', get(ax1, 'YTick'));
    % c = colorbar(ax2,'Visible','on');
    % OppTickLabels = {'a' 'b' 'c' 'd' 'e' 'f' 'g' 'h' 'i' 'j' 'k'}
    for i = 1:(size(out.wv_s,2) + 2)
        if i == 1 || i == (size(out.wv_s,2) + 2)
            OppTickLabels{i} = '';
        else
            OppTickLabels{i} = ['[' num2str(out.wv_s(1,i-1),3) ',' num2str(out.wv_s(2,i-1),3) ']'];
        end
    end
    % Set the x-tick and y-tick  labels for the second axes
    set(ax2, 'XTickLabel', OppTickLabels,'YTickLabel',OppTickLabels);
    set(ax2,'YDir','reverse');
end
