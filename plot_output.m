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
    hold on
    surf(out.X_e,out.Y_e,out.Z_e)
    scatter3(reshape(out.X_s,1,[]),reshape(out.Y_s,1,[]),reshape(out.Z_s,1,[]),'MarkerFaceColor','r')
    title('True')
    view(3)
    if isSavePlots
        fig = fix_pdf_border(fig);
        save_in_all_formats(fig,['true_dispersion_' save_appendage],plot_folder,true)
    end
    
    fig = figure2();
    hold on
    surf(out.X_e,out.Y_e,out.Z_pred)
    scatter3(reshape(out.X_s,1,[]),reshape(out.Y_s,1,[]),reshape(out.Z_s,1,[]),'MarkerFaceColor','r')
    title(['Predicted - ' sqexp_or_empir newline 'L^2 error: ' num2str(out.e_L2) ' || H^1 error: ' num2str(out.e_H1)])
    view(3)
    if isSavePlots
        fig = fix_pdf_border(fig);
        save_in_all_formats(fig,['predicted_dispersion_' save_appendage],plot_folder,false)
    end
    
    fig = figure2();
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
%     add_top_labels(gca,out)
    if isSavePlots
        fig = fix_pdf_border(fig);
        save_in_all_formats(fig,['covariance_' save_appendage],plot_folder,false)
    end
    
end