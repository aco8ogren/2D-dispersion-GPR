clear; close all;

quantile_of_interest = 0.95;
cmap = containers.Map;
cmap('GPR') = 'cool';
cmap('linear') = 'autumn';
model_names = {'GPR','linear'};
err_names = {'L2','H1'};

Q = containers.Map;
N_samples = containers.Map;
N_eig = containers.Map;

% load GPR quantile data
% data = load('C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\2D-dispersion-GPR\OUTPUT\ground_truth1 quantile data\quantile_data_GPR.mat');
% data = load('C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\2D-dispersion-GPR\OUTPUT\ground_truth3 quantile data\quantile_data_GPR_confusing.mat');
data = load('C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\2D-dispersion-GPR\OUTPUT\ground_truth3 quantile data\quantile_data_GPR_ttt.mat');


model_name = 'GPR';

quantile_idx = find(data.q == quantile_of_interest);
model_idx = find(strcmp(data.model_names,'GPR'));

Q([model_name ' L2']) = squeeze(data.Q_L2(:,quantile_idx,model_idx,:));
Q([model_name ' H1']) = squeeze(data.Q_H1(:,quantile_idx,model_idx,:));

N_samples('GPR') = data.N_samples;
N_eig('GPR') = size(data.Q_L2,1);
% N_eig('GPR') = 4;

% load linear interpolation quantile data
data = load('C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\2D-dispersion-GPR\OUTPUT\ground_truth1 quantile data\quantile_data_linear.mat');
% data = load('C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\2D-dispersion-GPR\OUTPUT\ground_truth3 quantile data\quantile_data_linear.mat');

model_name = 'linear';

quantile_idx = find(data.q == quantile_of_interest);
model_idx = find(strcmp(data.model_names,model_name));

Q([model_name ' L2']) = squeeze(data.Q_L2(:,quantile_idx,model_idx,:));
Q([model_name ' H1']) = squeeze(data.Q_H1(:,quantile_idx,model_idx,:));

N_samples('linear') = data.N_samples;
N_eig('linear') = size(data.Q_L2,1);
% N_eig('linear') = 4;

% N_colors = max(cell2mat(values(N_eig)));

cm = containers.Map;
for model_idx = 1:length(model_names)
    model_name = model_names{model_idx};
    cm(model_name) = get_N_colors(cmap(model_name),N_eig(model_name));
end

% mm = {'s','d'};
mm = {'none','none'};

legend_location = 'northeastoutside';

% make plots of original quantile data

for err_idx = 1:2
    figure
    hold on
    for model_idx = 1:2
        for eig_idx = 1:min(cell2mat(values(N_eig)))
            model_name = model_names{model_idx};
            Q_temp = Q([model_name ' ' err_names{err_idx}]);
            cm_temp = cm(model_name);
            p_temp = plot(N_samples(model_name),Q_temp(eig_idx,:),'color',cm_temp(eig_idx,:),'marker',mm{model_idx});
            p(eig_idx,model_idx) = p_temp(1);
            p(eig_idx,model_idx).DisplayName = [model_name ' eig\_idx = ' num2str(eig_idx)];
        end
    end
    set(gca,'yscale','log')
    axis tight
    legend(reshape(p,[],1),'location',legend_location)
    title([err_names{err_idx} ' error quantiles ' num2str(quantile_of_interest)])
    xlabel('N\_sample')
    ylabel([err_names{err_idx} ' error'])
    grid on
    grid minor
end

% make plots of smoothed quantile data
for key_cell = Q.keys
    key = char(key_cell);
    Q([key ' smooth']) = smoothdata(Q(key),2,'movmedian',10);
end

for err_idx = 1:2
    figure
    hold on
    for model_idx = 1:2
        for eig_idx = 1:min(cell2mat(values(N_eig)))
            model_name = model_names{model_idx};
            Q_temp = Q([model_name ' ' err_names{err_idx} ' smooth']);
            cm_temp = cm(model_name);
            p_temp = plot(N_samples(model_name),Q_temp(eig_idx,:),'color',cm_temp(eig_idx,:),'marker',mm{model_idx});
            p(eig_idx,model_idx) = p_temp(1);
            p(eig_idx,model_idx).DisplayName = [model_name ' eig\_idx = ' num2str(eig_idx)];
        end
    end
    set(gca,'yscale','log')
    axis tight
    legend(reshape(p,[],1),'location',legend_location)
    title([err_names{err_idx} ' error quantiles ' num2str(quantile_of_interest) newline 'smoothed'])
    xlabel('N\_sample')
    ylabel([err_names{err_idx} ' error'])
    grid on
    grid minor
end

% make plots of distances to GPR errors at specific GPR sample resolutions
N_samples_of_interest = [7 13 22];

for err_idx = 1:2
    for N_samples_of_interest_idx = 1:length(N_samples_of_interest)
        clear p
        figure
        hold on
        for eig_idx = 1:min(cell2mat(values(N_eig)))
            N_sample_GPR = N_samples_of_interest(N_samples_of_interest_idx);
            Q_temp_GPR = Q(['GPR ' err_names{err_idx} ' smooth']);
            idx = find(N_samples('GPR') == N_sample_GPR);
            Q_temp_GPR = Q_temp_GPR(:,idx);
            Q_temp_linear = Q(['linear ' err_names{err_idx} ' smooth']);
            distances = abs(Q_temp_linear(eig_idx,:) - Q_temp_GPR(eig_idx));
            N_samp_lin = N_samples('linear');
            cm_temp = cm('linear');
            p(eig_idx) = plot(N_samp_lin,distances,'color',cm_temp(eig_idx,:),'marker',mm{find(strcmp(model_names,'linear'))});
            p(eig_idx).DisplayName = ['eig\_idx = ' num2str(eig_idx)];
            [pks,locs] = findpeaks(-log10(distances),'SortStr','descend');
            
            N_samp_lin_match(eig_idx) = N_samp_lin(locs(1));
            scatter(N_samp_lin(locs(1)),10^(-pks(1)),'ko');
        end
        set(gca,'yscale','log')
        axis tight
        legend(reshape(p,[],1),'location',legend_location)
        title(['Distance of linear ' err_names{err_idx} ' error quantile ' num2str(quantile_of_interest) ' to GPR error quantile ' num2str(quantile_of_interest)...
            newline 'GPR N\_sample = ' num2str(N_sample_GPR)...
            newline 'Matching linear N\_sample = ' regexprep(num2str(N_samp_lin_match),' +',' ')])
        xlabel('N\_sample for linear model')
        ylabel(['Distance from linear ' err_names{err_idx} ' error quantile' newline 'to fixed GPR error quantile'])
        grid on
        grid minor
    end
end