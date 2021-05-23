clear; close all;

N_sample = 10;
N = 100;
N_interps = [10 25 40 50 70 80 100 125];

X_grid_vec = linspace(0,1,N);

setups = {'CPU','GPU (w/o bus)','GPU (w/ bus)'};
N_setup = length(setups);

times = zeros(N_sample,length(N_interps),N_setup);

figure
hold on

for setup_idx = 1:N_setup
    setup = setups{setup_idx};
    for N_interp_idx = 1:length(N_interps)
        N_interp = N_interps(N_interp_idx);
        X_grid_vec_interp = linspace(0,1,N_interp);
        for sample_idx = 1:N_sample
            if strcmp(setup,'CPU')
                A = rand(N,N,N,N);
                tic
                B = interpn(X_grid_vec,X_grid_vec,X_grid_vec,X_grid_vec,A,X_grid_vec_interp',X_grid_vec_interp,X_grid_vec_interp,X_grid_vec_interp);
                times(sample_idx,N_interp_idx,setup_idx) = toc;
            elseif strcmp(setup,'GPU (w/o bus)')
                A = gpuArray(rand(N,N,N,N));
                tic
                B = interpn(X_grid_vec,X_grid_vec,X_grid_vec,X_grid_vec,A,X_grid_vec_interp',X_grid_vec_interp,X_grid_vec_interp,X_grid_vec_interp);
                times(sample_idx,N_interp_idx,setup_idx) = toc;
            elseif strcmp(setup,'GPU (w/ bus)')
                A = rand(N,N,N,N);
                tic
                A_gpu = gpuArray(A);
                B = interpn(X_grid_vec,X_grid_vec,X_grid_vec,X_grid_vec,A_gpu,X_grid_vec_interp',X_grid_vec_interp,X_grid_vec_interp,X_grid_vec_interp);
                times(sample_idx,N_interp_idx,setup_idx) = toc;
            end
        end
    end
    
    p(setup_idx) = plot(N_interps,mean(times(:,:,setup_idx),1),'*');
    p(setup_idx).DisplayName = setup;
    title(['Time elapsed vs size of 4D interpolation grid vector' newline setup])
    xlabel(['Size of 4D interp grid vec' newline '(i.e. axis label = 10 --> 10^4 points interpolated)'])
    ylabel(['time elapsed (averaged over ' num2str(N_sample) ' samples'])
end

legend('location','northwest')