clear; close all;
gpuDevice(1);

X_cpu = rand(5e5,21^2);
x = whos('X_cpu');
disp('Data GB')
x.bytes/1e9

C_cpu = cov(X_cpu);
c = whos('C_cpu');
disp('Covariance GB')
c.bytes/1e9

X_gpu = gpuArray(X_cpu);
F_cpu = @() cov(X_cpu);
F_gpu = @() cov(X_gpu);
disp('CPU time')
timeit(F_cpu)
disp('GPU time')
gputimeit(F_gpu)
