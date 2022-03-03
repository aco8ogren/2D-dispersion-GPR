clear; close all;
addpath('../2D-dispersion')

% The first dataset's frequency, wavevector, eigenvector are sorted such
% that wavevector's second component varies first.
% dataset_path = "C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\2D-dispersion\OUTPUT\debug output 24-Feb-2022 17-05-46\data.mat";

% The second dataset's frequency, wavevector, eigenvector are sorted such
% that wavevector's first component varies first
dataset_path = "C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\2D-dispersion\OUTPUT\debug2 output 24-Feb-2022 17-29-32\data.mat";

data = load(dataset_path);
dispersion_dataset = data.dispersion_dataset;
disp('wavevector:')
disp(dispersion_dataset.wavevector(1:10,:))

ecfc = EmpiricalCovarianceFunctionConstructor(dataset_path);

ecfc = ecfc.run;

% We can see from the plots that the ordering in the second dataset is the
% proper ordering, given the permute-less covariance array computation in
% EmpiricalCovarianceFunctionConstructor.
% This ordering produces plots that show that the 4D covariance array is
% continuous when slicing it along (all?) of its axes.

figure
plot(squeeze(ecfc.covariance_array{1}(:,1,1,1)))

figure
plot(squeeze(ecfc.covariance_array{1}(1,:,1,1)))

figure
plot(squeeze(ecfc.covariance_array{1}(1,1,:,1)))

figure
plot(squeeze(ecfc.covariance_array{1}(1,1,1,:)))

figure
imagesc(squeeze(ecfc.covariance_array{1}(1,1,:,:)))

figure
imagesc(squeeze(ecfc.covariance_array{1}(:,1,1,:)))

figure
imagesc(squeeze(ecfc.covariance_array{1}(:,:,1,1)))

figure
imagesc(squeeze(ecfc.covariance_array{1}(1,:,:,1)))

figure
imagesc(squeeze(ecfc.covariance_array{1}(1,:,1,:)))

figure
imagesc(squeeze(ecfc.covariance_array{1}(:,1,:,1)))
