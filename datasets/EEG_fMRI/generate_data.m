clear; clc; close all

% Code for generating synthetic fMRI and EEG data.

% According to the paper:
% Correa NM, Li YO, Adal? T, Calhoun VD (2008) Canonical correlation analysis for feature-based
% fusion of biomedical imaging modalities and its application to detection of associative networks
% in schizophrenia. IEEE J Sel Topics Signal Process 2(6):998?1007

% multi-modal CCA: 1 = fMRI; 2 = EEG
D = 5; % intrinsic dimension (or the number of components)
V1 = 2500; % voxel dimension, for example 60x60 images as fMRI components
V2 = 6360; % time points for EEG signal time course

N_all = [200, 500, 700, 1000, 1300, 2000, 5000];

%% loading of fMRI components
load synth_fMRI_M1024_C20_ni250_d2500
C1 = A(:, 1:D)';

%% loading of EEG components
load EEGIFT_ica
C2 = icasig;

%% generation of the modulation profiles
for n_id = 1:length(N_all)
       
    N = N_all(n_id); % number of simulated subjects

    % generate correlated orthonormal basis
    A1 = RandOrthMat(N, 1e-10);
    A1 = A1(:, 1:D);

    A2 = RandOrthMat(N, 1e-10);
    A2 = A2(:, 1:D);

    A1n = zeros(N, D);
    A2n = zeros(N, D);
    corr_orig = [0.9, 0.7, 0.5, 0.3, 0.1];
    corr_est = zeros(D,1);
    for d = 1:D
        R = [1 corr_orig(d); corr_orig(d) 1];
        M = [A1(:, d), A2(:, d)];
        L = chol(R);
        M = M * L;
        A1n(:, d) = M(:, 1);
        A2n(:, d) = M(:, 2);
        corr_est(d) = corr(A1n(:, d), A2n(:, d));
    end

    A1_mean = mean(A1n, 1);
    A1 = A1n - ones(N, 1) * A1_mean;
    A1 = normc(A1);

    A2_mean = mean(A2n, 1);
    A2 = A2n - ones(N,1) * A2_mean;
    A2 = normc(A2);

    for d = 1:D
        corr_est(d) = corr(A1(:, d), A2(:, d));
    end

% % check
% figure; imagesc(A1' * A1)
% 
% figure; imagesc(A2' * A2)
% 
% figure; imagesc(A1' * A2)



    % generate signals
    X1 = A1 * C1;
    X2 = A2 * C2;


    % dimension reduction
    cov1 = (1/N) * (X1' * X1);
    cov2 = (1/N) * (X2' * X2);

    [u1, s1, v1] = svd(cov1);
    [u2, s2, v2] = svd(cov2);

    u1 = u1(:, 1:D);
    u2 = u2(:, 1:D);

    Y1 = X1 * u1;
    Y2 = X2 * u2;
    
    % normalize data columns
    nrm1 = zeros(N,1);
    nrm2 = zeros(N,1);
    for n = 1:N
        nrm1(n) = norm(Y1(n, :));
        nrm2(n) = norm(Y2(n, :));
    end
    
    max_nrm1 = max(nrm1);
    max_nrm2 = max(nrm2);
    
    for n = 1:N
        Y1(n, :) = Y1(n, :)/max_nrm1;
        Y2(n, :) = Y2(n, :)/max_nrm2;
    end

    save(['synth_fmri_eeg_D5_N', num2str(N)], 'X1', 'X2', 'Y1', 'Y2', 'A1', 'A2', 'C1', 'C2', 'corr_orig', 'u1', 'u2', 'corr_est')
end



