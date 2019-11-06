clear all;clc;close all
% This code generates synthetic samples for CCA experiments. It generates
% two sets of obervations: x1 and x2, which are the two views. View 1 is
% Gaussian and View 2 is Laplace. They will have the same mean and
% variance. Means of different classes are different, but the variance
% across all the classes is the same. See "Multi View Clustering via CCA"
% by Kamalika Chaudhury for the data generation details
% N_all = [30 50 100 200 500 800]*1000';

N_all = [5 10 20 30 50 100 200 500 800]*1000';

d = 20;                     % dimension of data
K = 5;                      % number of sources
w = (1/K) * ones(K,1);      % probability of selecting source i

% defining the distributions
true_dim = 20;
D = {};
D{1} = diag([1 * ones(true_dim,1); 1/100 * ones(d-true_dim,1)]);
D{2} = diag([2 * ones(true_dim,1); 1/100 * ones(d-true_dim,1)]);
D{3} = diag([3 * ones(true_dim,1); 1/100 * ones(d-true_dim,1)]);
D{4} = diag([4 * ones(true_dim,1); 1/100 * ones(d-true_dim,1)]);
D{5} = diag([5 * ones(true_dim,1); 1/100 * ones(d-true_dim,1)]);

M = {};
M{1} = 1 * ones(d,1);
M{2} = 2 * ones(d,1);
M{3} = 3 * ones(d,1);
M{4} = 4 * ones(d,1);
M{5} = 5 * ones(d,1);

% generating samples
for nid = 1:length(N_all)
    N = N_all(nid);                                     % number of samples
    x = zeros(2*d, N);                                  % samples from both views together
    src_id = randsample(K, N, true, w);         % randomly selecting sources for all the samples
    
    for n = 1:N
        i = src_id(n); % select a source
        
        this_M = M{i};
        this_D = D{i};
        
        x1 = my_rand_data_generator(this_M, this_D, 1);         % sample from view 1 - Gaussian
        x2 = my_rand_data_generator(2 * this_M - 5, this_D, 1); % sample from view 2 - Laplace
        x(:,n) = [x1;x2];
    end
    
    mu = mean(x,2);
    
    % normalize data columns
    nrm = zeros(N,1);
    for n = 1:N
        nrm(n) = norm(x(:,n));
    end
    max_nrm = max(nrm);
    for n = 1:N
        x(:,n) = x(:,n)/max_nrm;
    end
    save(['synth_data_d',num2str(d),'_K',num2str(K),'_N',num2str(N/1000),'k_new.mat'],...
        'N','d','K','mu','D','M','w','src_id','x')
end