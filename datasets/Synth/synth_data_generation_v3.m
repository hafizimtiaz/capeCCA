clear all;clc;close all
% This code generates synthetic samples for CCA experiments. It generates
% two sets of obervations: x1 and x2, which are the two views. View 1 is
% Gaussian and View 2 is Laplace. They will have the same mean and
% variance. Means of different classes are different, but the variance
% across all the classes is the same. See "Multi View Clustering via CCA"
% by Kamalika Chaudhury for the data generation details
% N_all = [30 50 100 200 500 800]*1000';
N_all = [5 10 20 30 50 100 200 500 800]*1000';

for nid = 1:length(N_all)
    N = N_all(nid);     % 800e3;      % number of samples
    num_class = 5;             % number of different sources
    d = 20;             % dimension of data
    K = 5;              % true dimension of data
    
    mu = linspace(0,10,num_class);  % means of different classes
    
    sigma1 = 1;              % variances in view 1 and 2
    sigma2 = 2;
    D1 = diag([sigma1 * ones(K,1); sigma1/100 * ones(d-K,1)]);
    D2 = diag([sigma2 * ones(K,1); sigma2/100 * ones(d-K,1)]);
    
    w = (1/num_class) * ones(num_class,1);      % probability of selecting source i
    x = zeros(2*d, N);                           % samples from both views together
    
    src_id = randsample(num_class, N, true, w);        % randomly selecting sources for all the samples
    
    for n = 1:N
        i = src_id(n); % select a source
        M = ones(d,1)*mu(i); % mean of source i
        
        x1 = my_rand_data_generator(M,D1,1); % sample from view 1 - Gaussian
        x2 = my_rand_data_generator(2*M-5,D1,1); % sample from view 2 - Laplace
        x(:,n) = [x1;x2];
    end
    
    % making samples mean-centered
    mu = mean(x,2);
%     x = x - mu*ones(1,N);
    
    % normalize data columns
    nrm = zeros(N,1);
    for n = 1:N
        nrm(n) = norm(x(:,n));
    end
    max_nrm = max(nrm);
    for n = 1:N
        x(:,n) = x(:,n)/max_nrm;
    end
    save(['synth_data_d',num2str(d),'_K',num2str(K),'_N',num2str(N/1000),'k.mat'],...
        'N','d','K','mu','sigma1','sigma2','w','src_id','x', 'num_class')
end