clear all;clc;close all

load('subj1_MFCC_articulatory_data.mat');
[~,N] = size(X);
src_id1 = ones(N,1);        % true digit identifier

X = X - mean(X,2)*ones(1,N);
Y = Y - mean(Y,2)*ones(1,N);

% COVX = (1/N) * (X * X');
% [U,~,~] = svd(COVX);
% U = U(:,1:25);
% X = U'*X;
% 
% COVY = (1/N) * (Y * Y');
% [U,~,~] = svd(COVY);
% U = U(:,1:25);
% Y = U'*Y;

x1 = [X;Y];
x1 = normc(x1);
% nrm = zeros(N,1);
% for i = 1:N
%     nrm(i) = norm(x1(:,i));
% end
% max_nrm = max(nrm);
% for i = 1:N
%     x1(:,i) = x1(:,i)/max_nrm;
% end

load('subj2_MFCC_articulatory_data.mat');
[~,N] = size(X);
src_id2 = 2*ones(N,1);        % true digit identifier

X = X - mean(X,2)*ones(1,N);
Y = Y - mean(Y,2)*ones(1,N);

% COVX = (1/N) * (X * X');
% [U,~,~] = svd(COVX);
% U = U(:,1:25);
% X = U'*X;
% 
% COVY = (1/N) * (Y * Y');
% [U,~,~] = svd(COVY);
% U = U(:,1:25);
% Y = U'*Y;

x2 = [X;Y];
x2 = normc(x2);
% nrm = zeros(N,1);
% for i = 1:N
%     nrm(i) = norm(x2(:,i));
% end
% max_nrm = max(nrm);
% for i = 1:N
%     x2(:,i) = x2(:,i)/max_nrm;
% end

x = [x1 x2];
src_id = [src_id1;src_id2];
N = length(src_id);
K = 2;
dx = 7*36;
dy = 7*16;

% dx = 25;
% dy = 25;

p = 50;
x = repmat(x,1,p);
src_id = repmat(src_id,p,1);
N = N*p;

save('XRMB_preprocessed_p50_new.mat','N','K','dx','dy','src_id','x','p');
