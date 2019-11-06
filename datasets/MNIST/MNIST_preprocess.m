clear all;clc;close all

load('~/OneDrive/Backup from Dropbox/Rutgers Research - Others/MNIST Database/MNIST_train.mat');
[r,c] = size(images);
N = 50000;                  % total samples I will use
n = 5000;                   % samples per class
K = 10;                     % number of classes
d = r/2;                    % dimension of each view samples
X = zeros(d,N);             % view 1 - left half of the image
Y = zeros(d,N);             % view 2 - right half of the image
src_id = zeros(N,1);        % true digit identifier

for dig = 1:K
    ids = find(labels==(dig-1));
    temp = randperm(length(ids),n);
    ids = ids(temp);
    X(:,(1+(dig-1)*n):n*dig) = images(1:r/2,ids);
    Y(:,(1+(dig-1)*n):n*dig) = images(r/2+1:r,ids);
    src_id((1+(dig-1)*n):n*dig) = dig-1;
end
X = X - mean(X,2)*ones(1,N);
Y = Y - mean(Y,2)*ones(1,N);

% COVX = (1/N) * (X * X');
% [U,~,~] = svd(COVX);
% U = U(:,1:50);
% X = U'*X;
% 
% COVY = (1/N) * (Y * Y');
% [U,~,~] = svd(COVY);
% U = U(:,1:50);
% Y = U'*Y;

x = [X;Y];
x = normc(x);
% max_nrm = sqrt(max(diag(x'*x)));
% D = diag((1/max_nrm)*ones(N,1));
% x = x*D;
% nrm = zeros(N,1);
% for i = 1:N
%     nrm(i) = norm(x(:,i));
% end
% max_nrm = max(nrm);
% for i = 1:N
%     x(:,i) = x(:,i)/max_nrm;
% end

d = 784/2;
p = 10;
x = repmat(x,1,p);
src_id = repmat(src_id,p,1);
N = N*p;
save('MNIST_preprocessed_d784_N_50k_new.mat','N','K','d','src_id','x','p','-v7.3');
