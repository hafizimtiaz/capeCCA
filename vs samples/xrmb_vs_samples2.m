clear all;clc;close all

% path = '../Data/XRMB/';
path = '';

nensemble = 100;    % number of independent runs
p_all = [10, 20, 30, 40, 50]';
epsilon = 0.1;
delta = 0.01;


CHindex_np_pool = zeros(nensemble,length(p_all));
perc_err_np_pool = zeros(nensemble,length(p_all));
mutInfo_np_pool = zeros(nensemble,length(p_all));

CHindex_dp_local = zeros(nensemble,length(p_all));
perc_err_dp_local = zeros(nensemble,length(p_all));
mutInfo_dp_local = zeros(nensemble,length(p_all));

CHindex_dp_conv1 = zeros(nensemble,length(p_all));
perc_err_dp_conv1 = zeros(nensemble,length(p_all));
mutInfo_dp_conv1 = zeros(nensemble,length(p_all));

CHindex_dp_conv2 = zeros(nensemble,length(p_all));
perc_err_dp_conv2 = zeros(nensemble,length(p_all));
mutInfo_dp_conv2 = zeros(nensemble,length(p_all));

CHindex_dp_cape = zeros(nensemble,length(p_all));
perc_err_dp_cape = zeros(nensemble,length(p_all));
mutInfo_dp_cape = zeros(nensemble,length(p_all));


score1_np_pool = zeros(nensemble,length(p_all));
score2_np_pool = zeros(nensemble,length(p_all));

score1_dp_local = zeros(nensemble,length(p_all));
score2_dp_local = zeros(nensemble,length(p_all));

score1_dp_conv1 = zeros(nensemble,length(p_all));
score2_dp_conv1 = zeros(nensemble,length(p_all));

score1_dp_conv2 = zeros(nensemble,length(p_all));
score2_dp_conv2 = zeros(nensemble,length(p_all));

score1_dp_cape = zeros(nensemble,length(p_all));
score2_dp_cape = zeros(nensemble,length(p_all));


%% run experiments
for itr = 1:nensemble
    disp(['Iteration: ',num2str(itr),' of ',num2str(nensemble)])
    
    for n_id = 1:length(p_all)
        p = p_all(n_id);
        load([path, 'XRMB_preprocessed_d50_p', num2str(p),'_new.mat'])
        K = 10;
        num_class = 2;
        
        d = 25;
        r_x = eps;             % to avoid ill-conditioned matrices
        r_y = eps;
        S = 10;
        R = 2*K;
        
        %% make data distributed
        xd = x;
        Ns = N / S;
        st_id = 1 + ([1:S] - 1) * Ns;
        en_id = [1:S] * Ns;
        site_members = zeros(S, Ns);
        x_sites = cell(S);
        for s = 1:S
            site_members(s,:) = st_id(s) : en_id(s);
            x_sites{s} = xd(:, st_id(s) : en_id(s));
        end
        src_id_local = src_id(st_id(s) : en_id(s));

    
        %% non-private pooled CCA - gold standard
        C = (1/N)*(x*x');
        C11 = C(1:d,1:d) + r_x*eye(d);
        C12 = C(1:d,d+1:2*d);
        C21 = C(d+1:2*d,1:d);
        C22 = C(d+1:2*d,d+1:2*d) + r_y*eye(d);
        
        [U1,~,~] = svd((C11\C12)*(C22\C21));
        [U2,~,~] = svd((C22\C21)*(C11\C12));
        
        U1 = U1(:,1:K);
        U2 = U2(:,1:K);
        
        x1_new = U1'*x(1:d,:);
        x2_new = U2'*x(d+1:2*d,:);
        
        [score1_np_pool(itr, n_id), score2_np_pool(itr, n_id)] = myCCAscore(x1_new,x2_new);
        [perc_err_np_pool(itr, n_id), CHindex_np_pool(itr, n_id), mutInfo_np_pool(itr, n_id)] = myClusterPerf(x1_new, src_id);
    
        
        %% using capeCCA
        tau = (1/(epsilon*N))*sqrt(2*log(1.25/delta));
        temp = normrnd(0,tau,2*d,2*d);
        temp2 = triu(temp);
        temp3 = temp2';
        temp4 = tril(temp3,-1);
        E = temp2+temp4;
        
        C_hat = C + E;
        C11 = C_hat(1:d,1:d);
        C12 = C_hat(1:d,d+1:2*d);
        C21 = C_hat(d+1:2*d,1:d);
        C22 = C_hat(d+1:2*d,d+1:2*d);
        
        [U1,~,~] = svd((C11\C12)*(C22\C21));
        [U2,~,~] = svd((C22\C21)*(C11\C12));
        
        U1 = U1(:,1:K);
        U2 = U2(:,1:K);
        
        x1_new = U1'*x(1:d,:);
        x2_new = U2'*x(d+1:2*d,:);
        
        [score1_dp_cape(itr,n_id), score2_dp_cape(itr,n_id)] = myCCAscore(x1_new,x2_new);
        [perc_err_dp_cape(itr,n_id), CHindex_dp_cape(itr,n_id), mutInfo_dp_cape(itr,n_id)] = myClusterPerf(x1_new, src_id);

        %% conv approach 1 : send Cs
        approach = 1;
        [U1_dist_dp1, U2_dist_dp1, ~, ~] = myDistCCA(x_sites, K, R, approach, epsilon, delta);
        x1_new = U1_dist_dp1'*x(1:d,:);
        x2_new = U2_dist_dp1'*x(d+1:2*d,:);
        
        [score1_dp_conv1(itr,n_id), score2_dp_conv1(itr,n_id)] = myCCAscore(x1_new,x2_new);
        [perc_err_dp_conv1(itr,n_id), CHindex_dp_conv1(itr,n_id), mutInfo_dp_conv1(itr,n_id)] = myClusterPerf(x1_new, src_id);
        
        %% conv approach 2 : send Ps
        approach = 2;
        [U1_dist_dp2, U2_dist_dp2, ~, ~] = myDistCCA(x_sites, K, R, approach, epsilon, delta);
        x1_new = U1_dist_dp2'*x(1:d,:);
        x2_new = U2_dist_dp2'*x(d+1:2*d,:);
        
        [score1_dp_conv2(itr,n_id), score2_dp_conv2(itr,n_id)] = myCCAscore(x1_new,x2_new);
        [perc_err_dp_conv2(itr,n_id), CHindex_dp_conv2(itr,n_id), mutInfo_dp_conv2(itr,n_id)] = myClusterPerf(x1_new, src_id);
        
        %% local
        tau = (S/(epsilon*N))*sqrt(2*log(1.25/delta));
        temp = normrnd(0,tau,2*d,2*d);
        temp2 = triu(temp);
        temp3 = temp2';
        temp4 = tril(temp3,-1);
        E = temp2+temp4;
        
        Clocal = (S/N)*(x_sites{S} * x_sites{S}');
        
        C_hat = Clocal + E;
        C11 = C_hat(1:d,1:d);
        C12 = C_hat(1:d,d+1:2*d);
        C21 = C_hat(d+1:2*d,1:d);
        C22 = C_hat(d+1:2*d,d+1:2*d);
        
        [U1,~,~] = svd((C11\C12)*(C22\C21));
        [U2,~,~] = svd((C22\C21)*(C11\C12));
        
        U1 = U1(:,1:K);
        U2 = U2(:,1:K);
        
        x1_new = U1' * x_sites{S}(1:d,:);
        x2_new = U2' * x_sites{S}(d+1:2*d,:);
        
        [score1_dp_local(itr,n_id), score2_dp_local(itr,n_id)] = myCCAscore(x1_new,x2_new);
        [perc_err_dp_local(itr,n_id), CHindex_dp_local(itr,n_id), mutInfo_dp_local(itr,n_id)] = myClusterPerf(x1_new, src_id_local);

    end
end
res_filename = 'results_vs_samples_xrmb_d50_K10_eps_0_1.mat';
save(res_filename,...
    'p_all','epsilon', 'delta','score1_np_pool', 'score2_np_pool', 'perc_err_np_pool', 'CHindex_np_pool','mutInfo_np_pool', ...
    'score1_dp_local', 'score2_dp_local', 'perc_err_dp_local', 'CHindex_dp_local','mutInfo_dp_local', ...
    'score1_dp_conv1', 'score2_dp_conv1', 'perc_err_dp_conv1', 'CHindex_dp_conv1','mutInfo_dp_conv1', ...
    'score1_dp_conv2', 'score2_dp_conv2', 'perc_err_dp_conv2', 'CHindex_dp_conv2','mutInfo_dp_conv2', ...
    'score1_dp_cape', 'score2_dp_cape', 'perc_err_dp_cape', 'CHindex_dp_cape','mutInfo_dp_cape');


% %% plots
% FS = 24;
% load(res_filename)
% ids = 1:5;
% 
% % for CHindex
% subplot(131);
% 
% CHindex_avg_np_pool = mean(CHindex_np_pool(:,ids));
% mutInfo_avg_np_pool = mean(mutInfo_np_pool(:,ids));
% perc_err_avg_np_pool = mean(perc_err_np_pool(:,ids));
% 
% CHindex_avg_dp_local = mean(CHindex_dp_local(:,ids));
% mutInfo_avg_dp_local = mean(mutInfo_dp_local(:,ids));
% perc_err_avg_dp_local = mean(perc_err_dp_local(:,ids));
% 
% CHindex_avg_dp_conv1 = mean(CHindex_dp_conv1(:,ids));
% mutInfo_avg_dp_conv1 = mean(mutInfo_dp_conv1(:,ids));
% perc_err_avg_dp_conv1 = mean(perc_err_dp_conv1(:,ids));
% 
% CHindex_avg_dp_cape = mean(CHindex_dp_cape(:,ids));
% mutInfo_avg_dp_cape = mean(mutInfo_dp_cape(:,ids));
% perc_err_avg_dp_cape = mean(perc_err_dp_cape(:,ids));
% 
% loglog(p_all, CHindex_avg_np_pool, 'ro-','LineWidth',3,'MarkerSize',10); hold on
% loglog(p_all, CHindex_avg_dp_local, 'r*-','LineWidth',3,'MarkerSize',10); hold on
% loglog(p_all, CHindex_avg_dp_conv1, 'bsq-','LineWidth',3,'MarkerSize',10); hold on
% loglog(p_all, CHindex_avg_dp_cape, 'kp-','LineWidth',3,'MarkerSize',10); hold on
% 
% % axis([p_all(idsDP(1))/10 p_all(end)*10 1e-4 1])
% set(gca,'FontSize',FS,'FontWeight','bold')
% xlabel('\epsilon','FontSize',FS,'FontWeight','bold');
% ylabel('CHindex','FontSize',FS,'FontWeight','bold');
% title('Synth (N = 200k, \delta = 0.01)','FontSize',FS,'FontWeight','bold')
% legend('non-priv','local','conv1','conv2','capeCCA','Location','best')
% 
% % for mutInfo
% subplot(132);
% 
% loglog(p_all, mutInfo_avg_np_pool, 'ro-','LineWidth',3,'MarkerSize',10); hold on
% loglog(p_all, mutInfo_avg_dp_local, 'r*-','LineWidth',3,'MarkerSize',10); hold on
% loglog(p_all, mutInfo_avg_dp_conv1, 'bsq-','LineWidth',3,'MarkerSize',10); hold on
% loglog(p_all, mutInfo_avg_dp_cape, 'kp-','LineWidth',3,'MarkerSize',10); hold on
% 
% % axis([p_all(idsDP(1))/10 p_all(end)*10 1e-4 1])
% set(gca,'FontSize',FS,'FontWeight','bold')
% xlabel('\epsilon','FontSize',FS,'FontWeight','bold');
% ylabel('mutInfo','FontSize',FS,'FontWeight','bold');
% title('Synth (N = 200k, \delta = 0.01)','FontSize',FS,'FontWeight','bold')
% 
% % for avgErr
% subplot(133);
% 
% loglog(p_all, perc_err_avg_np_pool, 'ro-','LineWidth',3,'MarkerSize',10); hold on
% loglog(p_all, perc_err_avg_dp_local, 'r*-','LineWidth',3,'MarkerSize',10); hold on
% loglog(p_all, perc_err_avg_dp_conv1, 'bsq-','LineWidth',3,'MarkerSize',10); hold on
% loglog(p_all, perc_err_avg_dp_cape, 'kp-','LineWidth',3,'MarkerSize',10); hold on
% 
% % axis([p_all(idsDP(1))/10 p_all(end)*10 1e-4 1])
% set(gca,'FontSize',FS,'FontWeight','bold')
% xlabel('\epsilon','FontSize',FS,'FontWeight','bold');
% ylabel('perc_err','FontSize',FS,'FontWeight','bold');
% title('Synth (N = 200k, \delta = 0.01)','FontSize',FS,'FontWeight','bold')