clear; clc; close all

FS = 24;

%% vs delta
figure

% % synth
% load(['results_vs_del_synth_d20_K5_N50k_eps_1.mat'])
% 
% % for s1
% subplot(151);
% ids = 1:8;
% 
% score1_np_pool = repmat(score1_np_pool, 1, length(delta_all));
% score1_avg_np_pool = mean(score1_np_pool(:,ids));
% 
% score2_np_pool = repmat(score2_np_pool, 1, length(delta_all));
% score2_avg_np_pool = mean(score2_np_pool(:,ids));
% 
% score1_avg_dp_local = mean(score1_dp_local(:,ids));
% score2_avg_dp_local = mean(score2_dp_local(:,ids));
% 
% score1_avg_dp_conv1 = mean(score1_dp_conv1(:,ids));
% score2_avg_dp_conv1 = mean(score2_dp_conv1(:,ids));
% 
% score1_avg_dp_cape = mean(score1_dp_cape(:,ids));
% score2_avg_dp_cape = mean(score2_dp_cape(:,ids));
% 
% semilogx(delta_all, score1_avg_np_pool, 'ro--','LineWidth',3,'MarkerSize',10); hold on
% semilogx(delta_all, score1_avg_dp_local, 'r*--','LineWidth',3,'MarkerSize',10); hold on
% semilogx(delta_all, score1_avg_dp_conv1, 'bsq:','LineWidth',3,'MarkerSize',10); hold on
% semilogx(delta_all, score1_avg_dp_cape, 'kp-.','LineWidth',3,'MarkerSize',10); hold on
% 
% axis([delta_all(1)/10 delta_all(end)*10 1e-4 8e-3])
% xticks([1e-5, 1e-4, 1e-3, 1e-2, 0.05])
% xticklabels({'1e^{-5}','1e^{-4}','1e^{-3}', '1e^{-2}', '5e^{-2}'})
% set(gca,'FontSize',FS,'FontWeight','bold')
% xlabel('(a) Privacy param (\delta)','FontSize',FS,'FontWeight','bold');
% ylabel('s_1','FontSize',FS,'FontWeight','bold');
% title('Synth (N = 50k, \epsilon = 1.0)','FontSize',FS,'FontWeight','bold')
% legend('non-priv','local','conv','capeCCA','Location','best')
% 
% % for s2
% subplot(152);
% semilogx(delta_all, score2_avg_np_pool, 'ro--','LineWidth',3,'MarkerSize',10); hold on
% semilogx(delta_all, score2_avg_dp_local, 'r*--','LineWidth',3,'MarkerSize',10); hold on
% semilogx(delta_all, score2_avg_dp_conv1, 'bsq:','LineWidth',3,'MarkerSize',10); hold on
% semilogx(delta_all, score2_avg_dp_cape, 'kp-.','LineWidth',3,'MarkerSize',10); hold on
% 
% axis([delta_all(1)/10 delta_all(end)*10 5e-3 1e-1])
% xticks([1e-5, 1e-4, 1e-3, 1e-2, 0.05])
% xticklabels({'1e^{-5}','1e^{-4}','1e^{-3}', '1e^{-2}', '5e^{-2}'})
% set(gca,'FontSize',FS,'FontWeight','bold')
% xlabel('(b) Privacy param (\delta)','FontSize',FS,'FontWeight','bold');
% ylabel('s_2','FontSize',FS,'FontWeight','bold');
% title('Synth (N = 50k, \epsilon = 1.0)','FontSize',FS,'FontWeight','bold')

% MNIST
load(['results_vs_del_mnist_d100_K50_N50k_eps_0_05.mat'])
ids = 1:length(delta_all)

% for CHindex
subplot(131)
CHindex_np_pool = repmat(CHindex_np_pool, 1, length(delta_all));
CHindex_avg_np_pool = mean(CHindex_np_pool(:,ids));

CHindex_avg_dp_local = mean(CHindex_dp_local(:,ids));

CHindex_avg_dp_conv1 = mean(CHindex_dp_conv1(:,ids));

CHindex_avg_dp_cape = mean(CHindex_dp_cape(:,ids));

semilogx(delta_all, CHindex_avg_np_pool, 'ro--','LineWidth',3,'MarkerSize',10); hold on
semilogx(delta_all, CHindex_avg_dp_local, 'r*--','LineWidth',3,'MarkerSize',10); hold on
semilogx(delta_all, CHindex_avg_dp_conv1, 'bsq:','LineWidth',3,'MarkerSize',10); hold on
semilogx(delta_all, CHindex_avg_dp_cape, 'kp-.','LineWidth',3,'MarkerSize',10); hold on

axis([delta_all(1)/10 delta_all(end)*10 1e-4 3e4])
xticks([1e-5, 1e-4, 1e-3, 1e-2, 0.05])
xticklabels({'1e^{-5}','1e^{-4}','1e^{-3}', '1e^{-2}', '5e^{-2}'})
set(gca,'FontSize',FS,'FontWeight','bold')
xlabel('(a) Privacy param (\delta)','FontSize',FS,'FontWeight','bold');
ylabel('CHindex','FontSize',FS,'FontWeight','bold');
title('MNIST (N = 50k, \epsilon = 0.05)','FontSize',FS,'FontWeight','bold')
legend('non-priv','local','conv','capeCCA','Location','best')

% XRMB
load(['results_vs_del_xrmb_d50_K20_p30_eps_0_2.mat'])
ids = 1:length(delta_all)

% for CHindex
subplot(132)
CHindex_np_pool = repmat(CHindex_np_pool, 1, length(delta_all));
CHindex_avg_np_pool = mean(CHindex_np_pool(:,ids));

CHindex_avg_dp_local = mean(CHindex_dp_local(:,ids));

CHindex_avg_dp_conv1 = mean(CHindex_dp_conv1(:,ids));

CHindex_avg_dp_cape = mean(CHindex_dp_cape(:,ids));

semilogx(delta_all, CHindex_avg_np_pool, 'ro--','LineWidth',3,'MarkerSize',10); hold on
semilogx(delta_all, CHindex_avg_dp_local, 'r*--','LineWidth',3,'MarkerSize',10); hold on
semilogx(delta_all, CHindex_avg_dp_conv1, 'bsq:','LineWidth',3,'MarkerSize',10); hold on
semilogx(delta_all, CHindex_avg_dp_cape, 'kp-.','LineWidth',3,'MarkerSize',10); hold on

axis([delta_all(1)/10 delta_all(end)*10 1e-4 5e4])
xticks([1e-5, 1e-4, 1e-3, 1e-2, 0.05])
xticklabels({'1e^{-5}','1e^{-4}','1e^{-3}', '1e^{-2}', '5e^{-2}'})
set(gca,'FontSize',FS,'FontWeight','bold')
xlabel('(b) Privacy param (\delta)','FontSize',FS,'FontWeight','bold');
ylabel('CHindex','FontSize',FS,'FontWeight','bold');
title('XRMB (p = 30, \epsilon = 0.2)','FontSize',FS,'FontWeight','bold')

% for err_corr of fMRI+EEG
load results_vs_del_synth_D5_eps_0_5_N2000

% for err_corr
subplot(133);
ids = 1:length(delta_all);

err_np_pool = repmat(err_np_pool, 1, length(delta_all));
err_avg_np_pool = mean(err_np_pool(:,ids));

err_avg_dp_local = mean(err_dp_local(:,ids));

err_avg_dp_conv1 = mean(err_dp_conv1(:,ids));

err_avg_dp_cape = mean(err_dp_cape(:,ids));

semilogx(delta_all, err_avg_np_pool, 'ro--','LineWidth',3,'MarkerSize',10); hold on
semilogx(delta_all, err_avg_dp_local, 'r*--','LineWidth',3,'MarkerSize',10); hold on
semilogx(delta_all, err_avg_dp_conv1, 'bsq:','LineWidth',3,'MarkerSize',10); hold on
semilogx(delta_all, err_avg_dp_cape, 'kp-.','LineWidth',3,'MarkerSize',10); hold on

xticks([1e-5, 1e-4, 1e-3, 1e-2, 0.05])
xticklabels({'1e^{-5}','1e^{-4}','1e^{-3}', '1e^{-2}', '5e^{-2}'})
axis([delta_all(1)/10 delta_all(end)*10 1e-4 0.2])
set(gca,'FontSize',FS,'FontWeight','bold')
xlabel('(c) Privacy param (\delta)','FontSize',FS,'FontWeight','bold');
ylabel('err_{corr}','FontSize',FS,'FontWeight','bold');
title('fMRI+EEG (N = 2000, \epsilon = 0.5)','FontSize',FS,'FontWeight','bold')

%% for synth vs all
figure
% vs epsilon
load('results_vs_eps_synth_d20_K5_N10k.mat')
ids = 1:7;

% for s1
subplot(141)
score1_np_pool = repmat(score1_np_pool, 1, length(epsilon_all));
score1_avg_np_pool = mean(score1_np_pool(:,ids));

score2_np_pool = repmat(score2_np_pool, 1, length(epsilon_all));
score2_avg_np_pool = mean(score2_np_pool(:,ids));

score1_avg_dp_local = mean(score1_dp_local(:,ids));
score2_avg_dp_local = mean(score2_dp_local(:,ids));

score1_avg_dp_conv1 = mean(score1_dp_conv1(:,ids));
score2_avg_dp_conv1 = mean(score2_dp_conv1(:,ids));

score1_avg_dp_cape = mean(score1_dp_cape(:,ids));
score2_avg_dp_cape = mean(score2_dp_cape(:,ids));

semilogx(epsilon_all, score1_avg_np_pool, 'ro--','LineWidth',3,'MarkerSize',10); hold on
semilogx(epsilon_all, score1_avg_dp_local, 'r*--','LineWidth',3,'MarkerSize',10); hold on
semilogx(epsilon_all, score1_avg_dp_conv1, 'bsq:','LineWidth',3,'MarkerSize',10); hold on
semilogx(epsilon_all, score1_avg_dp_cape, 'kp-.','LineWidth',3,'MarkerSize',10); hold on

xticks([1e-3, 1e-2, 1e-1, 1, 10])
xticklabels({'1e^{-3}','1e^{-2}','1e^{-1}', '1', '10'})
axis([epsilon_all(1)/10 epsilon_all(end)*10 0.002 0.018])
set(gca,'FontSize',FS,'FontWeight','bold')
xlabel('(a) Privacy param (\epsilon)','FontSize',FS,'FontWeight','bold');
ylabel('s_1','FontSize',FS,'FontWeight','bold');
title('Synth (N = 10k)','FontSize',FS,'FontWeight','bold')
legend('non-priv','local','conv','capeCCA','Location','best')

% for s2
subplot(142);
semilogx(epsilon_all, score2_avg_np_pool, 'ro--','LineWidth',3,'MarkerSize',10); hold on
semilogx(epsilon_all, score2_avg_dp_local, 'r*--','LineWidth',3,'MarkerSize',10); hold on
semilogx(epsilon_all, score2_avg_dp_conv1, 'bsq:','LineWidth',3,'MarkerSize',10); hold on
semilogx(epsilon_all, score2_avg_dp_cape, 'kp-.','LineWidth',3,'MarkerSize',10); hold on

xticks([1e-3, 1e-2, 1e-1, 1, 10])
xticklabels({'1e^{-3}','1e^{-2}','1e^{-1}', '1', '10'})
axis([epsilon_all(1)/10 epsilon_all(end)*10 0.01 0.08])
set(gca,'FontSize',FS,'FontWeight','bold')
xlabel('(b) Privacy param (\epsilon)','FontSize',FS,'FontWeight','bold');
ylabel('s_2','FontSize',FS,'FontWeight','bold');
title('Synth (N = 10k)','FontSize',FS,'FontWeight','bold')


% vs samples
load('results_vs_samples_synth_d20_K5_eps_0_5.mat')
ids = 1:8;

% for s1
subplot(143)
score1_np_pool = repmat(score1_np_pool, 1, length(N_all));
score1_avg_np_pool = mean(score1_np_pool(:,ids));

score2_np_pool = repmat(score2_np_pool, 1, length(N_all));
score2_avg_np_pool = mean(score2_np_pool(:,ids));

score1_avg_dp_local = mean(score1_dp_local(:,ids));
score2_avg_dp_local = mean(score2_dp_local(:,ids));

score1_avg_dp_conv1 = mean(score1_dp_conv1(:,ids));
score2_avg_dp_conv1 = mean(score2_dp_conv1(:,ids));

score1_avg_dp_cape = mean(score1_dp_cape(:,ids));
score2_avg_dp_cape = mean(score2_dp_cape(:,ids));

semilogx(N_all, score1_avg_np_pool, 'ro--','LineWidth',3,'MarkerSize',10); hold on
semilogx(N_all, score1_avg_dp_local, 'r*--','LineWidth',3,'MarkerSize',10); hold on
semilogx(N_all, score1_avg_dp_conv1, 'bsq:','LineWidth',3,'MarkerSize',10); hold on
semilogx(N_all, score1_avg_dp_cape, 'kp-.','LineWidth',3,'MarkerSize',10); hold on

axis([N_all(1)/10 N_all(end)*10 1e-4 0.018])
xticks([1e4, 1e5, 1e6])
xticklabels({'10k','100k','1000k'})
set(gca,'FontSize',FS,'FontWeight','bold')
xlabel('(c) Total samples (N)','FontSize',FS,'FontWeight','bold');
ylabel('s_1','FontSize',FS,'FontWeight','bold');
title('Synth (\epsilon = 0.5)','FontSize',FS,'FontWeight','bold')

% for s2
subplot(144);
semilogx(N_all, score2_avg_np_pool, 'ro--','LineWidth',3,'MarkerSize',10); hold on
semilogx(N_all, score2_avg_dp_local, 'r*--','LineWidth',3,'MarkerSize',10); hold on
semilogx(N_all, score2_avg_dp_conv1, 'bsq:','LineWidth',3,'MarkerSize',10); hold on
semilogx(N_all, score2_avg_dp_cape, 'kp-.','LineWidth',3,'MarkerSize',10); hold on

axis([N_all(1)/10 N_all(end)*10 0.01 0.08])
xticks([1e4, 1e5, 1e6])
xticklabels({'10k','100k','1000k'})
set(gca,'FontSize',FS,'FontWeight','bold')
xlabel('(d) Total samples (N)','FontSize',FS,'FontWeight','bold');
ylabel('s_2','FontSize',FS,'FontWeight','bold');
title('Synth (\epsilon = 0.5)','FontSize',FS,'FontWeight','bold')


 
%% for MNIST and XRMB vs all
figure

% for CHindex vs epsilon
load('results_vs_eps_MNIST_d100_K50_N30k.mat')
subplot(261);
ids = 1:7;

CHindex_np_pool = repmat(CHindex_np_pool, 1, length(epsilon_all));
CHindex_avg_np_pool = mean(CHindex_np_pool(:,ids));

CHindex_avg_dp_local = mean(CHindex_dp_local(:,ids));

CHindex_avg_dp_conv1 = mean(CHindex_dp_conv1(:,ids));

CHindex_avg_dp_cape = mean(CHindex_dp_cape(:,ids));

semilogx(epsilon_all, CHindex_avg_np_pool, 'ro--','LineWidth',3,'MarkerSize',10); hold on
semilogx(epsilon_all, CHindex_avg_dp_local, 'r*--','LineWidth',3,'MarkerSize',10); hold on
semilogx(epsilon_all(2:end), CHindex_avg_dp_conv1(2:end), 'bsq:','LineWidth',3,'MarkerSize',10); hold on
semilogx(epsilon_all(2:end), CHindex_avg_dp_cape(2:end), 'kp-.','LineWidth',3,'MarkerSize',10); hold on

xticks([1e-3, 1e-2, 1e-1, 10])
xticklabels({'1e^{-3}','1e^{-2}','1e^{-1}', '10'})
axis([epsilon_all(1)/10 epsilon_all(end)*10 1e-4 2e4])
set(gca,'FontSize',FS,'FontWeight','bold')
xlabel('(a) Privacy param (\epsilon)','FontSize',FS,'FontWeight','bold');
ylabel('CHindex','FontSize',FS,'FontWeight','bold');
title('MNIST (N = 30k)','FontSize',FS,'FontWeight','bold')
legend('non-priv','local','conv','capeCCA','Location','best')

load('results_vs_eps_MNIST_d100_K50_N50k.mat')
subplot(262);

CHindex_np_pool = repmat(CHindex_np_pool, 1, length(epsilon_all));
CHindex_avg_np_pool = mean(CHindex_np_pool(:,ids));

CHindex_avg_dp_local = mean(CHindex_dp_local(:,ids));

CHindex_avg_dp_conv1 = mean(CHindex_dp_conv1(:,ids));

CHindex_avg_dp_cape = mean(CHindex_dp_cape(:,ids));

semilogx(epsilon_all, CHindex_avg_np_pool, 'ro--','LineWidth',3,'MarkerSize',10); hold on
semilogx(epsilon_all, CHindex_avg_dp_local, 'r*--','LineWidth',3,'MarkerSize',10); hold on
semilogx(epsilon_all(2:end), CHindex_avg_dp_conv1(2:end), 'bsq:','LineWidth',3,'MarkerSize',10); hold on
semilogx(epsilon_all(2:end), CHindex_avg_dp_cape(2:end), 'kp-.','LineWidth',3,'MarkerSize',10); hold on

xticks([1e-3, 1e-2, 1e-1, 10])
xticklabels({'1e^{-3}','1e^{-2}','1e^{-1}', '10'})
axis([epsilon_all(1)/10 epsilon_all(end)*10 1e-4 3e4])
set(gca,'FontSize',FS,'FontWeight','bold')
xlabel('(b) Privacy param (\epsilon)','FontSize',FS,'FontWeight','bold');
ylabel('CHindex','FontSize',FS,'FontWeight','bold');
title('MNIST (N = 50k)','FontSize',FS,'FontWeight','bold')

load('results_vs_eps_XRMB_d50_K20_p30.mat')
subplot(263);

CHindex_np_pool = repmat(CHindex_np_pool, 1, length(epsilon_all));
CHindex_avg_np_pool = mean(CHindex_np_pool(:,ids));

CHindex_avg_dp_local = mean(CHindex_dp_local(:,ids));

CHindex_avg_dp_conv1 = mean(CHindex_dp_conv1(:,ids));

CHindex_avg_dp_cape = mean(CHindex_dp_cape(:,ids));

semilogx(epsilon_all, CHindex_avg_np_pool, 'ro--','LineWidth',3,'MarkerSize',10); hold on
semilogx(epsilon_all, CHindex_avg_dp_local, 'r*--','LineWidth',3,'MarkerSize',10); hold on
semilogx(epsilon_all(2:end), CHindex_avg_dp_conv1(2:end), 'bsq:','LineWidth',3,'MarkerSize',10); hold on
semilogx(epsilon_all(2:end), CHindex_avg_dp_cape(2:end), 'kp-.','LineWidth',3,'MarkerSize',10); hold on

xticks([1e-3, 1e-2, 1e-1, 10])
xticklabels({'1e^{-3}','1e^{-2}','1e^{-1}', '10'})
axis([epsilon_all(1)/10 epsilon_all(end)*10 1e-4 2.5e4])
set(gca,'FontSize',FS,'FontWeight','bold')
xlabel('(c) Privacy param (\epsilon)','FontSize',FS,'FontWeight','bold');
ylabel('CHindex','FontSize',FS,'FontWeight','bold');
title('XRMB (p = 30)','FontSize',FS,'FontWeight','bold')

load('results_vs_eps_XRMB_d50_K20_p50.mat')
subplot(264);

CHindex_np_pool = repmat(CHindex_np_pool, 1, length(epsilon_all));
CHindex_avg_np_pool = mean(CHindex_np_pool(:,ids));

CHindex_avg_dp_local = mean(CHindex_dp_local(:,ids));

CHindex_avg_dp_conv1 = mean(CHindex_dp_conv1(:,ids));

CHindex_avg_dp_cape = mean(CHindex_dp_cape(:,ids));

semilogx(epsilon_all, CHindex_avg_np_pool, 'ro--','LineWidth',3,'MarkerSize',10); hold on
semilogx(epsilon_all, CHindex_avg_dp_local, 'r*--','LineWidth',3,'MarkerSize',10); hold on
semilogx(epsilon_all(2:end), CHindex_avg_dp_conv1(2:end), 'bsq:','LineWidth',3,'MarkerSize',10); hold on
semilogx(epsilon_all(2:end), CHindex_avg_dp_cape(2:end), 'kp-.','LineWidth',3,'MarkerSize',10); hold on

xticks([1e-3, 1e-2, 1e-1, 10])
xticklabels({'1e^{-3}','1e^{-2}','1e^{-1}', '10'})
axis([epsilon_all(1)/10 epsilon_all(end)*10 1e-4 4.5e4])
set(gca,'FontSize',FS,'FontWeight','bold')
xlabel('(d) Privacy param (\epsilon)','FontSize',FS,'FontWeight','bold');
ylabel('CHindex','FontSize',FS,'FontWeight','bold');
title('XRMB (p = 50)','FontSize',FS,'FontWeight','bold')

% for err_corr
subplot(265)
load results_vs_eps_synth_D5_N500
ids = 1:length(epsilon_all);

err_np_pool = repmat(err_np_pool, 1, length(epsilon_all));
err_avg_np_pool = mean(err_np_pool(:,ids));

err_avg_dp_local = mean(err_dp_local(:,ids));

err_avg_dp_conv1 = mean(err_dp_conv1(:,ids));

err_avg_dp_cape = mean(err_dp_cape(:,ids));

semilogx(epsilon_all, err_avg_np_pool, 'ro--','LineWidth',3,'MarkerSize',10); hold on
semilogx(epsilon_all, err_avg_dp_local, 'r*--','LineWidth',3,'MarkerSize',10); hold on
semilogx(epsilon_all, err_avg_dp_conv1, 'bsq:','LineWidth',3,'MarkerSize',10); hold on
semilogx(epsilon_all, err_avg_dp_cape, 'kp-.','LineWidth',3,'MarkerSize',10); hold on

xticks([1e-3, 1e-2, 1e-1, 10])
xticklabels({'1e^{-3}','1e^{-2}','1e^{-1}', '10'})
axis([epsilon_all(1)/10 epsilon_all(end)*10 1e-4 0.2])
set(gca,'FontSize',FS,'FontWeight','bold')
xlabel('(e) Privacy param (\epsilon)','FontSize',FS,'FontWeight','bold');
ylabel('err_{corr}','FontSize',FS,'FontWeight','bold');
title('fMRI+EEG (N = 500)','FontSize',FS,'FontWeight','bold')

subplot(266)
load results_vs_eps_synth_D5_N1000
ids = 1:length(epsilon_all);

err_np_pool = repmat(err_np_pool, 1, length(epsilon_all));
err_avg_np_pool = mean(err_np_pool(:,ids));

err_avg_dp_local = mean(err_dp_local(:,ids));

err_avg_dp_conv1 = mean(err_dp_conv1(:,ids));

err_avg_dp_cape = mean(err_dp_cape(:,ids));

semilogx(epsilon_all, err_avg_np_pool, 'ro--','LineWidth',3,'MarkerSize',10); hold on
semilogx(epsilon_all, err_avg_dp_local, 'r*--','LineWidth',3,'MarkerSize',10); hold on
semilogx(epsilon_all, err_avg_dp_conv1, 'bsq:','LineWidth',3,'MarkerSize',10); hold on
semilogx(epsilon_all, err_avg_dp_cape, 'kp-.','LineWidth',3,'MarkerSize',10); hold on

xticks([1e-3, 1e-2, 1e-1, 10])
xticklabels({'1e^{-3}','1e^{-2}','1e^{-1}', '10'})
axis([epsilon_all(1)/10 epsilon_all(end)*10 1e-4 0.2])
set(gca,'FontSize',FS,'FontWeight','bold')
xlabel('(f) Privacy param (\epsilon)','FontSize',FS,'FontWeight','bold');
ylabel('err_{corr}','FontSize',FS,'FontWeight','bold');
title('fMRI+EEG (N = 1000)','FontSize',FS,'FontWeight','bold')


% for CHindex vs samples
load('results_vs_samples_mnist_d100_K50_eps_0_05.mat')
subplot(267);
ids = 1:length(N_all);

CHindex_np_pool = repmat(CHindex_np_pool, 1, length(N_all));
CHindex_avg_np_pool = mean(CHindex_np_pool(:,ids));

CHindex_avg_dp_local = mean(CHindex_dp_local(:,ids));

CHindex_avg_dp_conv1 = mean(CHindex_dp_conv1(:,ids));

CHindex_avg_dp_cape = mean(CHindex_dp_cape(:,ids));

semilogx(N_all, CHindex_avg_np_pool, 'ro--','LineWidth',3,'MarkerSize',10); hold on
semilogx(N_all, CHindex_avg_dp_local, 'r*--','LineWidth',3,'MarkerSize',10); hold on
semilogx(N_all, CHindex_avg_dp_conv1, 'bsq:','LineWidth',3,'MarkerSize',10); hold on
semilogx(N_all, CHindex_avg_dp_cape, 'kp-.','LineWidth',3,'MarkerSize',10); hold on

axis([N_all(1)/2 N_all(end)*2 1e-4 3e4])
xticks([10e3, 30e3, 50e3])
xticklabels({'10k','30k', '50k'})
set(gca,'FontSize',FS,'FontWeight','bold')
xlabel('(g) Total samples (N)','FontSize',FS,'FontWeight','bold');
ylabel('CHindex','FontSize',FS,'FontWeight','bold');
title('MNIST (\epsilon = 0.05)','FontSize',FS,'FontWeight','bold')
legend('non-priv','local','conv','capeCCA','Location','best')

load('results_vs_samples_mnist_d100_K50_eps_0_1.mat')
subplot(268);
ids = 1:length(N_all);

CHindex_np_pool = repmat(CHindex_np_pool, 1, length(N_all));
CHindex_avg_np_pool = mean(CHindex_np_pool(:,ids));

CHindex_avg_dp_local = mean(CHindex_dp_local(:,ids));

CHindex_avg_dp_conv1 = mean(CHindex_dp_conv1(:,ids));

CHindex_avg_dp_cape = mean(CHindex_dp_cape(:,ids));

semilogx(N_all, CHindex_avg_np_pool, 'ro--','LineWidth',3,'MarkerSize',10); hold on
semilogx(N_all, CHindex_avg_dp_local, 'r*--','LineWidth',3,'MarkerSize',10); hold on
semilogx(N_all, CHindex_avg_dp_conv1, 'bsq:','LineWidth',3,'MarkerSize',10); hold on
semilogx(N_all, CHindex_avg_dp_cape, 'kp-.','LineWidth',3,'MarkerSize',10); hold on

axis([N_all(1)/2 N_all(end)*2 1e-4 3e4])
xticks([10e3, 30e3, 50e3])
xticklabels({'10k','30k', '50k'})
set(gca,'FontSize',FS,'FontWeight','bold')
xlabel('(h) Total samples (N)','FontSize',FS,'FontWeight','bold');
ylabel('CHindex','FontSize',FS,'FontWeight','bold');
title('MNIST (\epsilon = 0.1)','FontSize',FS,'FontWeight','bold')

load('results_vs_samples_xrmb_d50_K20_eps_0_1.mat')
ids = 1:length(p_all);
subplot(269);

CHindex_np_pool = repmat(CHindex_np_pool, 1, length(p_all));
CHindex_avg_np_pool = mean(CHindex_np_pool(:,ids));

CHindex_avg_dp_local = mean(CHindex_dp_local(:,ids));

CHindex_avg_dp_conv1 = mean(CHindex_dp_conv1(:,ids));

CHindex_avg_dp_cape = mean(CHindex_dp_cape(:,ids));

semilogx(p_all, CHindex_avg_np_pool, 'ro--','LineWidth',3,'MarkerSize',10); hold on
semilogx(p_all, CHindex_avg_dp_local, 'r*--','LineWidth',3,'MarkerSize',10); hold on
semilogx(p_all, CHindex_avg_dp_conv1, 'bsq:','LineWidth',3,'MarkerSize',10); hold on
semilogx(p_all, CHindex_avg_dp_cape, 'kp-.','LineWidth',3,'MarkerSize',10); hold on

axis([p_all(1)/2 p_all(end)*2 1e-4 4e4])
xticks([10, 20, 30, 50])
xticklabels({'10','20','30', '50'})
set(gca,'FontSize',FS,'FontWeight','bold')
xlabel('(i) Replication param (p)','FontSize',FS,'FontWeight','bold');
ylabel('CHindex','FontSize',FS,'FontWeight','bold');
title('XRMB (\epsilon = 0.1)','FontSize',FS,'FontWeight','bold')


load('results_vs_samples_xrmb_d50_K20_eps_0_5.mat')
ids = 1:length(p_all);
subplot(2,6,10);

CHindex_np_pool = repmat(CHindex_np_pool, 1, length(p_all));
CHindex_avg_np_pool = mean(CHindex_np_pool(:,ids));

CHindex_avg_dp_local = mean(CHindex_dp_local(:,ids));

CHindex_avg_dp_conv1 = mean(CHindex_dp_conv1(:,ids));

CHindex_avg_dp_cape = mean(CHindex_dp_cape(:,ids));

semilogx(p_all, CHindex_avg_np_pool, 'ro--','LineWidth',3,'MarkerSize',10); hold on
semilogx(p_all, CHindex_avg_dp_local, 'r*--','LineWidth',3,'MarkerSize',10); hold on
semilogx(p_all, CHindex_avg_dp_conv1, 'bsq:','LineWidth',3,'MarkerSize',10); hold on
semilogx(p_all, CHindex_avg_dp_cape, 'kp-.','LineWidth',3,'MarkerSize',10); hold on

axis([p_all(1)/2 p_all(end)*2 1e-4 4e4])
xticks([10, 20, 30, 50])
xticklabels({'10','20','30', '50'})
set(gca,'FontSize',FS,'FontWeight','bold')
xlabel('(j) Replication param (p)','FontSize',FS,'FontWeight','bold');
ylabel('CHindex','FontSize',FS,'FontWeight','bold');
title('XRMB (\epsilon = 0.5)','FontSize',FS,'FontWeight','bold')

% for err_corr
load results_vs_samples_synth_D5_eps_0_5
ids = 1:length(N_all);

subplot(2,6,11)

err_avg_np_pool = mean(err_np_pool(:,ids));

err_avg_dp_local = mean(err_dp_local(:,ids));

err_avg_dp_conv1 = mean(err_dp_conv1(:,ids));

err_avg_dp_cape = mean(err_dp_cape(:,ids));

semilogx(N_all, err_avg_np_pool, 'ro--','LineWidth',3,'MarkerSize',10); hold on
semilogx(N_all, err_avg_dp_local, 'r*--','LineWidth',3,'MarkerSize',10); hold on
semilogx(N_all, err_avg_dp_conv1, 'bsq:','LineWidth',3,'MarkerSize',10); hold on
semilogx(N_all, err_avg_dp_cape, 'kp-.','LineWidth',3,'MarkerSize',10); hold on

% axis([N_all(idsDP(1))/10 N_all(end)*10 1e-4 1])
set(gca,'FontSize',FS,'FontWeight','bold')
xlabel('(k) Total samples (N)','FontSize',FS,'FontWeight','bold');
ylabel('err_{corr}','FontSize',FS,'FontWeight','bold');
title('fMRI+EEG (\epsilon = 0.5)','FontSize',FS,'FontWeight','bold')

load results_vs_samples_synth_D5_eps_1
ids = 1:length(N_all);

subplot(2,6,12)

err_avg_np_pool = mean(err_np_pool(:,ids));

err_avg_dp_local = mean(err_dp_local(:,ids));

err_avg_dp_conv1 = mean(err_dp_conv1(:,ids));

err_avg_dp_cape = mean(err_dp_cape(:,ids));

semilogx(N_all, err_avg_np_pool, 'ro--','LineWidth',3,'MarkerSize',10); hold on
semilogx(N_all, err_avg_dp_local, 'r*--','LineWidth',3,'MarkerSize',10); hold on
semilogx(N_all, err_avg_dp_conv1, 'bsq:','LineWidth',3,'MarkerSize',10); hold on
semilogx(N_all, err_avg_dp_cape, 'kp-.','LineWidth',3,'MarkerSize',10); hold on

% axis([N_all(idsDP(1))/10 N_all(end)*10 1e-4 1])
set(gca,'FontSize',FS,'FontWeight','bold')
xlabel('(l) Total samples (N)','FontSize',FS,'FontWeight','bold');
ylabel('err_{corr}','FontSize',FS,'FontWeight','bold');
title('fMRI+EEG (\epsilon = 1.0)','FontSize',FS,'FontWeight','bold')
