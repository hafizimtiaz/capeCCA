clear all;clc;close all
ids = [1:6];

FS = 24;


%% for d = 100 - 200k samples
load results_vs_delta_synth_d_100_K_20_N_500k.mat
ids = 1:6;

score1CCA = repmat(score1CCA,1,length(delta_all));
score1avgCCA = mean(score1CCA(:,ids));
score1avgDPCCAAG = mean(score1DPCCAAG(:,ids));

figure(1) % for s1 vs epsilon
subplot(231);
loglog(delta_all(ids),score1avgCCA(ids),'gs:','LineWidth',6,'MarkerSize',10); hold on
loglog(delta_all(ids),score1avgDPCCAAG(ids),'rs:','LineWidth',6,'MarkerSize',10); hold on
% axis([1e-4 10 5e-5 1e-3])

set(gca,'FontSize',FS,'FontWeight','bold')
xlabel('\delta','FontSize',FS,'FontWeight','bold');
ylabel('s_1','FontSize',FS,'FontWeight','bold');
title('D_x = D_y = 100, N = 500k, \epsilon = 1','FontSize',FS,'FontWeight','bold')
legend('CCA','DPCCAG','Location','NE')

figure(1) % for s2 vs epsilon
subplot(232)
ids = [1:6];

score2CCA = repmat(score2CCA,1,length(ids));
score2avgCCA = mean(score2CCA(:,ids));
score2avgDPCCAAG = mean(score2DPCCAAG(:,ids));

loglog(delta_all(ids),score2avgCCA(ids),'gs:','LineWidth',6,'MarkerSize',10); hold on
loglog(delta_all(ids),score2avgDPCCAAG(ids),'rs:','LineWidth',6,'MarkerSize',10); hold on
% axis([1e-2 1e2 1e-5 1])

set(gca,'FontSize',FS,'FontWeight','bold')
xlabel('\delta','FontSize',FS,'FontWeight','bold');
ylabel('s_2','FontSize',FS,'FontWeight','bold');
title('D_x = D_y = 100, N = 500k, \epsilon = 1','FontSize',FS,'FontWeight','bold')
legend('CCA','DPCCAG','Location','SE')


figure(1) % for CHindex vs epsilon
subplot(233)

CHindexCCA = repmat(CHindexCCA,1,length(ids));
CHavgCCA = mean(CHindexCCA(:,ids));
CHavgDPCCAAG = mean(CHindexDPCCAAG(:,ids));


loglog(delta_all(ids),CHavgCCA(ids),'gs:','LineWidth',6,'MarkerSize',10); hold on
loglog(delta_all(ids),CHavgDPCCAAG(ids),'rs:','LineWidth',6,'MarkerSize',10); hold on
% axis([1e-3 1e2 1e2 1e6])

set(gca,'FontSize',FS,'FontWeight','bold')
xlabel('\delta','FontSize',FS,'FontWeight','bold');
ylabel('CH Index','FontSize',FS,'FontWeight','bold');
title('D_x = D_y = 100, N = 500k, \epsilon = 1','FontSize',FS,'FontWeight','bold')
% legend('CCA','PCA','DP-CCA','DP-PCA', 'Location','NW')
legend('CCA','DPCCAG','Location','SE')


%% for d = 100 - 800k samples
load results_vs_delta_synth_d_100_K_20_N_800k.mat
ids = 1:6;

score1CCA = repmat(score1CCA,1,length(delta_all));
score1avgCCA = mean(score1CCA(:,ids));
score1avgDPCCAAG = mean(score1DPCCAAG(:,ids));

figure(1) % for s1 vs epsilon
subplot(234);
loglog(delta_all(ids),score1avgCCA(ids),'gs:','LineWidth',6,'MarkerSize',10); hold on
loglog(delta_all(ids),score1avgDPCCAAG(ids),'rs:','LineWidth',6,'MarkerSize',10); hold on
% axis([1e-4 10 5e-5 1e-3])

set(gca,'FontSize',FS,'FontWeight','bold')
xlabel('\delta','FontSize',FS,'FontWeight','bold');
ylabel('s_1','FontSize',FS,'FontWeight','bold');
title('D_x = D_y = 100, N = 800k, \epsilon = 2','FontSize',FS,'FontWeight','bold')
legend('CCA','DPCCAG','Location','NE')

figure(1) % for s2 vs epsilon
subplot(235)
ids = [1:6];

score2CCA = repmat(score2CCA,1,length(ids));
score2avgCCA = mean(score2CCA(:,ids));
score2avgDPCCAAG = mean(score2DPCCAAG(:,ids));

loglog(delta_all(ids),score2avgCCA(ids),'gs:','LineWidth',6,'MarkerSize',10); hold on
loglog(delta_all(ids),score2avgDPCCAAG(ids),'rs:','LineWidth',6,'MarkerSize',10); hold on
% axis([1e-2 1e2 1e-5 1])

set(gca,'FontSize',FS,'FontWeight','bold')
xlabel('\delta','FontSize',FS,'FontWeight','bold');
ylabel('s_2','FontSize',FS,'FontWeight','bold');
title('D_x = D_y = 100, N = 800k, \epsilon = 2','FontSize',FS,'FontWeight','bold')
legend('CCA','DPCCAG','Location','SE')


figure(1) % for CHindex vs epsilon
subplot(236)

CHindexCCA = repmat(CHindexCCA,1,length(ids));
CHavgCCA = mean(CHindexCCA(:,ids));
CHavgDPCCAAG = mean(CHindexDPCCAAG(:,ids));


loglog(delta_all(ids),CHavgCCA(ids),'gs:','LineWidth',6,'MarkerSize',10); hold on
loglog(delta_all(ids),CHavgDPCCAAG(ids),'rs:','LineWidth',6,'MarkerSize',10); hold on
% axis([1e-3 1e2 1e2 1e6])

set(gca,'FontSize',FS,'FontWeight','bold')
xlabel('\delta','FontSize',FS,'FontWeight','bold');
ylabel('CH Index','FontSize',FS,'FontWeight','bold');
title('D_x = D_y = 100, N = 800k, \epsilon = 2','FontSize',FS,'FontWeight','bold')
% legend('CCA','PCA','DP-CCA','DP-PCA', 'Location','NW')
legend('CCA','DPCCAG','Location','SE')
