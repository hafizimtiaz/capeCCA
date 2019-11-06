clear all;clc;close all
ids = [1:6];

FS = 24;

%% for d = 100 - eps = 2
load results_vs_samples_synth_d_100_K_20_eps_2.mat

score1avgCCA = mean(score1CCA(:,ids));
score2avgCCA = mean(score2CCA(:,ids));
score1avgDPCCAAG = mean(score1DPCCAAG(:,ids));
score2avgDPCCAAG = mean(score2DPCCAAG(:,ids));

figure(1)
subplot(231);
ids_AG = [1:4];
loglog(N_all(ids_AG),score1avgCCA(ids_AG),'gs:','LineWidth',6,'MarkerSize',10); hold on
% loglog(N_all,STavgPCA,'rs:','LineWidth',3,'MarkerSize',10); hold on
% loglog(N_all,score1avgDPCCA,'ks:','LineWidth',6,'MarkerSize',10); hold on
% loglog(N_all,STavgDPPCA,'bs:','LineWidth',3,'MarkerSize',10); hold on

loglog(N_all(ids_AG),score1avgDPCCAAG(ids_AG),'rs:','LineWidth',6,'MarkerSize',10); hold on

axis([1e4 5e5 1e-5 1e-3])
set(gca,'FontSize',FS,'FontWeight','bold')
xlabel('No. of samples','FontSize',FS,'FontWeight','bold');
ylabel('s_1','FontSize',FS,'FontWeight','bold');
title('D_x = D_y = 100, \epsilon = 2, \delta = 0.01','FontSize',FS,'FontWeight','bold')
% legend('CCA','PCA','DP-CCA','DP-PCA', 'Location','NW')
legend('CCA','DPCCAG','Location','NE')

figure(1)
subplot(232);
loglog(N_all,score2avgCCA,'gs:','LineWidth',6,'MarkerSize',10); hold on
% loglog(N_all,STavgPCA,'rs:','LineWidth',3,'MarkerSize',10); hold on
% loglog(N_all,score2avgDPCCA,'ks:','LineWidth',6,'MarkerSize',10); hold on
% loglog(N_all,STavgDPPCA,'bs:','LineWidth',3,'MarkerSize',10); hold on
loglog(N_all,score2avgDPCCAAG,'rs:','LineWidth',6,'MarkerSize',10); hold on

set(gca,'FontSize',FS,'FontWeight','bold')
xlabel('No. of samples','FontSize',FS,'FontWeight','bold');
ylabel('s_2','FontSize',FS,'FontWeight','bold');
title('D_x = D_y = 100, \epsilon = 2, \delta = 0.01','FontSize',FS,'FontWeight','bold')
% legend('CCA','PCA','DP-CCA','DP-PCA', 'Location','NW')
legend('CCA','DPCCAG','Location','SE')


CHavgCCA = mean(CHindexCCA(:,ids));
CHavgPCA = mean(CHindexPCA(:,ids));

CHavgDPCCAAG = mean(CHindexDPCCAAG(:,ids));
CHavgDPPCA = mean(CHindexDPPCA(:,ids));

figure(1)
subplot(233)
% subplot(121);
loglog(N_all,CHavgCCA,'gs:','LineWidth',6,'MarkerSize',10); hold on
% loglog(N_all,CHavgPCA,'rs:','LineWidth',3,'MarkerSize',10); hold on
% loglog(N_all,CHavgDPCCA,'ks:','LineWidth',6,'MarkerSize',10); hold on
% loglog(N_all,CHavgDPPCA,'bs:','LineWidth',3,'MarkerSize',10); hold on
loglog(N_all,CHavgDPCCAAG,'rs:','LineWidth',6,'MarkerSize',10); hold on

set(gca,'FontSize',FS,'FontWeight','bold')
xlabel('No. of samples','FontSize',FS,'FontWeight','bold');
ylabel('CH Index','FontSize',FS,'FontWeight','bold');
title('D_x = D_y = 100, \epsilon = 2, \delta = 0.01','FontSize',FS,'FontWeight','bold')
% legend('CCA','PCA','DP-CCA','DP-PCA', 'Location','NW')
legend('CCA','DPCCAG','Location','NW')

%% for d = 100 - eps = 5
load results_vs_samples_synth_d_100_K_20_eps_5.mat

score1avgCCA = mean(score1CCA(:,ids));
score2avgCCA = mean(score2CCA(:,ids));

score1avgDPCCAAG = mean(score1DPCCAAG(:,ids));
score2avgDPCCAAG = mean(score2DPCCAAG(:,ids));

figure(1)
subplot(234);
loglog(N_all,score1avgCCA,'gs:','LineWidth',6,'MarkerSize',10); hold on
% loglog(N_all,STavgPCA,'rs:','LineWidth',3,'MarkerSize',10); hold on
% loglog(N_all,score1avgDPCCA,'ks:','LineWidth',6,'MarkerSize',10); hold on
% loglog(N_all,STavgDPPCA,'bs:','LineWidth',3,'MarkerSize',10); hold on
ids = [1:3,5:6];
loglog(N_all(ids),score1avgDPCCAAG(ids),'rs:','LineWidth',6,'MarkerSize',10); hold on

set(gca,'FontSize',FS,'FontWeight','bold')
xlabel('No. of samples','FontSize',FS,'FontWeight','bold');
ylabel('s_1','FontSize',FS,'FontWeight','bold');
title('D_x = D_y = 100, \epsilon = 5, \delta = 0.01','FontSize',FS,'FontWeight','bold')
% legend('CCA','PCA','DP-CCA','DP-PCA', 'Location','NW')
legend('CCA','DPCCAG','Location','NE')

figure(1)
subplot(235);
loglog(N_all,score2avgCCA,'gs:','LineWidth',6,'MarkerSize',10); hold on
% loglog(N_all,STavgPCA,'rs:','LineWidth',3,'MarkerSize',10); hold on
% loglog(N_all,score2avgDPCCA,'ks:','LineWidth',6,'MarkerSize',10); hold on
% loglog(N_all,STavgDPPCA,'bs:','LineWidth',3,'MarkerSize',10); hold on
loglog(N_all,score2avgDPCCAAG,'rs:','LineWidth',6,'MarkerSize',10); hold on

set(gca,'FontSize',FS,'FontWeight','bold')
xlabel('No. of samples','FontSize',FS,'FontWeight','bold');
ylabel('s_2','FontSize',FS,'FontWeight','bold');
title('D_x = D_y = 100, \epsilon = 5, \delta = 0.01','FontSize',FS,'FontWeight','bold')
% legend('CCA','PCA','DP-CCA','DP-PCA', 'Location','NW')
legend('CCA','DPCCAG','Location','SE')

ids = 1:6;
CHavgCCA = mean(CHindexCCA(:,ids));
CHavgPCA = mean(CHindexPCA(:,ids));

CHavgDPCCAAG = mean(CHindexDPCCAAG(:,ids));
CHavgDPPCA = mean(CHindexDPPCA(:,ids));

figure(1)
subplot(236)
% subplot(121);
loglog(N_all,CHavgCCA,'gs:','LineWidth',6,'MarkerSize',10); hold on
% loglog(N_all,CHavgPCA,'rs:','LineWidth',3,'MarkerSize',10); hold on
% loglog(N_all,CHavgDPCCA,'ks:','LineWidth',6,'MarkerSize',10); hold on
% loglog(N_all,CHavgDPPCA,'bs:','LineWidth',3,'MarkerSize',10); hold on
loglog(N_all,CHavgDPCCAAG,'rs:','LineWidth',6,'MarkerSize',10); hold on

set(gca,'FontSize',FS,'FontWeight','bold')
xlabel('No. of samples','FontSize',FS,'FontWeight','bold');
ylabel('CH Index','FontSize',FS,'FontWeight','bold');
title('D_x = D_y = 100, \epsilon = 5, \delta = 0.01','FontSize',FS,'FontWeight','bold')
% legend('CCA','PCA','DP-CCA','DP-PCA', 'Location','NW')
legend('CCA','DPCCAG','Location','NW')