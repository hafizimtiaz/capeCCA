clear all;clc;close all
ids = [1:6];

FS = 24;


%% for p = 30
load results_vs_delta_XRMB_d50_K20_p30.mat
ids = 1:6;

figure(1) % for CHindex vs epsilon
subplot(231)

CHindexCCA = repmat(CHindexCCA,1,length(ids));
CHavgCCA = mean(CHindexCCA(:,ids));
CHavgDPCCAAG = mean(CHindexDPCCAAG(:,ids));


loglog(delta_all(ids),CHavgCCA(ids),'gs:','LineWidth',6,'MarkerSize',10); hold on
loglog(delta_all(ids),CHavgDPCCAAG(ids),'rs:','LineWidth',6,'MarkerSize',10); hold on
% axis([1e-3 1e2 1e2 1e6])

set(gca,'FontSize',FS,'FontWeight','bold')
xlabel('\delta','FontSize',FS,'FontWeight','bold');
ylabel('CH Index','FontSize',FS,'FontWeight','bold');
title('XRMB data, p = 30, \epsilon = 1','FontSize',FS,'FontWeight','bold')
% legend('CCA','PCA','DP-CCA','DP-PCA', 'Location','NW')
legend('CCA','DPCCAG','Location','SE')

figure(1) % for avgErr vs epsilon
subplot(232)
ids = 1:6;

perc_errCCA = repmat(perc_errCCA',1,length(ids));
PEavgCCA = mean(perc_errCCA(:,ids));
PEavgDPCCAAG = mean(perc_errDPCCAAG(:,ids));

loglog(delta_all(ids),PEavgCCA(ids),'gs:','LineWidth',6,'MarkerSize',10); hold on
loglog(delta_all(ids),PEavgDPCCAAG(ids),'rs:','LineWidth',6,'MarkerSize',10); hold on
% axis([1e-3 1e2 1e-1 1e-0])

set(gca,'FontSize',FS,'FontWeight','bold')
xlabel('\delta','FontSize',FS,'FontWeight','bold');
ylabel('avgErr','FontSize',FS,'FontWeight','bold');
title('XRMB data, p = 30, \epsilon = 1','FontSize',FS,'FontWeight','bold')
legend('CCA','DPCCAG','Location','SE')

figure(1) % for avgErr vs epsilon
subplot(233)
ids = 1:6;

mutInfoCCA = repmat(mutInfoCCA',1,length(ids));
MIavgCCA = mean(mutInfoCCA(:,ids));
MIavgDPCCAAG = mean(mutInfoDPCCAAG(:,ids));

loglog(delta_all(ids),MIavgCCA(ids),'gs:','LineWidth',6,'MarkerSize',10); hold on
loglog(delta_all(ids),MIavgDPCCAAG(ids),'rs:','LineWidth',6,'MarkerSize',10); hold on
% axis([1e-3 1e2 1e-3 1e-1])

set(gca,'FontSize',FS,'FontWeight','bold')
xlabel('\delta','FontSize',FS,'FontWeight','bold');
ylabel('NMI','FontSize',FS,'FontWeight','bold');
title('XRMB data, p = 30, \epsilon = 1','FontSize',FS,'FontWeight','bold')
legend('CCA','DPCCAG','Location','SE')


%% for p = 50
load results_vs_delta_XRMB_d50_K20_p50.mat
ids = 1:6;

figure(1) % for CHindex vs epsilon
subplot(234)

CHindexCCA = repmat(CHindexCCA,1,length(ids));
CHavgCCA = mean(CHindexCCA(:,ids));
CHavgDPCCAAG = mean(CHindexDPCCAAG(:,ids));


loglog(delta_all(ids),CHavgCCA(ids),'gs:','LineWidth',6,'MarkerSize',10); hold on
loglog(delta_all(ids),CHavgDPCCAAG(ids),'rs:','LineWidth',6,'MarkerSize',10); hold on
% axis([1e-3 1e2 1e2 1e6])

set(gca,'FontSize',FS,'FontWeight','bold')
xlabel('\delta','FontSize',FS,'FontWeight','bold');
ylabel('CH Index','FontSize',FS,'FontWeight','bold');
title('XRMB data, p = 50, \epsilon = 0.5','FontSize',FS,'FontWeight','bold')
legend('CCA','DPCCAG','Location','SE')

figure(1) % for avgErr vs epsilon
subplot(235)
ids = 1:6;

perc_errCCA = repmat(perc_errCCA',1,length(ids));
PEavgCCA = mean(perc_errCCA(:,ids));
PEavgDPCCAAG = mean(perc_errDPCCAAG(:,ids));

loglog(delta_all(ids),PEavgCCA(ids),'gs:','LineWidth',6,'MarkerSize',10); hold on
loglog(delta_all(ids),PEavgDPCCAAG(ids),'rs:','LineWidth',6,'MarkerSize',10); hold on
% axis([1e-3 1e2 1e-1 1e-0])

set(gca,'FontSize',FS,'FontWeight','bold')
xlabel('\delta','FontSize',FS,'FontWeight','bold');
ylabel('avgErr','FontSize',FS,'FontWeight','bold');
title('XRMB data, p = 50, \epsilon = 0.5','FontSize',FS,'FontWeight','bold')
legend('CCA','DPCCAG','Location','SE')

figure(1) % for avgErr vs epsilon
subplot(236)
ids = 1:6;

mutInfoCCA = repmat(mutInfoCCA',1,length(ids));
MIavgCCA = mean(mutInfoCCA(:,ids));
MIavgDPCCAAG = mean(mutInfoDPCCAAG(:,ids));

loglog(delta_all(ids),MIavgCCA(ids),'gs:','LineWidth',6,'MarkerSize',10); hold on
loglog(delta_all(ids),MIavgDPCCAAG(ids),'rs:','LineWidth',6,'MarkerSize',10); hold on
% axis([1e-3 1e2 1e-3 1e-1])

set(gca,'FontSize',FS,'FontWeight','bold')
xlabel('\delta','FontSize',FS,'FontWeight','bold');
ylabel('NMI','FontSize',FS,'FontWeight','bold');
title('XRMB data, p = 50, \epsilon = 0.5','FontSize',FS,'FontWeight','bold')
legend('CCA','DPCCAG','Location','SE')

%% for paper
FS = 46;
load results_vs_delta_synth_d_100_K_20_N_800k.mat
figure(2)
subplot(141)
ids = [1:6];

score2CCA = repmat(score2CCA,1,length(ids));
score2avgCCA = mean(score2CCA(:,ids));
score2avgDPCCAAG = mean(score2DPCCAAG(:,ids));

loglog(delta_all(ids),score2avgCCA(ids),'bo-','LineWidth',6,'MarkerSize',10); hold on
loglog(delta_all(ids),score2avgDPCCAAG(ids),'ro:','LineWidth',6,'MarkerSize',10); hold on
% axis([1e-2 1e2 1e-5 1])

set(gca,'FontSize',FS,'FontWeight','bold')
xlabel('\delta','FontSize',FS,'FontWeight','bold');
ylabel('s_2','FontSize',FS,'FontWeight','bold');
title('Synthetic (N = 800k, \epsilon = 2)','FontSize',FS,'FontWeight','bold')
legend('CCA','DPCCAG','Location','SE')


figure(2) % for CHindex vs epsilon
subplot(142)

CHindexCCA = repmat(CHindexCCA,1,length(ids));
CHavgCCA = mean(CHindexCCA(:,ids));
CHavgDPCCAAG = mean(CHindexDPCCAAG(:,ids));


loglog(delta_all(ids),CHavgCCA(ids),'bo-','LineWidth',6,'MarkerSize',10); hold on
loglog(delta_all(ids),CHavgDPCCAAG(ids),'ro:','LineWidth',6,'MarkerSize',10); hold on
axis([0.5e-4 0.1 1e5 3.5e5])

set(gca,'FontSize',FS,'FontWeight','bold')
xlabel('\delta','FontSize',FS,'FontWeight','bold');
ylabel('CH Index','FontSize',FS,'FontWeight','bold');
title('Synthetic (N = 800k, \epsilon = 2)','FontSize',FS,'FontWeight','bold')
% legend('CCA','PCA','DP-CCA','DP-PCA', 'Location','NW')
legend('CCA','DPCCAG','Location','SE')

figure(2)
subplot(143)
load results_vs_delta_XRMB_d50_K20_p30.mat
ids = 1:5;

perc_errCCA = repmat(perc_errCCA',1,length(ids));
PEavgCCA = mean(perc_errCCA(:,ids));
PEavgDPCCAAG = mean(perc_errDPCCAAG(:,ids));

loglog(delta_all(ids),PEavgCCA(ids),'bo-','LineWidth',6,'MarkerSize',10); hold on
loglog(delta_all(ids),PEavgDPCCAAG(ids),'ro:','LineWidth',6,'MarkerSize',10); hold on
axis([0.5e-4 0.1 0.3 0.5])

set(gca,'FontSize',FS,'FontWeight','bold')
xlabel('\delta','FontSize',FS,'FontWeight','bold');
ylabel('avgErr','FontSize',FS,'FontWeight','bold');
title('XRMB (p = 30, \epsilon = 1)','FontSize',FS,'FontWeight','bold')
legend('CCA','DPCCAG','Location','SE')

figure(2) % for avgErr vs epsilon
subplot(144)
ids = 1:5;

mutInfoCCA = repmat(mutInfoCCA',1,length(ids));
MIavgCCA = mean(mutInfoCCA(:,ids));
MIavgDPCCAAG = mean(mutInfoDPCCAAG(:,ids));

loglog(delta_all(ids),MIavgCCA(ids),'bo-','LineWidth',6,'MarkerSize',10); hold on
loglog(delta_all(ids),MIavgDPCCAAG(ids),'ro:','LineWidth',6,'MarkerSize',10); hold on
axis([0.5e-4 0.1 0.03 0.05])

set(gca,'FontSize',FS,'FontWeight','bold')
xlabel('\delta','FontSize',FS,'FontWeight','bold');
ylabel('NMI','FontSize',FS,'FontWeight','bold');
title('XRMB (p = 30, \epsilon = 1)','FontSize',FS,'FontWeight','bold')
legend('CCA','DPCCAG','Location','SE')
