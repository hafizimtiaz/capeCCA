clear all;clc;close all
ids = [1:7];

FS = 24;

%% for d = 100 - 200k samples
load results_vs_eps_synth_d100_K20_N200k_v2.mat


ids = 1:7;
idsDP = 1:7;
score1CCA_pool = repmat(score1CCA_pool,1,length(epsilon_all));
score1avg_pool = mean(score1CCA_pool(:,ids));

score1CCA_dist = repmat(score1CCA_dist,1,length(epsilon_all));
score1avg_dist = mean(score1CCA_dist(:,ids));

score1avgDP_pool = mean(score1DPCCAAG_pool(:,ids));
score1avgDP_dist1 = mean(score1DPCCAAG_dist1(:,ids));
score1avgDP_dist2 = mean(score1DPCCAAG_dist2(:,ids));

figure(1) % for s1 vs epsilon
subplot(231);
loglog(epsilon_all(idsDP),score1avg_pool(idsDP),'ro-','LineWidth',3,'MarkerSize',10); hold on
loglog(epsilon_all(idsDP),score1avgDP_pool(idsDP),'bo-','LineWidth',3,'MarkerSize',10); hold on

loglog(epsilon_all(idsDP),score1avg_dist(idsDP),'r*-','LineWidth',3,'MarkerSize',10); hold on
loglog(epsilon_all(idsDP),score1avgDP_dist1(idsDP),'b*-','LineWidth',3,'MarkerSize',10); hold on
loglog(epsilon_all(idsDP),score1avgDP_dist2(idsDP),'k*-','LineWidth',3,'MarkerSize',10); hold on


set(gca,'FontSize',FS,'FontWeight','bold')
xlabel('\epsilon','FontSize',FS,'FontWeight','bold');
ylabel('s_1','FontSize',FS,'FontWeight','bold');
title('Synth (N = 200k, \delta = 0.01)','FontSize',FS,'FontWeight','bold')
legend('Pooled CCA','Pooled dpCCA','Dist CCA','Dist dpCCA1','Dist dpCCA2','Location','best')

figure(1) % for s2 vs epsilon
subplot(232)
ids = 1:7;
idsDP = 1:7;

score2CCA_pool = repmat(score2CCA_pool,1,length(epsilon_all));
score2avg_pool = mean(score2CCA_pool(:,ids));

score2CCA_dist = repmat(score2CCA_dist,1,length(epsilon_all));
score2avg_dist = mean(score2CCA_dist(:,ids));

score2avgDP_pool = mean(score2DPCCAAG_pool(:,ids));
score2avgDP_dist1 = mean(score2DPCCAAG_dist1(:,ids));
score2avgDP_dist2 = mean(score2DPCCAAG_dist2(:,ids));

loglog(epsilon_all(idsDP),score2avg_pool(idsDP),'ro-','LineWidth',3,'MarkerSize',10); hold on
loglog(epsilon_all(idsDP),score2avgDP_pool(idsDP),'bo-','LineWidth',3,'MarkerSize',10); hold on

loglog(epsilon_all(idsDP),score2avg_dist(idsDP),'r*-','LineWidth',3,'MarkerSize',10); hold on
loglog(epsilon_all(idsDP),score2avgDP_dist1(idsDP),'b*-','LineWidth',3,'MarkerSize',10); hold on
loglog(epsilon_all(idsDP),score2avgDP_dist2(idsDP),'k*-','LineWidth',3,'MarkerSize',10); hold on


set(gca,'FontSize',FS,'FontWeight','bold')
xlabel('\epsilon','FontSize',FS,'FontWeight','bold');
ylabel('s_2','FontSize',FS,'FontWeight','bold');
title('Synth (N = 200k, \delta = 0.01)','FontSize',FS,'FontWeight','bold')
legend('Pooled CCA','Pooled dpCCA','Dist CCA','Dist dpCCA1','Dist dpCCA2','Location','best')

% for CH index
ids = 1:7;
idsDP = 1:7;
CHindexCCA_pool = repmat(CHindexCCA_pool,1,length(epsilon_all));
CHindexavg_pool = mean(CHindexCCA_pool(:,ids));

CHindexCCA_dist = repmat(CHindexCCA_dist,1,length(epsilon_all));
CHindexavg_dist = mean(CHindexCCA_dist(:,ids));

CHindexavgDP_pool = mean(CHindexDPCCAAG_pool(:,ids));
CHindexavgDP_dist1 = mean(CHindexDPCCAAG_dist1(:,ids));
CHindexavgDP_dist2 = mean(CHindexDPCCAAG_dist2(:,ids));

figure(1)
subplot(233);
loglog(epsilon_all(idsDP),CHindexavg_pool(idsDP),'ro-','LineWidth',3,'MarkerSize',10); hold on
loglog(epsilon_all(idsDP),CHindexavgDP_pool(idsDP),'bo-','LineWidth',3,'MarkerSize',10); hold on

loglog(epsilon_all(ids),CHindexavg_dist(ids),'r*-','LineWidth',3,'MarkerSize',10); hold on
loglog(epsilon_all(idsDP),CHindexavgDP_dist1(idsDP),'b*-','LineWidth',3,'MarkerSize',10); hold on
loglog(epsilon_all(idsDP),CHindexavgDP_dist2(idsDP),'k*-','LineWidth',3,'MarkerSize',10); hold on

% axis([epsilon_all(idsDP(1))/10 epsilon_all(end)*10 1e4 1.8e4])
set(gca,'FontSize',FS,'FontWeight','bold')
xlabel('\epsilon','FontSize',FS,'FontWeight','bold');
ylabel('CH Index','FontSize',FS,'FontWeight','bold');
title('Synth (N = 200k, \delta = 0.01)','FontSize',FS,'FontWeight','bold')
legend('Pooled CCA','Pooled dpCCA','Dist dpCCA','Location','best')

%% for d = 100 - 800k samples
load results_vs_eps_synth_d100_K20_N800k_v2.mat


ids = 1:7;
idsDP = 1:7;
score1CCA_pool = repmat(score1CCA_pool,1,length(epsilon_all));
score1avg_pool = mean(score1CCA_pool(:,ids));

score1CCA_dist = repmat(score1CCA_dist,1,length(epsilon_all));
score1avg_dist = mean(score1CCA_dist(:,ids));

score1avgDP_pool = mean(score1DPCCAAG_pool(:,ids));
score1avgDP_dist1 = mean(score1DPCCAAG_dist1(:,ids));
score1avgDP_dist2 = mean(score1DPCCAAG_dist2(:,ids));

figure(1) % for s1 vs epsilon
subplot(234);
loglog(epsilon_all(idsDP),score1avg_pool(idsDP),'ro-','LineWidth',3,'MarkerSize',10); hold on
loglog(epsilon_all(idsDP),score1avgDP_pool(idsDP),'bo-','LineWidth',3,'MarkerSize',10); hold on

loglog(epsilon_all(idsDP),score1avg_dist(idsDP),'r*-','LineWidth',3,'MarkerSize',10); hold on
loglog(epsilon_all(idsDP),score1avgDP_dist1(idsDP),'b*-','LineWidth',3,'MarkerSize',10); hold on
loglog(epsilon_all(idsDP),score1avgDP_dist2(idsDP),'k*-','LineWidth',3,'MarkerSize',10); hold on


set(gca,'FontSize',FS,'FontWeight','bold')
xlabel('\epsilon','FontSize',FS,'FontWeight','bold');
ylabel('s_1','FontSize',FS,'FontWeight','bold');
title('Synth (N = 800k, \delta = 0.01)','FontSize',FS,'FontWeight','bold')
legend('Pooled CCA','Pooled dpCCA','Dist CCA','Dist dpCCA1','Dist dpCCA2','Location','best')

figure(1) % for s2 vs epsilon
subplot(235)
ids = 1:7;
idsDP = 1:7;

score2CCA_pool = repmat(score2CCA_pool,1,length(epsilon_all));
score2avg_pool = mean(score2CCA_pool(:,ids));

score2CCA_dist = repmat(score2CCA_dist,1,length(epsilon_all));
score2avg_dist = mean(score2CCA_dist(:,ids));

score2avgDP_pool = mean(score2DPCCAAG_pool(:,ids));
score2avgDP_dist1 = mean(score2DPCCAAG_dist1(:,ids));
score2avgDP_dist2 = mean(score2DPCCAAG_dist2(:,ids));

loglog(epsilon_all(idsDP),score2avg_pool(idsDP),'ro-','LineWidth',3,'MarkerSize',10); hold on
loglog(epsilon_all(idsDP),score2avgDP_pool(idsDP),'bo-','LineWidth',3,'MarkerSize',10); hold on

loglog(epsilon_all(idsDP),score2avg_dist(idsDP),'r*-','LineWidth',3,'MarkerSize',10); hold on
loglog(epsilon_all(idsDP),score2avgDP_dist1(idsDP),'b*-','LineWidth',3,'MarkerSize',10); hold on
loglog(epsilon_all(idsDP),score2avgDP_dist2(idsDP),'k*-','LineWidth',3,'MarkerSize',10); hold on


set(gca,'FontSize',FS,'FontWeight','bold')
xlabel('\epsilon','FontSize',FS,'FontWeight','bold');
ylabel('s_2','FontSize',FS,'FontWeight','bold');
title('Synth (N = 800k, \delta = 0.01)','FontSize',FS,'FontWeight','bold')
legend('Pooled CCA','Pooled dpCCA','Dist CCA','Dist dpCCA1','Dist dpCCA2','Location','best')

% for CH index
ids = 1:7;
idsDP = 1:7;
CHindexCCA_pool = repmat(CHindexCCA_pool,1,length(epsilon_all));
CHindexavg_pool = mean(CHindexCCA_pool(:,ids));

CHindexCCA_dist = repmat(CHindexCCA_dist,1,length(epsilon_all));
CHindexavg_dist = mean(CHindexCCA_dist(:,ids));

CHindexavgDP_pool = mean(CHindexDPCCAAG_pool(:,ids));
CHindexavgDP_dist1 = mean(CHindexDPCCAAG_dist1(:,ids));
CHindexavgDP_dist2 = mean(CHindexDPCCAAG_dist2(:,ids));

figure(1)
subplot(236);
loglog(epsilon_all(idsDP),CHindexavg_pool(idsDP),'ro-','LineWidth',3,'MarkerSize',10); hold on
loglog(epsilon_all(idsDP),CHindexavgDP_pool(idsDP),'bo-','LineWidth',3,'MarkerSize',10); hold on

loglog(epsilon_all(ids),CHindexavg_dist(ids),'r*-','LineWidth',3,'MarkerSize',10); hold on
loglog(epsilon_all(idsDP),CHindexavgDP_dist1(idsDP),'b*-','LineWidth',3,'MarkerSize',10); hold on
loglog(epsilon_all(idsDP),CHindexavgDP_dist2(idsDP),'k*-','LineWidth',3,'MarkerSize',10); hold on

% axis([epsilon_all(idsDP(1))/10 epsilon_all(end)*10 1e4 1.8e4])
set(gca,'FontSize',FS,'FontWeight','bold')
xlabel('\epsilon','FontSize',FS,'FontWeight','bold');
ylabel('CH Index','FontSize',FS,'FontWeight','bold');
title('Synth (N = 800k, \delta = 0.01)','FontSize',FS,'FontWeight','bold')
legend('Pooled CCA','Pooled dpCCA','Dist CCA','Dist dpCCA1','Dist dpCCA2','Location','best')