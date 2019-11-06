clear all;clc;close all
ids = [1:7];

FS = 24;

%% for N = 200k samples
load results_vs_eps_synth_d100_K20_N200k_v2.mat

% for s2
figure(1)
subplot(231);
ids = 1:7;
idsDP = 2:7;

score2CCA_pool = repmat(score2CCA_pool,1,length(epsilon_all));
score2avg_pool = mean(score2CCA_pool(:,ids));

score2CCA_dist = repmat(score2CCA_dist,1,length(epsilon_all));
% score2avg_dist = mean(score2CCA_dist(:,ids));

score2avgDP_pool = mean(score2DPCCAAG_pool(:,ids));
score2avgDP_dist1 = mean(score2DPCCAAG_dist1(:,ids));
% score2avgDP_dist2 = mean(score2DPCCAAG_dist2(:,ids));

loglog(epsilon_all(idsDP),score2avg_pool(idsDP),'ro-','LineWidth',3,'MarkerSize',10); hold on
loglog(epsilon_all(idsDP),score2avgDP_pool(idsDP),'bo-','LineWidth',3,'MarkerSize',10); hold on

% loglog(epsilon_all(idsDP),score2avg_dist(idsDP),'r*-','LineWidth',3,'MarkerSize',10); hold on
loglog(epsilon_all(idsDP),score2avgDP_dist1(idsDP),'b*-','LineWidth',3,'MarkerSize',10); hold on
% loglog(epsilon_all(idsDP),score2avgDP_dist2(idsDP),'k*-','LineWidth',3,'MarkerSize',10); hold on

axis([epsilon_all(idsDP(1))/10 epsilon_all(end)*10 1e-4 1])
set(gca,'FontSize',FS,'FontWeight','bold')
xlabel('\epsilon','FontSize',FS,'FontWeight','bold');
ylabel('s_2','FontSize',FS,'FontWeight','bold');
title('Synth (N = 200k, \delta = 0.01)','FontSize',FS,'FontWeight','bold')
legend('Pooled CCA','Pooled dpCCA','Dist dpCCA','Location','best')

% for CH index
ids = [1:7];
CHindexCCA_pool = repmat(CHindexCCA_pool,1,length(epsilon_all));
CHindexavg_pool = mean(CHindexCCA_pool(:,ids));

CHindexCCA_dist = repmat(CHindexCCA_dist,1,length(epsilon_all));
% CHindexavg_dist = mean(CHindexCCA_dist(:,ids));

CHindexavgDP_pool = mean(CHindexDPCCAAG_pool(:,ids));
CHindexavgDP_dist1 = mean(CHindexDPCCAAG_dist1(:,ids));
% CHindexavgDP_dist2 = mean(CHindexDPCCAAG_dist2(:,ids));

figure(1)
subplot(232);
idsDP = 2:7;
loglog(epsilon_all(idsDP),CHindexavg_pool(idsDP),'ro-','LineWidth',3,'MarkerSize',10); hold on
loglog(epsilon_all(idsDP),CHindexavgDP_pool(idsDP),'bo-','LineWidth',3,'MarkerSize',10); hold on

% loglog(epsilon_all(ids),CHindexavg_dist(ids),'r*-','LineWidth',3,'MarkerSize',10); hold on
loglog(epsilon_all(idsDP),CHindexavgDP_dist1(idsDP),'b*-','LineWidth',3,'MarkerSize',10); hold on
% loglog(epsilon_all(idsDP),CHindexavgDP_dist2(idsDP),'k*-','LineWidth',3,'MarkerSize',10); hold on

axis([epsilon_all(idsDP(1))/10 epsilon_all(end)*10 1e4 1e6])
set(gca,'FontSize',FS,'FontWeight','bold')
xlabel('\epsilon','FontSize',FS,'FontWeight','bold');
ylabel('CH Index','FontSize',FS,'FontWeight','bold');
title('MNIST (N = 30k, \delta = 0.01)','FontSize',FS,'FontWeight','bold')
legend('Pooled CCA','Pooled dpCCA','Dist dpCCA','Location','best')

% for NMI
ids = [1:7];
idsDP = 3:7;
mutInfoCCA_pool = repmat(mutInfoCCA_pool,1,length(epsilon_all));
mutInfoavg_pool = mean(mutInfoCCA_pool(:,ids));

mutInfoCCA_dist = repmat(mutInfoCCA_dist,1,length(epsilon_all));
% mutInfoavg_dist = mean(mutInfoCCA_dist(:,ids));

mutInfoavgDP_pool = mean(mutInfoDPCCAAG_pool(:,ids));
mutInfoavgDP_dist1 = mean(mutInfoDPCCAAG_dist1(:,ids));
% mutInfoavgDP_dist2 = mean(mutInfoDPCCAAG_dist2(:,ids));

figure(1)
subplot(233);
loglog(epsilon_all(idsDP),mutInfoavg_pool(idsDP),'ro-','LineWidth',3,'MarkerSize',10); hold on
loglog(epsilon_all(idsDP),mutInfoavgDP_pool(idsDP),'bo-','LineWidth',3,'MarkerSize',10); hold on

% loglog(epsilon_all(ids),mutInfoavg_dist(ids),'r*-','LineWidth',3,'MarkerSize',10); hold on
loglog(epsilon_all(idsDP),mutInfoavgDP_dist1(idsDP),'b*-','LineWidth',3,'MarkerSize',10); hold on
% loglog(epsilon_all(idsDP),mutInfoavgDP_dist2(idsDP),'k*-','LineWidth',3,'MarkerSize',10); hold on

axis([epsilon_all(idsDP(1))/10 epsilon_all(end)*10 0.45 0.7])
set(gca,'FontSize',FS,'FontWeight','bold')
xlabel('\epsilon','FontSize',FS,'FontWeight','bold');
ylabel('NMI','FontSize',FS,'FontWeight','bold');
title('MNIST (N = 30k, \delta = 0.01)','FontSize',FS,'FontWeight','bold')
legend('Pooled CCA','Pooled dpCCA','Dist dpCCA','Location','best')


%% for N = 800k
load results_vs_eps_synth_d100_K20_N800k_v2.mat

% for s2
figure(1)
subplot(234);
ids = 1:7;
idsDP = 2:7;

score2CCA_pool = repmat(score2CCA_pool,1,length(epsilon_all));
score2avg_pool = mean(score2CCA_pool(:,ids));

score2CCA_dist = repmat(score2CCA_dist,1,length(epsilon_all));
% score2avg_dist = mean(score2CCA_dist(:,ids));

score2avgDP_pool = mean(score2DPCCAAG_pool(:,ids));
score2avgDP_dist1 = mean(score2DPCCAAG_dist1(:,ids));
% score2avgDP_dist2 = mean(score2DPCCAAG_dist2(:,ids));

loglog(epsilon_all(idsDP),score2avg_pool(idsDP),'ro-','LineWidth',3,'MarkerSize',10); hold on
loglog(epsilon_all(idsDP),score2avgDP_pool(idsDP),'bo-','LineWidth',3,'MarkerSize',10); hold on

% loglog(epsilon_all(idsDP),score2avg_dist(idsDP),'r*-','LineWidth',3,'MarkerSize',10); hold on
loglog(epsilon_all(idsDP),score2avgDP_dist1(idsDP),'b*-','LineWidth',3,'MarkerSize',10); hold on
% loglog(epsilon_all(idsDP),score2avgDP_dist2(idsDP),'k*-','LineWidth',3,'MarkerSize',10); hold on

axis([epsilon_all(idsDP(1))/10 epsilon_all(end)*10 5e-5 1])
set(gca,'FontSize',FS,'FontWeight','bold')
xlabel('\epsilon','FontSize',FS,'FontWeight','bold');
ylabel('s_2','FontSize',FS,'FontWeight','bold');
title('Synth (N = 800k, \delta = 0.01)','FontSize',FS,'FontWeight','bold')
legend('Pooled CCA','Pooled dpCCA','Dist dpCCA','Location','best')

% for CH index
ids = [1:7];
idsDP = 2:7;
CHindexCCA_pool = repmat(CHindexCCA_pool,1,length(epsilon_all));
CHindexavg_pool = mean(CHindexCCA_pool(:,ids));

CHindexCCA_dist = repmat(CHindexCCA_dist,1,length(epsilon_all));
% CHindexavg_dist = mean(CHindexCCA_dist(:,ids));

CHindexavgDP_pool = mean(CHindexDPCCAAG_pool(:,ids));
CHindexavgDP_dist1 = mean(CHindexDPCCAAG_dist1(:,ids));
% CHindexavgDP_dist2 = mean(CHindexDPCCAAG_dist2(:,ids));

figure(1)
subplot(235);
loglog(epsilon_all(idsDP),CHindexavg_pool(idsDP),'ro-','LineWidth',3,'MarkerSize',10); hold on
loglog(epsilon_all(idsDP),CHindexavgDP_pool(idsDP),'bo-','LineWidth',3,'MarkerSize',10); hold on

% loglog(epsilon_all(ids),CHindexavg_dist(ids),'r*-','LineWidth',3,'MarkerSize',10); hold on
loglog(epsilon_all(idsDP),CHindexavgDP_dist1(idsDP),'b*-','LineWidth',3,'MarkerSize',10); hold on
% loglog(epsilon_all(idsDP),CHindexavgDP_dist2(idsDP),'k*-','LineWidth',3,'MarkerSize',10); hold on

axis([epsilon_all(idsDP(1))/10 epsilon_all(end)*10 8e4 2e6])
set(gca,'FontSize',FS,'FontWeight','bold')
xlabel('\epsilon','FontSize',FS,'FontWeight','bold');
ylabel('CH Index','FontSize',FS,'FontWeight','bold');
title('MNIST (N = 50k, \delta = 0.01)','FontSize',FS,'FontWeight','bold')
legend('Pooled CCA','Pooled dpCCA','Dist dpCCA','Location','best')

% for NMI
ids = [1:7];
idsDP = 3:7;
mutInfoCCA_pool = repmat(mutInfoCCA_pool,1,length(epsilon_all));
mutInfoavg_pool = mean(mutInfoCCA_pool(:,ids));

mutInfoCCA_dist = repmat(mutInfoCCA_dist,1,length(epsilon_all));
% mutInfoavg_dist = mean(mutInfoCCA_dist(:,ids));

mutInfoavgDP_pool = mean(mutInfoDPCCAAG_pool(:,ids));
mutInfoavgDP_dist1 = mean(mutInfoDPCCAAG_dist1(:,ids));
% mutInfoavgDP_dist2 = mean(mutInfoDPCCAAG_dist2(:,ids));

figure(1)
subplot(236);
loglog(epsilon_all(idsDP),mutInfoavg_pool(idsDP),'ro-','LineWidth',3,'MarkerSize',10); hold on
loglog(epsilon_all(idsDP),mutInfoavgDP_pool(idsDP),'bo-','LineWidth',3,'MarkerSize',10); hold on

% loglog(epsilon_all(ids),mutInfoavg_dist(ids),'r*-','LineWidth',3,'MarkerSize',10); hold on
loglog(epsilon_all(4:end),mutInfoavgDP_dist1(4:end),'b*-','LineWidth',3,'MarkerSize',10); hold on
% loglog(epsilon_all(idsDP),mutInfoavgDP_dist2(idsDP),'k*-','LineWidth',3,'MarkerSize',10); hold on

axis([epsilon_all(idsDP(1))/10 epsilon_all(end)*10 0.45 0.7])
set(gca,'FontSize',FS,'FontWeight','bold')
xlabel('\epsilon','FontSize',FS,'FontWeight','bold');
ylabel('NMI','FontSize',FS,'FontWeight','bold');
title('MNIST (N = 50k, \delta = 0.01)','FontSize',FS,'FontWeight','bold')
legend('Pooled CCA','Pooled dpCCA','Dist dpCCA','Location','best')