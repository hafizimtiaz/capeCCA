clear all;clc;close all
ids = [1:7];

FS = 24;


%% for N = 30k
load results_vs_eps_mnist_d100_K50_N30k.mat

% for avgErr
ids = [1:7];
perc_errCCA_pool = repmat(perc_errCCA_pool,1,length(epsilon_all));
perc_erravg_pool = mean(perc_errCCA_pool(:,ids));

perc_errCCA_dist = repmat(perc_errCCA_dist,1,length(epsilon_all));
% perc_erravg_dist = mean(perc_errCCA_dist(:,ids));

perc_erravgDP_pool = mean(perc_errDPCCAAG_pool(:,ids));
perc_erravgDP_dist1 = mean(perc_errDPCCAAG_dist1(:,ids));
% perc_erravgDP_dist2 = mean(perc_errDPCCAAG_dist2(:,ids));

figure(1)
subplot(231);
ids = [1:7];
idsDP = 2:7;
loglog(epsilon_all(idsDP),perc_erravg_pool(idsDP),'ro-','LineWidth',3,'MarkerSize',10); hold on
loglog(epsilon_all(idsDP),perc_erravgDP_pool(idsDP),'bo-','LineWidth',3,'MarkerSize',10); hold on

% loglog(epsilon_all(ids),perc_erravg_dist(ids),'r*-','LineWidth',3,'MarkerSize',10); hold on
loglog(epsilon_all(3:end),perc_erravgDP_dist1(3:end),'b*-','LineWidth',3,'MarkerSize',10); hold on
% loglog(epsilon_all(idsDP),perc_erravgDP_dist2(idsDP),'k*-','LineWidth',3,'MarkerSize',10); hold on

axis([epsilon_all(idsDP(1))/10 epsilon_all(end)*10 0.5 0.7])
set(gca,'FontSize',FS,'FontWeight','bold')
xlabel('\epsilon','FontSize',FS,'FontWeight','bold');
ylabel('avgErr','FontSize',FS,'FontWeight','bold');
title('MNIST (N = 30k, \delta = 0.01)','FontSize',FS,'FontWeight','bold')
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

axis([epsilon_all(idsDP(1))/10 epsilon_all(end)*10 1e4 1.8e4])
set(gca,'FontSize',FS,'FontWeight','bold')
xlabel('\epsilon','FontSize',FS,'FontWeight','bold');
ylabel('CH Index','FontSize',FS,'FontWeight','bold');
title('MNIST (N = 30k, \delta = 0.01)','FontSize',FS,'FontWeight','bold')
legend('Pooled CCA','Pooled dpCCA','Dist dpCCA','Location','best')

% for NMI
ids = [1:7];
idsDP = 2:7;
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
loglog(epsilon_all(3:end),mutInfoavgDP_dist1(3:end),'b*-','LineWidth',3,'MarkerSize',10); hold on
% loglog(epsilon_all(idsDP),mutInfoavgDP_dist2(idsDP),'k*-','LineWidth',3,'MarkerSize',10); hold on

axis([epsilon_all(idsDP(1))/10 epsilon_all(end)*10 0.2 0.4])
set(gca,'FontSize',FS,'FontWeight','bold')
xlabel('\epsilon','FontSize',FS,'FontWeight','bold');
ylabel('NMI','FontSize',FS,'FontWeight','bold');
title('MNIST (N = 30k, \delta = 0.01)','FontSize',FS,'FontWeight','bold')
legend('Pooled CCA','Pooled dpCCA','Dist dpCCA','Location','best')


%% for N = 50k
load results_vs_eps_mnist_d100_K50_N50k.mat

% for avgErr
ids = [1:7];
idsDP = 2:7;
perc_errCCA_pool = repmat(perc_errCCA_pool,1,length(epsilon_all));
perc_erravg_pool = mean(perc_errCCA_pool(:,ids));

perc_errCCA_dist = repmat(perc_errCCA_dist,1,length(epsilon_all));
% perc_erravg_dist = mean(perc_errCCA_dist(:,ids));

perc_erravgDP_pool = mean(perc_errDPCCAAG_pool(:,ids));
perc_erravgDP_dist1 = mean(perc_errDPCCAAG_dist1(:,ids));
% perc_erravgDP_dist2 = mean(perc_errDPCCAAG_dist2(:,ids));

figure(1)
subplot(234);
loglog(epsilon_all(idsDP),perc_erravg_pool(idsDP),'ro-','LineWidth',3,'MarkerSize',10); hold on
loglog(epsilon_all(idsDP),perc_erravgDP_pool(idsDP),'bo-','LineWidth',3,'MarkerSize',10); hold on

% loglog(epsilon_all(ids),perc_erravg_dist(ids),'r*-','LineWidth',3,'MarkerSize',10); hold on
loglog(epsilon_all(idsDP),perc_erravgDP_dist1(idsDP),'b*-','LineWidth',3,'MarkerSize',10); hold on
% loglog(epsilon_all(idsDP),perc_erravgDP_dist2(idsDP),'k*-','LineWidth',3,'MarkerSize',10); hold on

axis([epsilon_all(idsDP(1))/10 epsilon_all(end)*10 0.5 0.7])
set(gca,'FontSize',FS,'FontWeight','bold')
xlabel('\epsilon','FontSize',FS,'FontWeight','bold');
ylabel('avgErr','FontSize',FS,'FontWeight','bold');
title('MNIST (N = 50k, \delta = 0.01)','FontSize',FS,'FontWeight','bold')
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

axis([epsilon_all(idsDP(1))/10 epsilon_all(end)*10 1.7e4 3e4])
set(gca,'FontSize',FS,'FontWeight','bold')
xlabel('\epsilon','FontSize',FS,'FontWeight','bold');
ylabel('CH Index','FontSize',FS,'FontWeight','bold');
title('MNIST (N = 50k, \delta = 0.01)','FontSize',FS,'FontWeight','bold')
legend('Pooled CCA','Pooled dpCCA','Dist dpCCA','Location','best')

% for NMI
ids = [1:7];
idsDP = 2:7;
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
loglog(epsilon_all(idsDP),mutInfoavgDP_dist1(idsDP),'b*-','LineWidth',3,'MarkerSize',10); hold on
% loglog(epsilon_all(idsDP),mutInfoavgDP_dist2(idsDP),'k*-','LineWidth',3,'MarkerSize',10); hold on

axis([epsilon_all(idsDP(1))/10 epsilon_all(end)*10 0.2 0.4])
set(gca,'FontSize',FS,'FontWeight','bold')
xlabel('\epsilon','FontSize',FS,'FontWeight','bold');
ylabel('NMI','FontSize',FS,'FontWeight','bold');
title('MNIST (N = 50k, \delta = 0.01)','FontSize',FS,'FontWeight','bold')
legend('Pooled CCA','Pooled dpCCA','Dist dpCCA','Location','best')