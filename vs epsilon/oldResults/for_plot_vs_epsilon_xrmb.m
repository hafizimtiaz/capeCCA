clear all;clc;close all
ids = [1:7];

FS = 24;


%% for p = 30
load results_vs_eps_xrmb_d50_K20_p30.mat

% for avgErr
ids = [1:7];
perc_errCCA_pool = repmat(perc_errCCA_pool,1,length(epsilon_all));
perc_erravg_pool = mean(perc_errCCA_pool(:,ids));

perc_errCCA_dist = repmat(perc_errCCA_dist,1,length(epsilon_all));
perc_erravg_dist = mean(perc_errCCA_dist(:,ids));

perc_erravgDP_pool = mean(perc_errDPCCAAG_pool(:,ids));
perc_erravgDP_dist1 = mean(perc_errDPCCAAG_dist1(:,ids));
% perc_erravgDP_dist2 = mean(perc_errDPCCAAG_dist2(:,ids));

figure(1)
subplot(231);
ids = [1:7];
idsDP = 3:7;
loglog(epsilon_all(idsDP),perc_erravg_pool(idsDP),'ro-','LineWidth',3,'MarkerSize',10); hold on
loglog(epsilon_all(idsDP),perc_erravgDP_pool(idsDP),'bo-','LineWidth',3,'MarkerSize',10); hold on

% loglog(epsilon_all(idsDP),perc_erravg_dist(idsDP),'r*-','LineWidth',3,'MarkerSize',10); hold on
loglog(epsilon_all(idsDP),perc_erravgDP_dist1(idsDP),'b*-','LineWidth',3,'MarkerSize',10); hold on
% loglog(epsilon_all(idsDP),perc_erravgDP_dist2(idsDP),'k*-','LineWidth',3,'MarkerSize',10); hold on

axis([epsilon_all(idsDP(1))/10 epsilon_all(end)*10 0.3 0.5])
set(gca,'FontSize',FS,'FontWeight','bold')
xlabel('\epsilon','FontSize',FS,'FontWeight','bold');
ylabel('avgErr','FontSize',FS,'FontWeight','bold');
title('XRMB (p = 30, \delta = 0.01)','FontSize',FS,'FontWeight','bold')
% legend('Pooled CCA','Pooled dpCCA','Dist CCA','Dist dpCCA1','Dist dpCCA2','Location','best')
legend('Pooled CCA','Pooled dpCCA','Dist dpCCA','Location','best')

% for CH index
ids = [1:7];
CHindexCCA_pool = repmat(CHindexCCA_pool,1,length(epsilon_all));
CHindexavg_pool = mean(CHindexCCA_pool(:,ids));

CHindexCCA_dist = repmat(CHindexCCA_dist,1,length(epsilon_all));
CHindexavg_dist = mean(CHindexCCA_dist(:,ids));

CHindexavgDP_pool = mean(CHindexDPCCAAG_pool(:,ids));
CHindexavgDP_dist1 = mean(CHindexDPCCAAG_dist1(:,ids));
% CHindexavgDP_dist2 = mean(CHindexDPCCAAG_dist2(:,ids));

figure(1)
subplot(232);
ids = [1:7];
idsDP = 3:7;
loglog(epsilon_all(idsDP),CHindexavg_pool(idsDP),'ro-','LineWidth',3,'MarkerSize',10); hold on
loglog(epsilon_all(idsDP),CHindexavgDP_pool(idsDP),'bo-','LineWidth',3,'MarkerSize',10); hold on

% loglog(epsilon_all(idsDP),CHindexavg_dist(idsDP),'r*-','LineWidth',3,'MarkerSize',10); hold on
loglog(epsilon_all(idsDP),CHindexavgDP_dist1(idsDP),'b*-','LineWidth',3,'MarkerSize',10); hold on
% loglog(epsilon_all(idsDP),CHindexavgDP_dist2(idsDP),'k*-','LineWidth',3,'MarkerSize',10); hold on

axis([epsilon_all(idsDP(1))/10 epsilon_all(end)*10 5e3 1e5])
set(gca,'FontSize',FS,'FontWeight','bold')
xlabel('\epsilon','FontSize',FS,'FontWeight','bold');
ylabel('CH Index','FontSize',FS,'FontWeight','bold');
title('XRMB (p = 30, \delta = 0.01)','FontSize',FS,'FontWeight','bold')
% legend('Pooled CCA','Pooled dpCCA','Dist CCA','Dist dpCCA1','Dist dpCCA2','Location','best')
legend('Pooled CCA','Pooled dpCCA','Dist dpCCA','Location','best')

% for NMI
ids = [1:7];
idsDP = 3:7;
mutInfoCCA_pool = repmat(mutInfoCCA_pool,1,length(epsilon_all));
mutInfoavg_pool = mean(mutInfoCCA_pool(:,ids));

mutInfoCCA_dist = repmat(mutInfoCCA_dist,1,length(epsilon_all));
mutInfoavg_dist = mean(mutInfoCCA_dist(:,ids));

mutInfoavgDP_pool = mean(mutInfoDPCCAAG_pool(:,ids));
mutInfoavgDP_dist1 = mean(mutInfoDPCCAAG_dist1(:,ids));
% mutInfoavgDP_dist2 = mean(mutInfoDPCCAAG_dist2(:,ids));

figure(1)
subplot(233);
ids = [1:7];
loglog(epsilon_all(idsDP),mutInfoavg_pool(idsDP),'ro-','LineWidth',3,'MarkerSize',10); hold on
loglog(epsilon_all(idsDP),mutInfoavgDP_pool(idsDP),'bo-','LineWidth',3,'MarkerSize',10); hold on

% loglog(epsilon_all(idsDP),mutInfoavg_dist(idsDP),'r*-','LineWidth',3,'MarkerSize',10); hold on
loglog(epsilon_all(idsDP),mutInfoavgDP_dist1(idsDP),'b*-','LineWidth',3,'MarkerSize',10); hold on
% loglog(epsilon_all(idsDP),mutInfoavgDP_dist2(idsDP),'k*-','LineWidth',3,'MarkerSize',10); hold on

axis([epsilon_all(idsDP(1))/10 epsilon_all(end)*10 1e-3 1e-1])
set(gca,'FontSize',FS,'FontWeight','bold')
xlabel('\epsilon','FontSize',FS,'FontWeight','bold');
ylabel('NMI','FontSize',FS,'FontWeight','bold');
title('XRMB (p = 30, \delta = 0.01)','FontSize',FS,'FontWeight','bold')
% legend('Pooled CCA','Pooled dpCCA','Dist CCA','Dist dpCCA1','Dist dpCCA2','Location','best')
legend('Pooled CCA','Pooled dpCCA','Dist dpCCA','Location','best')


%% for p = 50
load results_vs_eps_xrmb_d50_K20_p50.mat

% for avgErr
ids = [1:7];
idsDP = 3:7;
perc_errCCA_pool = repmat(perc_errCCA_pool,1,length(epsilon_all));
perc_erravg_pool = mean(perc_errCCA_pool(:,ids));

perc_errCCA_dist = repmat(perc_errCCA_dist,1,length(epsilon_all));
perc_erravg_dist = mean(perc_errCCA_dist(:,ids));

perc_erravgDP_pool = mean(perc_errDPCCAAG_pool(:,ids));
perc_erravgDP_dist1 = mean(perc_errDPCCAAG_dist1(:,ids));
% perc_erravgDP_dist2 = mean(perc_errDPCCAAG_dist2(:,ids));

figure(1)
subplot(234);
ids = [1:7];
loglog(epsilon_all(idsDP),perc_erravg_pool(idsDP),'ro-','LineWidth',3,'MarkerSize',10); hold on
loglog(epsilon_all(idsDP),perc_erravgDP_pool(idsDP),'bo-','LineWidth',3,'MarkerSize',10); hold on

% loglog(epsilon_all(idsDP),perc_erravg_dist(idsDP),'r*-','LineWidth',3,'MarkerSize',10); hold on
loglog(epsilon_all(idsDP),perc_erravgDP_dist1(idsDP),'b*-','LineWidth',3,'MarkerSize',10); hold on
% loglog(epsilon_all(idsDP),perc_erravgDP_dist2(idsDP),'k*-','LineWidth',3,'MarkerSize',10); hold on

axis([epsilon_all(idsDP(1))/10 epsilon_all(end)*10 0.3 0.5])
set(gca,'FontSize',FS,'FontWeight','bold')
xlabel('\epsilon','FontSize',FS,'FontWeight','bold');
ylabel('avgErr','FontSize',FS,'FontWeight','bold');
title('XRMB (p = 50, \delta = 0.01)','FontSize',FS,'FontWeight','bold')
% legend('Pooled CCA','Pooled dpCCA','Dist CCA','Dist dpCCA1','Dist dpCCA2','Location','best')
legend('Pooled CCA','Pooled dpCCA','Dist dpCCA','Location','best')

% for CH index
ids = [1:7];
idsDP = 3:7;
CHindexCCA_pool = repmat(CHindexCCA_pool,1,length(epsilon_all));
CHindexavg_pool = mean(CHindexCCA_pool(:,ids));

CHindexCCA_dist = repmat(CHindexCCA_dist,1,length(epsilon_all));
CHindexavg_dist = mean(CHindexCCA_dist(:,ids));

CHindexavgDP_pool = mean(CHindexDPCCAAG_pool(:,ids));
CHindexavgDP_dist1 = mean(CHindexDPCCAAG_dist1(:,ids));
% CHindexavgDP_dist2 = mean(CHindexDPCCAAG_dist2(:,ids));

figure(1)
subplot(235);
ids = [1:7];
loglog(epsilon_all(idsDP),CHindexavg_pool(idsDP),'ro-','LineWidth',3,'MarkerSize',10); hold on
loglog(epsilon_all(idsDP),CHindexavgDP_pool(idsDP),'bo-','LineWidth',3,'MarkerSize',10); hold on

% loglog(epsilon_all(idsDP),CHindexavg_dist(idsDP),'r*-','LineWidth',3,'MarkerSize',10); hold on
loglog(epsilon_all(idsDP),CHindexavgDP_dist1(idsDP),'b*-','LineWidth',3,'MarkerSize',10); hold on
% loglog(epsilon_all(idsDP),CHindexavgDP_dist2(idsDP),'k*-','LineWidth',3,'MarkerSize',10); hold on

axis([epsilon_all(idsDP(1))/10 epsilon_all(end)*10 1e4 5e4])
set(gca,'FontSize',FS,'FontWeight','bold')
xlabel('\epsilon','FontSize',FS,'FontWeight','bold');
ylabel('CH Index','FontSize',FS,'FontWeight','bold');
title('XRMB (p = 50, \delta = 0.01)','FontSize',FS,'FontWeight','bold')
% legend('Pooled CCA','Pooled dpCCA','Dist CCA','Dist dpCCA1','Dist dpCCA2','Location','best')
legend('Pooled CCA','Pooled dpCCA','Dist dpCCA','Location','best')

% for NMI
ids = [1:7];
idsDP = 3:7;
mutInfoCCA_pool = repmat(mutInfoCCA_pool,1,length(epsilon_all));
mutInfoavg_pool = mean(mutInfoCCA_pool(:,ids));

mutInfoCCA_dist = repmat(mutInfoCCA_dist,1,length(epsilon_all));
mutInfoavg_dist = mean(mutInfoCCA_dist(:,ids));

mutInfoavgDP_pool = mean(mutInfoDPCCAAG_pool(:,ids));
mutInfoavgDP_dist1 = mean(mutInfoDPCCAAG_dist1(:,ids));
% mutInfoavgDP_dist2 = mean(mutInfoDPCCAAG_dist2(:,ids));

figure(1)
subplot(236);
ids = [1:7];
loglog(epsilon_all(idsDP),mutInfoavg_pool(idsDP),'ro-','LineWidth',3,'MarkerSize',10); hold on
loglog(epsilon_all(idsDP),mutInfoavgDP_pool(idsDP),'bo-','LineWidth',3,'MarkerSize',10); hold on

% loglog(epsilon_all(idsDP),mutInfoavg_dist(idsDP),'r*-','LineWidth',3,'MarkerSize',10); hold on
loglog(epsilon_all(idsDP),mutInfoavgDP_dist1(idsDP),'b*-','LineWidth',3,'MarkerSize',10); hold on
% loglog(epsilon_all(idsDP),mutInfoavgDP_dist2(idsDP),'k*-','LineWidth',3,'MarkerSize',10); hold on

axis([epsilon_all(idsDP(1))/10 epsilon_all(end)*10 1e-3 1e-1])
set(gca,'FontSize',FS,'FontWeight','bold')
xlabel('\epsilon','FontSize',FS,'FontWeight','bold');
ylabel('NMI','FontSize',FS,'FontWeight','bold');
title('XRMB (p = 50, \delta = 0.01)','FontSize',FS,'FontWeight','bold')
% legend('Pooled CCA','Pooled dpCCA','Dist CCA','Dist dpCCA1','Dist dpCCA2','Location','best')
legend('Pooled CCA','Pooled dpCCA','Dist dpCCA','Location','best')