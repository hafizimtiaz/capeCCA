clear all;clc;close all
ids = [1:7];

FS = 24;

%% SYNTH DATA
% for N = 200k samples
load results_vs_eps_synth_d100_K20_N200k_v2.mat

% for s2
figure(1)
subplot(361);
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
subplot(363);
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
title('Synth (N = 200k, \delta = 0.01)','FontSize',FS,'FontWeight','bold')
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
subplot(365);
loglog(epsilon_all(idsDP),mutInfoavg_pool(idsDP),'ro-','LineWidth',3,'MarkerSize',10); hold on
loglog(epsilon_all(idsDP),mutInfoavgDP_pool(idsDP),'bo-','LineWidth',3,'MarkerSize',10); hold on

% loglog(epsilon_all(ids),mutInfoavg_dist(ids),'r*-','LineWidth',3,'MarkerSize',10); hold on
loglog(epsilon_all(idsDP),mutInfoavgDP_dist1(idsDP),'b*-','LineWidth',3,'MarkerSize',10); hold on
% loglog(epsilon_all(idsDP),mutInfoavgDP_dist2(idsDP),'k*-','LineWidth',3,'MarkerSize',10); hold on

axis([epsilon_all(idsDP(1))/10 epsilon_all(end)*10 0.45 0.7])
set(gca,'FontSize',FS,'FontWeight','bold')
xlabel('\epsilon','FontSize',FS,'FontWeight','bold');
ylabel('NMI','FontSize',FS,'FontWeight','bold');
title('Synth (N = 200k, \delta = 0.01)','FontSize',FS,'FontWeight','bold')
legend('Pooled CCA','Pooled dpCCA','Dist dpCCA','Location','best')


% for N = 800k
load results_vs_eps_synth_d100_K20_N800k_v2.mat

% for s2
figure(1)
subplot(362);
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
subplot(364);
loglog(epsilon_all(idsDP),CHindexavg_pool(idsDP),'ro-','LineWidth',3,'MarkerSize',10); hold on
loglog(epsilon_all(idsDP),CHindexavgDP_pool(idsDP),'bo-','LineWidth',3,'MarkerSize',10); hold on

% loglog(epsilon_all(ids),CHindexavg_dist(ids),'r*-','LineWidth',3,'MarkerSize',10); hold on
loglog(epsilon_all(idsDP),CHindexavgDP_dist1(idsDP),'b*-','LineWidth',3,'MarkerSize',10); hold on
% loglog(epsilon_all(idsDP),CHindexavgDP_dist2(idsDP),'k*-','LineWidth',3,'MarkerSize',10); hold on

axis([epsilon_all(idsDP(1))/10 epsilon_all(end)*10 8e4 2e6])
set(gca,'FontSize',FS,'FontWeight','bold')
xlabel('\epsilon','FontSize',FS,'FontWeight','bold');
ylabel('CH Index','FontSize',FS,'FontWeight','bold');
title('Synth (N = 800k, \delta = 0.01)','FontSize',FS,'FontWeight','bold')
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
subplot(366);
loglog(epsilon_all(idsDP),mutInfoavg_pool(idsDP),'ro-','LineWidth',3,'MarkerSize',10); hold on
loglog(epsilon_all(idsDP),mutInfoavgDP_pool(idsDP),'bo-','LineWidth',3,'MarkerSize',10); hold on

% loglog(epsilon_all(ids),mutInfoavg_dist(ids),'r*-','LineWidth',3,'MarkerSize',10); hold on
loglog(epsilon_all(4:end),mutInfoavgDP_dist1(4:end),'b*-','LineWidth',3,'MarkerSize',10); hold on
% loglog(epsilon_all(idsDP),mutInfoavgDP_dist2(idsDP),'k*-','LineWidth',3,'MarkerSize',10); hold on

axis([epsilon_all(idsDP(1))/10 epsilon_all(end)*10 0.45 0.7])
set(gca,'FontSize',FS,'FontWeight','bold')
xlabel('\epsilon','FontSize',FS,'FontWeight','bold');
ylabel('NMI','FontSize',FS,'FontWeight','bold');
title('Synth (N = 800k, \delta = 0.01)','FontSize',FS,'FontWeight','bold')
legend('Pooled CCA','Pooled dpCCA','Dist dpCCA','Location','best')

%% XRMB DATA
% for p = 30
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
subplot(367);
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
subplot(369);
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
subplot(3,6,11);
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


% for p = 50
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
subplot(368);
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
subplot(3,6,10);
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
subplot(3,6,12);
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

%% MNIST DATA
% for N = 30k
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
subplot(3,6,13);
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
subplot(3,6,15);
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
subplot(3,6,17);
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


% for N = 50k
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
subplot(3,6,14);
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
subplot(3,6,16);
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
subplot(3,6,18);
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