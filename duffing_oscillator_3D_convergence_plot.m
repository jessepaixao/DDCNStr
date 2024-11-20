%% Data-driven Computing in Nonlinear Dynamics Systems 
%  Description: Plot convergece analysis for the Duffing oscillator
%  Versions:
%   1.0 - 15/07/24 Jessé Paixão

close all
clear all
clc

% return

% Load library of functions
addpath('lib');

% Figure settings
set(groot,'defaultAxesFontSize',24)
set(groot,'defaultTextInterpreter','latex')
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');

%% DEFINE ACQUISITION PARAMETERS AND EXCITATION SIGNAL

% Load processed data
load("results\03_8_24-14_14_duffing_3D_conv")

Fs=5e3;  % Sampling frequency of generated dataset

N=[length([0:1/Fs:80]) length([0:1/Fs:200]) length([0:1/Fs:600]) length([0:1/Fs:1000])]

% Figure 8a
figure(1)
% Xtick=[2e5 1e6 6e6];
clear hLegend
set(gcf,'units','normalized','outerposition',[0 0 1 0.7])
ax = gca;
% Noise 30 db
SNR_jj=1;
% h=boxplot(RMSE_x(:,:,SNR_jj)','Positions',N,'PlotStyle','compact','Colors','k','symbol', '');hold on
% xticks(Xtick)
% xticklabels(num2cell(Xtick))
% xtickformat('%1.e')
CI=std(RMSE_x(:,:,SNR_jj)')*2;
x_CI =[N, fliplr(N)];
y_CI=[mean(RMSE_x(:,:,SNR_jj),2)'+CI(:)', flipud(mean(RMSE_x(:,:,SNR_jj),2)-CI(:))'];
fill(x_CI, y_CI, 1,'facecolor', 'k', 'edgecolor', 'none', 'facealpha', 0.4,'HandleVisibility','off'); hold on
plot(N,mean(RMSE_x(:,:,SNR_jj),2),'k-x','LineWidth',1.5,'MarkerSize',10); hold on
% box_vars = findall(gca,'Tag','Box');
% hLegend=box_vars;
% Noise 40 db
SNR_jj=3;
% h=boxplot(RMSE_x(:,:,SNR_jj)','Positions',N,'PlotStyle','compact','Colors','r','symbol', '');hold on
% box_vars = findall(gca,'Tag','Box');
% hLegend=box_vars;
% xticks(Xtick)
% xticklabels(num2cell(Xtick))
% xtickformat('%1.e')
CI=std(RMSE_x(:,:,SNR_jj)')*2;
x_CI =[N, fliplr(N)];
y_CI=[mean(RMSE_x(:,:,SNR_jj),2)'+CI(:)', flipud(mean(RMSE_x(:,:,SNR_jj),2)-CI(:))'];
fill(x_CI, y_CI, 1,'facecolor', 'r', 'edgecolor', 'none', 'facealpha', 0.4,'HandleVisibility','off'); hold on
plot(N,mean(RMSE_x(:,:,SNR_jj),2),'-or','LineWidth',1.5,'MarkerSize',10); hold on
% Noise 45 db
SNR_jj=4;
% h=boxplot(RMSE_x(:,:,SNR_jj)','Positions',N,'BoxStyle','filled','Widths',50,'Colors','b','symbol', '');hold on
% box_vars = findall(gca,'Tag','Box');
% hLegend=box_vars;
% xticks(Xtick)
% xticklabels(num2cell(Xtick))
% xtickformat('%1.e')
CI=std(RMSE_x(:,:,SNR_jj)')*2;
x_CI =[N, fliplr(N)];
y_CI=[mean(RMSE_x(:,:,SNR_jj),2)'+CI(:)', flipud(mean(RMSE_x(:,:,SNR_jj),2)-CI(:))'];
fill(x_CI, y_CI, 1,'facecolor', 'b', 'edgecolor', 'none', 'facealpha', 0.4,'HandleVisibility','off'); hold on
plot(N,mean(RMSE_x(:,:,SNR_jj),2),'-sb','LineWidth',1.5,'MarkerSize',10); hold on
% No Noise 
SNR_jj=5;
% h=boxplot(RMSE_x(:,:,SNR_jj)','Positions',N,'PlotStyle','compact','Colors','m','symbol', '');hold on
% box_vars = findall(gca,'Tag','Box');
% hLegend=box_vars;
% xticks(Xtick)
% xticklabels(num2cell(Xtick))
% xtickformat('%1.e')
CI=std(RMSE_x(:,:,SNR_jj)')*2;
x_CI =[N, fliplr(N)];
y_CI=[mean(RMSE_x(:,:,SNR_jj),2)'+CI(:)', flipud(mean(RMSE_x(:,:,SNR_jj),2)-CI(:))'];
fill(x_CI, y_CI, 1,'facecolor',"#D95319" , 'edgecolor', 'none', 'facealpha', 0.4,'HandleVisibility','off'); hold on
plot(N,mean(RMSE_x(:,:,SNR_jj),2),'-*','Color',"#D95319",'LineWidth',1.5,'MarkerSize',10); hold on
ax.XAxis.Scale ="log";
set(gca,'YLim',[0 4e1])
ylabel(' RMSE $x$ ($\%$)')
% xlabel('Dataset Size')
% grid minor
legend(flip({'SNR 30 dB (3.16 \%)','SNR 40 dB (1 \%)','SNR 45 dB (0.56 \%)','No Noise (0 \%)'}),'Position',[0.256177169474693 0.934362056709471 0.534832648075864 0.0636645976060666],...
    'Orientation','horizontal')

% Save figure
exportgraphics(gcf,'figures/duffing_3D_convergence_noise_a.pdf','ContentType','vector')
% box_h = findobj(gca,'Tag','Box');
% for i = 1:length(box_h)
%     hp(i) = patch([box_h(i).XData],[box_h(i).YData],[0.5 0.5 0.5],'FaceAlpha',.2);
% end
% uistack(box_h,'bottom')
%

% Figure 8b
figure(2)
Xtick=[2e5 1e6 6e6];
clear hLegend
set(gcf,'units','normalized','outerposition',[0 0 1 0.7])
ax = gca;
% Noise 30 db
SNR_jj=1;
CI=std(RMSE_dx(:,:,SNR_jj)')*2;
x_CI =[N, fliplr(N)];
y_CI=[mean(RMSE_dx(:,:,SNR_jj),2)'+CI(:)', flipud(mean(RMSE_dx(:,:,SNR_jj),2)-CI(:))'];
fill(x_CI, y_CI, 1,'facecolor', 'k', 'edgecolor', 'none', 'facealpha', 0.4,'HandleVisibility','off'); hold on
plot(N,mean(RMSE_dx(:,:,SNR_jj),2),'k-x','LineWidth',1.5,'MarkerSize',10); hold on
box_vars = findall(gca,'Tag','Box');
% Noise 40 db
SNR_jj=3;
CI=std(RMSE_dx(:,:,SNR_jj)')*2;
x_CI =[N, fliplr(N)];
y_CI=[mean(RMSE_dx(:,:,SNR_jj),2)'+CI(:)', flipud(mean(RMSE_dx(:,:,SNR_jj),2)-CI(:))'];
fill(x_CI, y_CI, 1,'facecolor', 'r', 'edgecolor', 'none', 'facealpha', 0.4,'HandleVisibility','off'); hold on
plot(N,mean(RMSE_dx(:,:,SNR_jj),2),'r-o','LineWidth',1.5,'MarkerSize',10); hold on
box_vars = findall(gca,'Tag','Box');
% Noise 45 db
SNR_jj=4;
CI=std(RMSE_dx(:,:,SNR_jj)')*2;
x_CI =[N, fliplr(N)];
y_CI=[mean(RMSE_dx(:,:,SNR_jj),2)'+CI(:)', flipud(mean(RMSE_dx(:,:,SNR_jj),2)-CI(:))'];
fill(x_CI, y_CI, 1,'facecolor', 'b', 'edgecolor', 'none', 'facealpha', 0.4,'HandleVisibility','off'); hold on
plot(N,mean(RMSE_dx(:,:,SNR_jj),2),'b-s','LineWidth',1.5,'MarkerSize',10); hold on
box_vars = findall(gca,'Tag','Box');
% No Noise
SNR_jj=5;
CI=std(RMSE_dx(:,:,SNR_jj)')*2;
x_CI =[N, fliplr(N)];
y_CI=[mean(RMSE_dx(:,:,SNR_jj),2)'+CI(:)', flipud(mean(RMSE_dx(:,:,SNR_jj),2)-CI(:))'];
fill(x_CI, y_CI, 1,'facecolor',"#D95319" , 'edgecolor', 'none', 'facealpha', 0.4,'HandleVisibility','off'); hold on
plot(N,mean(RMSE_dx(:,:,SNR_jj),2),'-*','Color',"#D95319",'LineWidth',1.5,'MarkerSize',10); hold on
ax.XAxis.Scale ="log";
set(gca,'YLim',[0 4e1])
ylabel(' RMSE $\dot x$ ($\%$)')
% xlabel('Dataset Size')
% grid minor

% Save figure
exportgraphics(gcf,'figures/duffing_3D_convergence_noise_b.pdf','ContentType','vector')

% Figure 8c
figure(3)
Xtick=[2e5 1e6 6e6];
clear hLegend
% Noise 30 db
set(gcf,'units','normalized','outerposition',[0 0 1 0.7])
ax = gca;
SNR_jj=1;
CI=std(RMSE_fr(:,:,SNR_jj)')*2;
x_CI =[N, fliplr(N)];
y_CI=[mean(RMSE_fr(:,:,SNR_jj),2)'+CI(:)', flipud(mean(RMSE_fr(:,:,SNR_jj),2)-CI(:))'];
fill(x_CI, y_CI, 1,'facecolor', 'k', 'edgecolor', 'none', 'facealpha', 0.4,'HandleVisibility','off'); hold on
plot(N,mean(RMSE_fr(:,:,SNR_jj),2),'k-x','LineWidth',1.5,'MarkerSize',10); hold on
box_vars = findall(gca,'Tag','Box');
% Noise 40 db
SNR_jj=3;
CI=std(RMSE_fr(:,:,SNR_jj)')*2;
x_CI =[N, fliplr(N)];
y_CI=[mean(RMSE_fr(:,:,SNR_jj),2)'+CI(:)', flipud(mean(RMSE_fr(:,:,SNR_jj),2)-CI(:))'];
fill(x_CI, y_CI, 1,'facecolor', 'r', 'edgecolor', 'none', 'facealpha', 0.4,'HandleVisibility','off'); hold on
plot(N,mean(RMSE_fr(:,:,SNR_jj),2),'r-o','LineWidth',1.5,'MarkerSize',10); hold on
box_vars = findall(gca,'Tag','Box');
% Noise 45 db
SNR_jj=4;
CI=std(RMSE_fr(:,:,SNR_jj)')*2;
x_CI =[N, fliplr(N)];
y_CI=[mean(RMSE_fr(:,:,SNR_jj),2)'+CI(:)', flipud(mean(RMSE_fr(:,:,SNR_jj),2)-CI(:))'];
fill(x_CI, y_CI, 1,'facecolor', 'b', 'edgecolor', 'none', 'facealpha', 0.4,'HandleVisibility','off'); hold on
plot(N,mean(RMSE_fr(:,:,SNR_jj),2),'b-s','LineWidth',1.5,'MarkerSize',10); hold on
box_vars = findall(gca,'Tag','Box');
% No Noise 
SNR_jj=5;
CI=std(RMSE_fr(:,:,SNR_jj)')*2;
x_CI =[N, fliplr(N)];
y_CI=[mean(RMSE_fr(:,:,SNR_jj),2)'+CI(:)', flipud(mean(RMSE_fr(:,:,SNR_jj),2)-CI(:))'];
fill(x_CI, y_CI, 1,'facecolor',"#D95319" , 'edgecolor', 'none', 'facealpha', 0.4,'HandleVisibility','off'); hold on
plot(N,mean(RMSE_fr(:,:,SNR_jj),2),'-*','Color',"#D95319",'LineWidth',1.5,'MarkerSize',10); hold on
ax.XAxis.Scale ="log";
% set(gca,'YLim',[0 4e1], ...
%         'XLim',[3.5e5 0.5e7],'XTick',10.^(5:7))
set(gca,'YLim',[0 4e1])
ylabel(' RMSE $f_r$ ($\%$)')
xlabel('Dataset Size')

% Save figure
exportgraphics(gcf,'figures/duffing_3D_convergence_noise_c.pdf','ContentType','vector')

