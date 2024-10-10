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
Xtick=[2e5 1e6 6e6];
clear hLegend
set(gcf,'units','normalized','outerposition',[0 0 1 0.7])
ax = gca;
% Noise 30 db
SNR_jj=1;
h=boxplot(RMSE_x(:,:,SNR_jj)','Positions',N,'PlotStyle','compact','Colors','k','symbol', '');hold on
xticks(Xtick)
xticklabels(num2cell(Xtick))
xtickformat('%1.e')
plot(N,mean(RMSE_x(:,:,SNR_jj),2),'--k','LineWidth',1.5); hold on
box_vars = findall(gca,'Tag','Box');
hLegend=box_vars;
% Noise 40 db
SNR_jj=3;
h=boxplot(RMSE_x(:,:,SNR_jj)','Positions',N,'PlotStyle','compact','Colors','r','symbol', '');hold on
box_vars = findall(gca,'Tag','Box');
hLegend=box_vars;
xticks(Xtick)
xticklabels(num2cell(Xtick))
xtickformat('%1.e')
plot(N,mean(RMSE_x(:,:,SNR_jj),2),'--r','LineWidth',1.5); hold on
% Noise 45 db
SNR_jj=4;
h=boxplot(RMSE_x(:,:,SNR_jj)','Positions',N,'PlotStyle','compact','Colors','b','symbol', '');hold on
box_vars = findall(gca,'Tag','Box');
hLegend=box_vars;
xticks(Xtick)
xticklabels(num2cell(Xtick))
xtickformat('%1.e')
plot(N,mean(RMSE_x(:,:,SNR_jj),2),'--b','LineWidth',1.5); hold on
% No Noise 
SNR_jj=5;
h=boxplot(RMSE_x(:,:,SNR_jj)','Positions',N,'PlotStyle','compact','Colors','m','symbol', '');hold on
box_vars = findall(gca,'Tag','Box');
hLegend=box_vars;
xticks(Xtick)
xticklabels(num2cell(Xtick))
xtickformat('%1.e')
plot(N,mean(RMSE_x(:,:,SNR_jj),2),'--m'); hold on
ax.XAxis.Scale ="log";
set(gca,'YLim',[0 4e1], ...
        'XLim',[3e5 1e7],'XTick',10.^(5:7))
ylabel(' RMSE $x$ ($\%$)')
% xlabel('Dataset Size')
% grid minor
legend(hLegend(1:4:16), {'SNR 30 dB','SNR 40 dB','SNR 45 dB','No Noise'})

% Save figure
exportgraphics(gcf,'figures/duffing_3D_convergence_noise_a.pdf','ContentType','vector')


% Figure 8b
figure(2)
Xtick=[2e5 1e6 6e6];
clear hLegend
set(gcf,'units','normalized','outerposition',[0 0 1 0.7])
ax = gca;
% Noise 30 db
SNR_jj=1;
h=boxplot(RMSE_dx(:,:,SNR_jj)','Positions',N,'PlotStyle','compact','Colors','k','symbol', '');hold on
xticks(Xtick)
xticklabels(num2cell(Xtick))
xtickformat('%1.e')
plot(N,mean(RMSE_dx(:,:,SNR_jj),2),'--k','LineWidth',1.5); hold on
box_vars = findall(gca,'Tag','Box');
hLegend=box_vars;
% Noise 40 db
SNR_jj=3;
h=boxplot(RMSE_dx(:,:,SNR_jj)','Positions',N,'PlotStyle','compact','Colors','r','symbol', '');hold on
box_vars = findall(gca,'Tag','Box');
hLegend=box_vars;
xticks(Xtick)
xticklabels(num2cell(Xtick))
xtickformat('%1.e')
plot(N,mean(RMSE_dx(:,:,SNR_jj),2),'--r','LineWidth',1.5); hold on
% Noise 45 db
SNR_jj=4;
h=boxplot(RMSE_dx(:,:,SNR_jj)','Positions',N,'PlotStyle','compact','Colors','b','symbol', '');hold on
box_vars = findall(gca,'Tag','Box');
hLegend=box_vars;
xticks(Xtick)
xticklabels(num2cell(Xtick))
xtickformat('%1.e')
plot(N,mean(RMSE_dx(:,:,SNR_jj),2),'--b','LineWidth',1.5); hold on
% No Noise
SNR_jj=5;
h=boxplot(RMSE_dx(:,:,SNR_jj)','Positions',N,'PlotStyle','compact','Colors','m','symbol', '');hold on
box_vars = findall(gca,'Tag','Box');
hLegend=box_vars;
xticks(Xtick)
xticklabels(num2cell(Xtick))
xtickformat('%1.e')
plot(N,mean(RMSE_dx(:,:,SNR_jj),2),'--m'); hold on
ax.XAxis.Scale ="log";
set(gca,'YLim',[0 4e1], ...
        'XLim',[3e5 1e7],'XTick',10.^(5:7))
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
h=boxplot(RMSE_fr(:,:,SNR_jj)','Positions',N,'PlotStyle','compact','Colors','k','symbol', '');hold on
xticks(Xtick)
xticklabels(num2cell(Xtick))
xtickformat('%1.e')
plot(N,mean(RMSE_fr(:,:,SNR_jj),2),'--k','LineWidth',1.5); hold on
box_vars = findall(gca,'Tag','Box');
hLegend=box_vars;
% Noise 40 db
SNR_jj=3;
h=boxplot(RMSE_fr(:,:,SNR_jj)','Positions',N,'PlotStyle','compact','Colors','r','symbol', '');hold on
box_vars = findall(gca,'Tag','Box');
hLegend=box_vars;
xticks(Xtick)
xticklabels(num2cell(Xtick))
xtickformat('%1.e')
plot(N,mean(RMSE_fr(:,:,SNR_jj),2),'--r','LineWidth',1.5); hold on
% Noise 45 db
SNR_jj=4;
h=boxplot(RMSE_fr(:,:,SNR_jj)','Positions',N,'PlotStyle','compact','Colors','b','symbol', '');hold on
box_vars = findall(gca,'Tag','Box');
hLegend=box_vars;
xticks(Xtick)
xticklabels(num2cell(Xtick))
xtickformat('%1.e')
plot(N,mean(RMSE_fr(:,:,SNR_jj),2),'--b','LineWidth',1.5); hold on
% No Noise 
SNR_jj=5;
h=boxplot(RMSE_fr(:,:,SNR_jj)','Positions',N,'PlotStyle','compact','Colors','m','symbol', '');hold on
box_vars = findall(gca,'Tag','Box');
hLegend=box_vars;
xticks(Xtick)
xticklabels(num2cell(Xtick))
xtickformat('%1.e')
plot(N,mean(RMSE_fr(:,:,SNR_jj),2),'--m'); hold on
ax.XAxis.Scale ="log";
set(gca,'YLim',[0 4e1], ...
        'XLim',[3e5 1e7],'XTick',10.^(5:7))
ylabel(' RMSE $f_r$ ($\%$)')
xlabel('Dataset Size')

% Save figure
exportgraphics(gcf,'figures/duffing_3D_convergence_noise_c.pdf','ContentType','vector')

