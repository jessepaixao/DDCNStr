%% Data-driven Computing in Nonlinear Dynamics Systems 
%  Description: Experimental demonstration case for the nonlinear vibration
%  absorber
%  Versions:
%   1.0 - 15/07/24 Jessé Paixão

close all
clear all
clc

% return

% Load library of functions
addpath('lib');
addpath('data')

% Figure settings
set(groot,'defaultAxesFontSize',20)
set(groot,'defaultTextInterpreter','latex')
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');

% return

%% DATA IMPORT

file_name='data/25_04_24_exp_nonlinear_absorber.mat'; 


data=load(file_name);
acc_NVA=data.('exp').Y(1).Data;
acc_str=data.('exp').Y(2).Data;
fa_bob_NVA=(data.('exp').Y(3).Data/6.3)*10.3;; % Conversion factor;
force_shaker=data.('exp').Y(4).Data;
vrel=data.('exp').Y(5).Data;
time_vector=data.('exp').X(1).Data;

% Get sampling frequency
dt=2e-4;
Fs=1/dt;


% Dataset signals
ind1=find(time_vector<=170 & time_vector>=0);
% ind1=find(time_vector>=0);
% ind1=find(temps<=100);
tspan_dataset=time_vector(ind1);
force_dataset=fa_bob_NVA(ind1);
acc_NVA_dataset=acc_NVA(ind1);

[b,a]=butter(2,[3 50]/(Fs/2),'Bandpass'); % Design filter
% [b,a]=butter(2,[200]/(Fs/2)); % Design filter
force_dataset=filtfilt(b,a,force_dataset)'; % Filtering
% force_dataset=trapeze(force_dataset,1/Fs,0)'; % Filtering
ddx_dataset=filtfilt(b,a,acc_NVA_dataset); % Filtering
dx_dataset=trapeze(ddx_dataset,1/Fs,0); % Integration
dx_dataset=filtfilt(b,a,dx_dataset); % Filtering
x_dataset=trapeze(dx_dataset,1/Fs,0); %Integration
x_dataset=filtfilt(b,a,x_dataset); % Filtering
x_dataset=x_dataset(1+Fs/2:1:end-Fs/2)';dx_dataset=dx_dataset(1+Fs/2:end-Fs/2)';ddx_dataset=ddx_dataset(1+Fs/2:1:end-Fs/2)';
% x_dataset=x_dataset';dx_dataset=dx_dataset';ddx_dataset=ddx_dataset';

% [b,a]=butter(2,[4 100]/(Fs/2),'Bandpass'); % Design filter
% force_dataset=filtfilt(b,a,force_dataset)'; % Filtering
% % ddx_dataset=filtfilt(b,a,acc_NVA_dataset); % Filtering
% ddx_dataset=detrend(acc_NVA_dataset);
% dx_dataset=trapeze(acc_NVA_dataset,1/Fs,0); % Integration
% % dx_dataset=filtfilt(b,a,dx_dataset); % Filtering
% dx_dataset=detrend(dx_dataset);
% x_dataset=trapeze(dx_dataset,1/Fs,0); %Integration
% % x_dataset=detrend(x_dataset);
% x_dataset=filtfilt(b,a,x_dataset); % Filtering
% x_dataset=x_dataset(1+Fs:1:end-2*Fs)';dx_dataset=dx_dataset(1+Fs:end-2*Fs)';ddx_dataset=ddx_dataset(1+Fs:1:end-2*Fs)';
% 


force_dataset=force_dataset(1+Fs/2:1:end-Fs/2);
tspan_dataset=tspan_dataset(1+Fs/2:1:end-Fs/2);

clear ind1 b a
% Test signals
ind1=find(time_vector<=170 & time_vector>=115);
tspan_test=time_vector(ind1);
force_test=fa_bob_NVA(ind1);
acc_NVA_test=acc_NVA(ind1);

% 
[b,a]=butter(2,[3 50]/(Fs/2),'Bandpass'); % Design filter
% [b,a]=butter(2,[100]/(Fs/2)); % Design filter
force_test=filtfilt(b,a,force_test)'; % Filtering
ddx_test=filtfilt(b,a,acc_NVA_test); % Filtering
dx_test=trapeze(ddx_test,1/Fs,0); % Integration
dx_test=filtfilt(b,a,dx_test); % Filtering
x_test=trapeze(dx_test,1/Fs,0); % Integration
x_test=filtfilt(b,a,x_test); % Filtering
x_test=x_test(1+Fs/2:1:end-Fs/2)';dx_test=dx_test(1+Fs/2:1:end-Fs/2)';ddx_test=ddx_test(1+Fs/2:1:end-Fs/2)';
% % 
force_test=force_test(1+Fs/2:1:end-Fs/2);
tspan_test=tspan_test(1+Fs/2:1:end-Fs/2);
% tspan_test=tspan_test-tspan_test(1);
% x_test=x_test';dx_test=dx_test';ddx_test=ddx_test';


% plot(acc_NVA_test);hold on;plot(ddx_test)
% plot(acc_NVA_dataset);hold on;plot(ddx_dataset)
% plot(x_dataset)

%% PLOT TIME SINGALS
aux=tspan_dataset(tspan_dataset<=70);
t_hline_dataset=aux(end);

% FIGURE 11a
figure()
set(gcf,'units','normalized','outerposition',[0 0 1 0.9])
subplot(311)
plot(tspan_dataset,x_dataset,'r-','LineWidth',1); hold on
xline(tspan_test(1),'--',{'Test','Data'},'LabelOrientation','horizontal','LabelHorizontalAlignment','left','FontSize',18,'Interpreter','latex')
xline(tspan_dataset(tspan_dataset==t_hline_dataset),'--',{'Dataset'},'LabelOrientation','horizontal','FontSize',18,'Interpreter','latex')
annotation('arrow',[0.64931650893796 0.661671924290221]+0.01,...
    [0.91556976744186 0.91453488372093]);
annotation('arrow',[0.448343848580442 0.437828601472134],...
    [0.914406976744186 0.913081395348837]);
% plot(tspan_test,x_test,'k-');hold on
ylabel('$x (m)$')
xlim([0 tspan_dataset(end)])
% legend('Dataset','Location','northeast')
subplot(312)
% plot(tspan_dataset,x_dataset,'r-'); hold on
plot(tspan_dataset,dx_dataset,'r-','LineWidth',1);hold on
xline(tspan_test(1),'--')
xline(tspan_dataset(tspan_dataset==t_hline_dataset),'--')
ylabel('$\dot{x} (m/s)$')
xlim([0 tspan_dataset(end)])
subplot(313)
plot(tspan_dataset,ddx_dataset,'r-','LineWidth',1); hold on
% plot(tspan_test,dx_test,'k-');hold on
xline(tspan_test(1),'--')
xline(tspan_dataset(tspan_dataset==t_hline_dataset),'--')
xlim([0 tspan_dataset(end)])
ylabel('$\ddot{x}$ ($m/s^2$)')
xlabel('t(s)')

% Save figure
exportgraphics(gcf,'figures/NVA_time_signals.pdf',"Resolution",600)


%% PLOT DENSITY OF POINTS

% FIGURE 13
figure()
set(gcf,'units','normalized','outerposition',[0 0 1 1])
h=scatterhist(x_dataset,dx_dataset,'MarkerSize',0.5,'Color','r'); hold on
xlabel('$x$ (m)','FontSize',24)
ylabel('$\dot{x}$ (m/s)','FontSize',24)

binscatter(x_dataset,dx_dataset,250,'FaceAlpha',0.8); hold on
cmp_back=colormap(gca,'sky')
cmp_init=[1 0 0];
cmp_max=[0.95    0.95    0.95];
dt_cmp=(cmp_max-cmp_init)/length(cmp_back);
cmp=[];
for ii=1:length(cmp_back)
    cmp=[cmp; cmp_init+ii*dt_cmp];
end
set(gca, 'ColorScale', 'log');
cb=colorbar();
colormap(gca,flip(cmp));
ylabel(cb,'Bin Counts','FontSize',24,'Rotation',270,'Interpreter','Latex')

% Save figure
exportgraphics(gcf,'figures/NVA_density_plot.pdf','ContentType','image',"Resolution",600)


%% PREDICTION ON PHYSICAL BASIS

m=1;
u_factor=-3.684024259131608;
fr_dataset=force_dataset*u_factor-m*ddx_dataset;
fr_test=force_test*u_factor-m*ddx_test;

% FIGURE 10b
figure(2)
subplot(1,2,1)
set(gcf,'units','normalized','outerposition',[0 0 1 0.8])
plot3(x_dataset,dx_dataset,fr_dataset,'ro','MarkerFaceColor','r','MarkerSize',0.5); hold on
[h,icons]=legend('Dataset','Location','northeast')
xlabel('$x$(m)')
ylabel('$\dot{x}$(m/s)')
zlabel('$f_r$ (N)')
icons = findobj(icons,'Type','line');
set(icons(1:2),'MarkerSize',2);
subplot(1,2,2)
plot(x_dataset,fr_dataset,'ro','MarkerFaceColor','r','MarkerSize',0.5); hold on
xlabel('$x$(m)')
ylabel('$f_r$ (N)')

% Save figure
exportgraphics(gcf,'figures/NVA_fr_dataset.pdf','ContentType','image',"Resolution",600)

%% DATA-DRIVEN SOLVER (CONSTANT HYPERPARAMETER)

Z_star=[x_dataset dx_dataset fr_dataset];

u_test=force_test*u_factor;
Z0_test=[x_test(2) dx_test(2) fr_test(2)];

% Define parameters of Newmark
datadriven_solver.newmark_gamma=0.5;
datadriven_solver.newmark_beta=0.25;

datadriven_solver.Fs=Fs;


% Define constants of objective function
datadriven_solver.A1=10*2/max(x_dataset)^2;
datadriven_solver.A2=1*2/max(dx_dataset)^2;
datadriven_solver.A3=2/(max(fr_dataset))^2;
datadriven_solver.m=m;

datadriven_solver.conv_tol=1e-15;
datadriven_solver.max_steps=100;

% Run the data-driven solver
[Z_rec,opt_data,elapsedTime] = DDCNSys_pred(Z_star,u_test,Z0_test,datadriven_solver);

x_rec=Z_rec(:,1);
dx_rec=Z_rec(:,2);
ddx_rec=Z_rec(:,3);
fr_rec=Z_rec(:,4);

RMSE_x=rms(x_rec-x_test)/rms(x_test)*100
RMSE_dx=rms(dx_rec-dx_test)/rms(dx_test)*100
RMSE_fr=rms(fr_rec-fr_test)/rms(fr_test)*100

% Save processed data
% save(strcat('results\duffing_3D_conv_N',num2str(N_dataset),'.mat'))
% save(strcat('results\case_2_swept_sine_25_T2_A200_',num2str(N_dataset),'.mat'))
% fileName = ['results/' datestr(now, 'dd-mmm-yyyy_HHMMSS') '_exp_case_NVA']
% save(fileName)

%% PLOT TESTE DATASET X RECONSTRUCTED 
close all
t_zoom=3;
jj=length(x_rec)-1;

% FIGURE 11a
figure()
set(gcf,'units','normalized','outerposition',[0 0 1 1])
% Subplot x
subplot(4,3,[1 2])
plot(tspan_test(1:jj+1),x_test(1:jj+1),'k-','MarkerSize',10,'Linewidth',1.5,'MarkerSize',12); hold on
plot(tspan_test(1:jj+1),x_rec(1:jj+1),'b--','MarkerSize',10,'Linewidth',1,'MarkerSize',12); hold on
plot(tspan_test(1:jj+1),(x_rec(1:jj+1)-x_test(1:jj+1)),'g-','MarkerSize',10,'Linewidth',1,'MarkerSize',12); hold on
ylabel('$x$(m)')
xlim([tspan_test(1) tspan_test(end)])
subplot(4,3,[3])
plot(tspan_test(1:jj+1),x_test(1:jj+1),'k-','MarkerSize',10,'Linewidth',1.5,'MarkerSize',12); hold on
plot(tspan_test(1:jj+1),x_rec(1:jj+1),'b--','MarkerSize',10,'Linewidth',1,'MarkerSize',12); hold on
plot(tspan_test(1:jj+1),(x_rec(1:jj+1)-x_test(1:jj+1)),'g-','MarkerSize',10,'Linewidth',1,'MarkerSize',12); hold on
xlim([tspan_test(1) tspan_test(1)+t_zoom])
title('Zoom')
legend('Reference Sol.','Data-driven Pred.','Abs. Error','Position',[0.855581709990948 0.879909955777893 0.122161016963136 0.106869838878143])
% Subplot dx
subplot(4,3,[4,5])
plot(tspan_test(1:jj+1),dx_test(1:jj+1),'k-','MarkerSize',10,'Linewidth',1.5,'MarkerSize',12); hold on
plot(tspan_test(1:jj+1),dx_rec(1:jj+1),'b--','MarkerSize',10,'Linewidth',1,'MarkerSize',12); hold on
plot(tspan_test(1:jj+1),(dx_rec(1:jj+1)-dx_test(1:jj+1)),'g-','MarkerSize',10,'Linewidth',1,'MarkerSize',12); hold on
ylabel('$\dot{x}$(m/s)')
xlim([tspan_test(1) tspan_test(end)])
subplot(4,3,[6])
plot(tspan_test(1:jj+1),dx_test(1:jj+1),'k-','MarkerSize',10,'Linewidth',1.5,'MarkerSize',12); hold on
plot(tspan_test(1:jj+1),dx_rec(1:jj+1),'b--','MarkerSize',10,'Linewidth',1,'MarkerSize',12); hold on
plot(tspan_test(1:jj+1),(dx_rec(1:jj+1)-dx_test(1:jj+1)),'g-','MarkerSize',10,'Linewidth',1,'MarkerSize',12); hold on
% ylabel('$\dot{x}$(m/s)')
xlim([tspan_test(1) tspan_test(1)+t_zoom])
% Subplot fr
subplot(4,3,[7,8])
plot(tspan_test(1:jj+1),fr_test(1:jj+1),'k-','MarkerSize',10,'Linewidth',1.5,'MarkerSize',12); hold on
plot(tspan_test(1:jj+1),fr_rec(1:jj+1),'b--','MarkerSize',10,'Linewidth',1,'MarkerSize',12); hold on
plot(tspan_test(1:jj+1),(fr_rec(1:jj+1)-fr_test(1:jj+1)),'g-','MarkerSize',10,'Linewidth',1,'MarkerSize',12); hold on
ylabel('$f_r$ (N)')
xlim([tspan_test(1) tspan_test(end)])
subplot(4,3,[9])
plot(tspan_test(1:jj+1),fr_test(1:jj+1),'k-','MarkerSize',10,'Linewidth',1.5,'MarkerSize',12); hold on
plot(tspan_test(1:jj+1),fr_rec(1:jj+1),'b--','MarkerSize',10,'Linewidth',1,'MarkerSize',12); hold on
plot(tspan_test(1:jj+1),(fr_rec(1:jj+1)-fr_test(1:jj+1)),'g-','MarkerSize',10,'Linewidth',1,'MarkerSize',12); hold on
xlim([tspan_test(1) tspan_test(1)+t_zoom])
% Subplot iterations
subplot(4,3,[10,11])
bar(tspan_test(2:jj+1),vertcat(opt_data(:).steps),'b'); hold on
xlabel('t(s)')
ylabel('$\#Iterations$')
xlim([tspan_test(1) tspan_test(end)])
subplot(4,3,[12])
bar(tspan_test(2:jj+1),vertcat(opt_data(:).steps),'b'); hold on
xlabel('t(s)')
xlim([tspan_test(1) tspan_test(1)+t_zoom])

% Save figure
exportgraphics(gcf,'figures/NVA_pred_time.pdf','ContentType','image','Resolution',600)


% FIGURE 11b
% Plot restoring force
figure(2)
set(gcf,'units','normalized','outerposition',[0 0 1 0.8])
plot3(x_dataset,dx_dataset,fr_dataset,'ro','MarkerFaceColor','r','MarkerSize',0.5); hold on
plot3(x_test,dx_test,fr_test,'k.','MarkerSize',10); hold on
plot3(x_rec,dx_rec,fr_rec,'bs','MarkerFaceColor','b','MarkerSize',4)
xlabel('$x$(m)')
ylabel('$\dot{x}$(m/s)')
zlabel('$f_r$ (N)')
% grid minor
% legend('Dataset','Reference Solution','Data-driven Pred.','Position',[0.805581709990948 0.879909955777893 0.122161016963136 0.106869838878143])
[h,icons]=legend('Dataset','Test Data (Ref.)','Data-driven Pred.','Position',[0.780853309541509 0.825465884143214 0.162154095464538 0.151928194913458])
icons = findobj(icons,'Type','line');
set(icons(1:2),'MarkerSize',4);

% Save figure
% exportgraphics(gcf,'figures/NVA_pred_rfs.pdf','ContentType','image','Resolution',600)