%% Data-driven Computing in Nonlinear Dynamics Systems 
%  Description: Numerical demonstration case in a Duffing oscillator
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

mm=0;

SNR_vec=[40 50 60];

for T=[50 100 500 1000]
    T
mm=mm+1;
pp=0;
for SNR=[40 50 60]
    SNR
    pp=pp+1;
    

close all
rng(123,'twister')

% Aqcquisition parameters
Fs=5000;                     % Sampling frequency

% % 2) Continuous Sweept-sine Exciation
% Excitation to generate training data
RMSu = 50;                        % RMS Amplitude
tspan_dataset  = [0:1/Fs:T]';      % Time vector
N_dataset = length(tspan_dataset);    % Number of samples
T_exc = tspan_dataset(end);         % Period of swept
f_ini = 1;                        % Initial frequency
f_end = 20;                       % Final frequency
u_exc=sin(2*pi*(f_ini*tspan_dataset+(f_end-f_ini)/(2*T_exc)*tspan_dataset.^2));
% u_exc=sin(2*pi*f_ini*tspan_train);
u_exc_dataset = u_exc/rms(u_exc(:,1))*RMSu;

% Excitation to generate test data
% RMSu =50;                   % RMS Amplitude
RMSu =60;                   % RMS Amplitude
tspan_test  = [0:1/Fs:2]';   % Time vector
N_test = length(tspan_test); % Number of sample
T_exc = tspan_test(end);     % Period of swept
f_ini = 5;                  % Initial frequency
f_end = 20;                  % Final frequency
f_exc = 24; 
% u_exc=sin(2*pi*f_exc*tspan_test);
u_exc=cos(2*pi*(f_ini*tspan_test+(f_end-f_ini)/(2*T_exc)*tspan_test.^2));
u_exc_test = u_exc/rms(u_exc(:,1))*RMSu;
% u_exc_test=u_exc_train(1:N_test);

% NUMERICAL SIMULATION OF DYNAMICAL SYSTEM

m = 1;                      % Mass [kg]
wn=2*pi*10;                 % Linear natural frequency 
k1 = wn^2*m;                % Linear stiffness [N/m]
k2 = 0;                     % Quadratic stiffness [N/m^2]
k3 = 1e5;                   % Cubic stiffness [N/m^3]
c = 4;                      % Damping [Ns/m]
eta=c/(2*m*(wn))*100;       % Damping ratio [%]

xc0=[0 0];  % Initial condition [displacement, velocity]
Gamma=0.5;  % Newmark parameter (1/2)
Beta=0.25;  % Newmark parameter (1/4 - average acc // 1/6 - linear acc) 
tol=1e-8;   % Newmark convergence tolerance

[x_dataset,dx_dataset,ddx_dataset] = newmark_duffing(m,c,k1,k2,k3,u_exc_dataset,tspan_dataset,xc0,Gamma,Beta,tol);
[x_test,dx_test,ddx_test] = newmark_duffing(m,c,k1,k2,k3,u_exc_test,tspan_test,xc0,Gamma,Beta,tol);

% % % Noise
% SNR=SNR_vec(pp);                     % Signal-to-noise ratio
sig=x_dataset;
noisy = randn(size(sig,1),1).*std(sig)/db2mag(SNR);
x_dataset=x_dataset+noisy;

sig=dx_dataset;
noisy = randn(size(sig,1),1).*std(sig)/db2mag(SNR);
dx_dataset=dx_dataset+noisy;

sig=ddx_dataset;
noisy = randn(size(sig,1),1).*std(sig)/db2mag(SNR);
ddx_dataset=ddx_dataset+noisy;

fr_dataset=u_exc_dataset-m*ddx_dataset;
fr_test=u_exc_test-m*ddx_test;


% DATA-DRIVEN SOLVER (CONSTANT HYPERPARAMETER)

Z_star=[x_dataset dx_dataset fr_dataset];

u_test=u_exc_test;
Z0_test=[x_test(1) dx_test(1) fr_test(1)];

% Define parameters of Newmark
datadriven_solver.newmark_gamma=0.5;
datadriven_solver.newmark_beta=0.25;

datadriven_solver.Fs=Fs;


% Define constants of objective function
datadriven_solver.A1=2/max(x_dataset)^2;
datadriven_solver.A2=2/max(dx_dataset)^2;
datadriven_solver.A3=2/(max(fr_dataset))^2;
% datadriven_solver.A1=1;
% datadriven_solver.A2=1;
% datadriven_solver.A3=1;
datadriven_solver.m=m;

datadriven_solver.conv_tol=1e-10;
datadriven_solver.max_steps=100;

for rr=11:30

% Run the data-driven solver
[Z_rec,opt_data,elapsedTime] = DDCNSys_pred(Z_star,u_test,Z0_test,datadriven_solver);

x_rec=Z_rec(:,1);
dx_rec=Z_rec(:,2);
ddx_rec=Z_rec(:,3);
fr_rec=Z_rec(:,4);

RMSE_x(mm,rr,pp)=rms(x_rec-x_test)/rms(x_test)*100
RMSE_dx(mm,rr,pp)=rms(dx_rec-dx_test)/rms(dx_test)*100
RMSE_fr(mm,rr,pp)=rms(fr_rec-fr_test)/rms(fr_test)*100

save(strcat('results\duffing_3D_conv_N',num2str(N_dataset),'_R_',num2str(rr),'_SNR_',num2str(SNR),'.mat'))
end

end
end


% save(strcat('results\case_2_swept_sine_25_T2_A200_',num2str(N_dataset),'.mat'))

%%
mm=0;

SNR_vec=[0 45 40 35];
% T_vec=[80 200 600 1000];
T_vec=[80 200];
R=10;

for T=T_vec
    T
    mm=mm+1;
    pp=0;
    
    % close all
    % rng(123,'twister')
    
    % Aqcquisition parameters
    Fs=5000;                     % Sampling frequency
    
    % % 2) Continuous Sweept-sine Exciation
    % Excitation to generate training data
    RMSu = 40;                        % RMS Amplitude
    tspan_dataset  = [0:1/Fs:T]';      % Time vector
    N_dataset = length(tspan_dataset);    % Number of samples
    T_exc = tspan_dataset(end);         % Period of swept
    f_ini = 1;                        % Initial frequency
    f_end = 25;                       % Final frequency
    u_exc=sin(2*pi*(f_ini*tspan_dataset+(f_end-f_ini)/(2*T_exc)*tspan_dataset.^2));
    % u_exc=sin(2*pi*f_ini*tspan_train);
    u_exc_dataset = u_exc/rms(u_exc(:,1))*RMSu;
    
    % Excitation to generate test data
    RMSu =60;                   % RMS Amplitude
    % RMSu =60;                   % RMS Amplitude
    tspan_test  = [0:1/Fs:2]';   % Time vector
    N_test = length(tspan_test); % Number of sample
    T_exc = tspan_test(end);     % Period of swept
    f_ini = 5;                  % Initial frequency
    f_end = 25;                  % Final frequency
    f_exc = 24; 
    % u_exc=sin(2*pi*f_exc*tspan_test);
    u_exc=cos(2*pi*(f_ini*tspan_test+(f_end-f_ini)/(2*T_exc)*tspan_test.^2));
    u_exc_test = u_exc/rms(u_exc(:,1))*RMSu;
    % u_exc_test=u_exc_train(1:N_test);
    
    % NUMERICAL SIMULATION OF DYNAMICAL SYSTEM
    
    m = 1;                      % Mass [kg]
    wn=2*pi*10;                 % Linear natural frequency 
    k1 = wn^2*m;                % Linear stiffness [N/m]
    k2 = 0;                     % Quadratic stiffness [N/m^2]
    k3 = 1e5;                   % Cubic stiffness [N/m^3]
    c = 4;                      % Damping [Ns/m]
    eta=c/(2*m*(wn))*100;       % Damping ratio [%]
    
    xc0=[0 0];  % Initial condition [displacement, velocity]
    Gamma=0.5;  % Newmark parameter (1/2)
    Beta=0.25;  % Newmark parameter (1/4 - average acc // 1/6 - linear acc) 
    tol=1e-8;   % Newmark convergence tolerance
    
    [x_dataset,dx_dataset,ddx_dataset] = newmark_duffing(m,c,k1,k2,k3,u_exc_dataset,tspan_dataset,xc0,Gamma,Beta,tol);
    [x_test,dx_test,ddx_test] = newmark_duffing(m,c,k1,k2,k3,u_exc_test,tspan_test,xc0,Gamma,Beta,tol);
    
    
    fr_dataset=u_exc_dataset-m*ddx_dataset;
    fr_test=u_exc_test-m*ddx_test;
    
    
    
    % DATA-DRIVEN SOLVER (CONSTANT HYPERPARAMETER)
       
    % Define parameters of Newmark
    datadriven_solver.newmark_gamma=0.5;
    datadriven_solver.newmark_beta=0.25;
    
    datadriven_solver.Fs=Fs;
    
    
    % Define constants of objective function
    datadriven_solver.A1=2/max(x_dataset)^2;
    datadriven_solver.A2=2/max(dx_dataset)^2;
    datadriven_solver.A3=2/(max(fr_dataset))^2;
    % datadriven_solver.A1=1;
    % datadriven_solver.A2=1;
    % datadriven_solver.A3=1;
    datadriven_solver.m=m;
    
    datadriven_solver.conv_tol=1e-10;
    datadriven_solver.max_steps=100;

    parfor pp=1:length(SNR_vec)
        SNR
        % pp  =pp+1;

        for rr=1:R
    
        % % % Noise
        % SNR=SNR_vec(pp);                     % Signal-to-noise ratio
        if SNR_vec(pp)~=0
        sig=x_dataset;
        noisy = randn(size(sig,1),1).*std(sig)/db2mag(SNR);
        x_dataset_noisy=x_dataset+noisy;
        
        sig=dx_dataset;
        noisy = randn(size(sig,1),1).*std(sig)/db2mag(SNR);
        dx_dataset_noisy=dx_dataset+noisy;
        
        sig=ddx_dataset;
        noisy = randn(size(sig,1),1).*std(sig)/db2mag(SNR);
        ddx_dataset_noisy=ddx_dataset+noisy;
    
    
    
        else
    
            x_dataset_noisy=x_dataset;
            dx_dataset_noisy=dx_dataset;
            ddx_dataset_noisy=ddx_dataset;
    
        end
        
        fr_dataset_noisy=u_exc_dataset-m*ddx_dataset_noisy;
    
    
        Z_star=[x_dataset_noisy dx_dataset_noisy fr_dataset_noisy];
        
        u_test=u_exc_test;
        Z0_test=[x_test(1) dx_test(1) fr_test(1)];
       
        
        % Run the data-driven solver
        [Z_rec,opt_data,elapsedTime] = DDCNSys_pred(Z_star,u_test,Z0_test,datadriven_solver);
        
        x_rec=Z_rec(:,1);
        dx_rec=Z_rec(:,2);
        ddx_rec=Z_rec(:,3);
        fr_rec=Z_rec(:,4);
        
        RMSE_x(mm,rr,pp)=rms(x_rec-x_test)/rms(x_test)*100
        RMSE_dx(mm,rr,pp)=rms(dx_rec-dx_test)/rms(dx_test)*100
        RMSE_fr(mm,rr,pp)=rms(fr_rec-fr_test)/rms(fr_test)*100
    
        elapsedTime_vec(mm,rr,pp)=elapsedTime;
        
        end

    save(strcat('results/duffing_3D_conv_N',num2str(N_dataset),'_SNR_',num2str(SNR_vec(pp)),'.mat'),'RMSE_fr','RMSE_x','RMSE_dx','T_vec','SNR_vec','elapsedTime_vec')
    
    end
end
%% PLOT TESTE DATASET X RECONSTRUCTED 
close all

jj=length(x_rec)-1;

figure()
set(gcf,'units','normalized','outerposition',[0 0 1 1])
subplot(411)
plot(tspan_test(1:jj+1),x_test(1:jj+1),'k-','MarkerSize',10,'Linewidth',1.5,'MarkerSize',12); hold on
plot(tspan_test(1:jj+1),x_rec(1:jj+1),'b--','MarkerSize',10,'Linewidth',1,'MarkerSize',12); hold on
plot(tspan_test(1:jj+1),(x_rec(1:jj+1)-x_test(1:jj+1)),'g-','MarkerSize',10,'Linewidth',1,'MarkerSize',12); hold on
ylabel('$x$(m)')
xlim([0 tspan_test(end)])
legend('Reference Sol.','Data-driven Pred.','Abs. Error','Position',[0.855581709990948 0.879909955777893 0.122161016963136 0.106869838878143])
subplot(412)
plot(tspan_test(1:jj+1),dx_test(1:jj+1),'k-','MarkerSize',10,'Linewidth',1.5,'MarkerSize',12); hold on
plot(tspan_test(1:jj+1),dx_rec(1:jj+1),'b--','MarkerSize',10,'Linewidth',1,'MarkerSize',12); hold on
plot(tspan_test(1:jj+1),(dx_rec(1:jj+1)-dx_test(1:jj+1)),'g-','MarkerSize',10,'Linewidth',1,'MarkerSize',12); hold on
ylabel('$\dot{x}$(m/s)')
xlim([0 tspan_test(end)])
subplot(413)
plot(tspan_test(1:jj+1),fr_test(1:jj+1),'k-','MarkerSize',10,'Linewidth',1.5,'MarkerSize',12); hold on
plot(tspan_test(1:jj+1),fr_rec(1:jj+1),'b--','MarkerSize',10,'Linewidth',1,'MarkerSize',12); hold on
plot(tspan_test(1:jj+1),(fr_rec(1:jj+1)-fr_test(1:jj+1)),'g-','MarkerSize',10,'Linewidth',1,'MarkerSize',12); hold on
xlabel('t(s)')
ylabel('$fr(x,\dot{x})$ (N)')
xlim([0 tspan_test(end)])
subplot(414)
bar(tspan_test(2:jj+1),vertcat(opt_data(:).steps),'b'); hold on
% ylim([0,20])
xlabel('t(s)')
ylabel('$\#Iterations$')
xlim([0 tspan_test(end)])
% exportgraphics(gcf,'figures/duffing_3D_pred_time_data_driven.png',"Resolution",600)

% Plot restoring force
figure(2)
set(gcf,'units','normalized','outerposition',[0 0 0.6 1])
plot3(x_dataset,dx_dataset,fr_dataset,'ro','MarkerSize',1); hold on
plot3(x_test,dx_test,fr_test,'k.','MarkerSize',10); hold on
plot3(x_rec,dx_rec,fr_rec,'bs','MarkerFaceColor','b','MarkerSize',4)
xlabel('$x$(m)')
ylabel('$\dot{x}$(m/s)')
zlabel('$fr(x,\dot{x})$ (N)')
% grid on
legend('Dataset','Reference Sol.','Data-driven Pred.','Position',[0.805581709990948 0.879909955777893 0.122161016963136 0.106869838878143])
% exportgraphics(gcf,'figures/duffing_3D_pred_rfs1_data_driven.png',"Resolution",600)

% % Plot restoring force
% figure(3)
% set(gcf,'units','normalized','outerposition',[0 0 1 1])
% % plot3(x_train,dx_train,fnl_train,'r.'); hold on
% plot(x_dataset,fr_dataset,'ro','MarkerSize',0.5); hold on
% plot(x_test,fr_test,'k.','MarkerSize',10); hold on
% plot(x_rec,fr_rec,'bs','MarkerFaceColor','b','MarkerSize',4)
% xlabel('$x$(m)')
% % ylabel('$\dot{x}$(m/s)')
% ylabel('$fr$ (N)')
% % grid on
% legend('Dataset','Reference Sol.','Data-driven Pred.','Position',[0.855581709990948 0.879909955777893 0.122161016963136 0.106869838878143])
% exportgraphics(gcf,'figures/duffing_3D_pred_rfs2_data_driven.png',"Resolution",600)

%% DATA 

RMSE_x_bakcup=RMSE_x;
RMSE_dx_bakcup=RMSE_dx;
RMSE_fr_bakcup=RMSE_fr;

%%

clear RMSE_x RMSE_dx RMSE_fr

RMSE_x=RMSE_x_bakcup(:,:,1:4);
RMSE_x(2,:,4)=RMSE_x_bakcup(2,:,5);
RMSE_x(3,:,4)=RMSE_x_bakcup(3,:,6);
RMSE_x(4,:,4)=RMSE_x_bakcup(4,:,7);

RMSE_dx=RMSE_dx_bakcup(:,:,1:4);
RMSE_dx(2,:,4)=RMSE_dx_bakcup(2,:,5);
RMSE_dx(3,:,4)=RMSE_dx_bakcup(3,:,6);
RMSE_dx(4,:,4)=RMSE_dx_bakcup(4,:,7);


RMSE_fr=RMSE_fr_bakcup(:,:,1:4);
RMSE_fr(2,:,4)=RMSE_fr_bakcup(2,:,5);
RMSE_fr(3,:,4)=RMSE_fr_bakcup(3,:,6);
RMSE_fr(4,:,4)=RMSE_fr_bakcup(4,:,7);

%% PLOT CONVERGENCE (Fr)

close all

N=[length([0:1/Fs:50]) length([0:1/Fs:100]) length([0:1/Fs:500]) length([0:1/Fs:1000])]

Xtick=[2e5 1e6 6e6];
clear hLegend
figure(1)
set(gcf,'units','normalized','outerposition',[0 0 1 0.7])
ax = gca;

SNR_jj=1;
h=boxplot(RMSE_x(:,:,SNR_jj)','Positions',N,'PlotStyle','compact','Colors','k');hold on
xticks(Xtick)
xticklabels(num2cell(Xtick))
xtickformat('%1.e')
plot(N,median(RMSE_x(:,:,SNR_jj),2),'--k','LineWidth',1.5); hold on
box_vars = findall(gca,'Tag','Box');
hLegend=box_vars;


SNR_jj=2;
h=boxplot(RMSE_x(:,:,SNR_jj)','Positions',N,'PlotStyle','compact','Colors','r');hold on
box_vars = findall(gca,'Tag','Box');
hLegend=box_vars;
xticks(Xtick)
xticklabels(num2cell(Xtick))
xtickformat('%1.e')
plot(N,median(RMSE_x(:,:,SNR_jj),2),'--r','LineWidth',1.5); hold on

SNR_jj=3;
h=boxplot(RMSE_x(:,:,SNR_jj)','Positions',N,'PlotStyle','compact','Colors','b');hold on
box_vars = findall(gca,'Tag','Box');
hLegend=box_vars;
xticks(Xtick)
xticklabels(num2cell(Xtick))
xtickformat('%1.e')
plot(N,median(RMSE_x(:,:,SNR_jj),2),'--b','LineWidth',1.5); hold on

% SNR_jj=4;
% h=boxplot(RMSE_x(:,:,SNR_jj)','Positions',N,'PlotStyle','compact','Colors','g');hold on
% box_vars = findall(gca,'Tag','Box');
% hLegend=box_vars;
% xticks(Xtick)
% xticklabels(num2cell(Xtick))
% xtickformat('%1.e')
% plot(N,median(RMSE_x(:,:,SNR_jj),2),'--g'); hold on


ax.YAxis.Scale ="log";
ax.XAxis.Scale ="log";
set(gca,'YLim',[5e0 1e2],'YTick',10.^(0:2), ...
        'XLim',[2e5 1e7],'XTick',10.^(5:7))
% set(gca,'XTick',[N])
ylabel(' RMSE $x$ ($\%$)')
% xlabel('Dataset Size')
% grid minor


legend(hLegend(1:4:12), {'SNR 60 dB','SNR 50 dB','SNR 40 dB','No Noise'})
% exportgraphics(gcf,'figures/duffing_3D_convergence_noise_a.png',"Resolution",600)

figure(2)
set(gcf,'units','normalized','outerposition',[0 0 1 0.7])
ax = gca;
% set(gcf,'position',[0.1300    0.3020    0.7750    0.6430])

SNR_jj=1;
h=boxplot(RMSE_dx(:,:,SNR_jj)','Positions',N,'PlotStyle','compact','Colors','k');hold on
xticks(Xtick)
xticklabels(num2cell(Xtick))
xtickformat('%1.e')
plot(N,median(RMSE_dx(:,:,SNR_jj),2),'--k','LineWidth',1.5); hold on
box_vars = findall(gca,'Tag','Box');
hLegend=box_vars;

SNR_jj=2;
h=boxplot(RMSE_dx(:,:,SNR_jj)','Positions',N,'PlotStyle','compact','Colors','r');hold on
box_vars = findall(gca,'Tag','Box');
hLegend=box_vars;
xticks(Xtick)
xticklabels(num2cell(Xtick))
xtickformat('%1.e')
plot(N,median(RMSE_dx(:,:,SNR_jj),2),'--r','LineWidth',1.5); hold on

SNR_jj=3;
h=boxplot(RMSE_dx(:,:,SNR_jj)','Positions',N,'PlotStyle','compact','Colors','b');hold on
box_vars = findall(gca,'Tag','Box');
hLegend=box_vars;
xticks(Xtick)
xticklabels(num2cell(Xtick))
xtickformat('%1.e')
plot(N,median(RMSE_dx(:,:,SNR_jj),2),'--b','LineWidth',1.5); hold on


% SNR_jj=4;
% plot(N,median(RMSE_dx(:,:,SNR_jj),2),'--g'); hold on
% 
% h=boxplot(RMSE_dx(:,:,SNR_jj)','Positions',N,'PlotStyle','compact','Colors','g');hold on
% % plot(N,median(RMSE_dx(:,:,SNR_jj),2),'--g'); hold on
% box_vars = findall(gca,'Tag','Box');
% hLegend=box_vars;
% xticks(Xtick)
% xticklabels(num2cell(Xtick))
% xtickformat('%1.e')

% ax.YAxis.Scale ="log";
ax.XAxis.Scale ="log";
set(gca,'YLim',[0 40], ...
        'XLim',[2e5 1e7],'XTick',10.^(5:7))
% set(gca,'XTick',[N])

ylabel(' RMSE $\dot{x}$ ($\%$)')
% xlabel('Dataset Size','Position',[443.4150000000001,-40.000000300407418,0],'VerticalAlignment','top','HorizontalAlignment','center')
% grid minor
% exportgraphics(gcf,'figures/duffing_3D_convergence_noise_b.png',"Resolution",600)

%%
figure(3)
set(gcf,'units','normalized','outerposition',[0 0 1 0.7])
ax = gca;
set(gca,'position',[0.13,0.17658878621654,0.775,0.765934578269442]);    

SNR_jj=1;
loglog(N,median(RMSE_fr(:,:,SNR_jj),2),'--k','LineWidth',1.5); hold on

h=boxplot(RMSE_fr(:,:,SNR_jj)','Positions',N,'PlotStyle','compact','Colors','k');hold on
xticks(Xtick)
xticklabels(num2cell(Xtick))
xtickformat('%1.e')
box_vars = findall(gca,'Tag','Box');
hLegend=box_vars;

SNR_jj=2;
h=boxplot(RMSE_fr(:,:,SNR_jj)','Positions',N,'PlotStyle','compact','Colors','r');hold on
box_vars = findall(gca,'Tag','Box');
hLegend=box_vars;
xticks(Xtick)
xticklabels(num2cell(Xtick))
xtickformat('%1.e')
loglog(N,median(RMSE_fr(:,:,SNR_jj),2),'--r','LineWidth',1.5); hold on
%
SNR_jj=3;
h=boxplot(RMSE_fr(:,:,SNR_jj)','Positions',N,'PlotStyle','compact','Colors','b');hold on
% h=boxchart(N,RMSE_fr(:,:,SNR_jj)');hold on
box_vars = findall(gca,'Tag','Box');
hLegend=box_vars;
xticks(Xtick)
xticklabels(num2cell(Xtick))
xtickformat('%1.e')
loglog(N,median(RMSE_dx(:,:,SNR_jj),2),'--b','LineWidth',1.5); hold on

% SNR_jj=4;
% h=boxplot(RMSE_dx(:,:,SNR_jj)','Positions',N,'PlotStyle','compact','Colors','g');hold on
% box_vars = findall(gca,'Tag','Box');
% hLegend=box_vars;
% xticks(Xtick)
% xticklabels(num2cell(Xtick))
% xtickformat('%1.e')
% plot(N,median(RMSE_dx(:,:,SNR_jj),2),'--g'); hold on
% set(gca,'position',pos)




ax.YAxis.Scale ="log";
ax.XAxis.Scale ="log";
set(gca,'YLim',[5e0 1e2],'YTick',10.^(0:2), ...
        'XLim',[2e5 1e7],'XTick',10.^(5:7))
% set(gca,'XTick',[N])
% grid minor
ylabel(' RMSE $f_r$ ($\%$)')
% xlabel('Dataset Size','Position',[442.2150000000001,-27.40000030040742,0],'VerticalAlignment','top','HorizontalAlignment','center')
% text('Dataset Size','Position',[1e7,10.40000030040742,0],'VerticalAlignment','top','HorizontalAlignment','center')
annotation('textbox',...
    [0.45683596214511,0.033878503558792,0.141481598317561,0.075934580553358],...
    'String',{'Dataset Size'},...
    'Interpreter','latex',...
    'EdgeColor','None','FontSize',24);
exportgraphics(gcf,'figures/duffing_3D_convergence_noise_c.png',"Resolution",600)


% xlabel('Dataset Size','Position',[263.655,-20,-1]);
% cleanfigure('targetResolution',50); 
% matlab2tikz('figures/duffing_3D_convergence_noise_c.tex','standalone',true,'width', '\figurewidth','height', '\figureheight');
% matlab2tikz('figures/duffing_3D_boxplot.tex','standalone',true,'width', '11cm','height', '6cm');

