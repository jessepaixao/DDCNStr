%% Data-driven Model-free Computing in Nonlinear Structural Dynamics (DMCNStr) 
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
set(groot,'defaultAxesFontSize',20)
set(groot,'defaultTextInterpreter','latex')
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');

%% DEFINE ACQUISITION PARAMETERS AND EXCITATION SIGNAL
close all
% rng(123,'twister')

% Aqcquisition parameters
Fs=5000;                     % Sampling frequency

% % 2) Continuous Sweept-sine Exciation
% Excitation to generate training data
% RMSu = 300;                        % RMS Amplitude
RMSu = 50;                        % RMS Amplitude
tspan_dataset  = [0:1/Fs:1000]';      % Time vector
N_dataset = length(tspan_dataset);    % Number of samples
T_exc = tspan_dataset(end);         % Period of swept
f_ini = 1;                        % Initial frequency
f_end = 20;                       % Final frequency
u_exc=sin(2*pi*(f_ini*tspan_dataset+(f_end-f_ini)/(2*T_exc)*tspan_dataset.^2));
% u_exc=normrnd(0,1,[length(tspan_dataset),1]);
% u_exc=rand(length(tspan_dataset),1)-0.5;
% u_exc=lowpass(u_exc,100,Fs);
% u_exc=sin(2*pi*f_ini*tspan_train);
u_exc_dataset = u_exc/rms(u_exc(:,1))*RMSu;

% Excitation to generate test data
RMSu =60;                   % RMS Amplitude
% RMSu =50;                   % RMS Amplitude
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



% Noise
% SNR=20;                     % Signal-to-noise ratio
% sig=x_dataset;
% noisy = randn(size(sig,1),1).*std(sig)/db2mag(SNR);
% x_dataset=x_dataset+noisy;
% 
% sig=dx_dataset;
% noisy = randn(size(sig,1),1).*std(sig)/db2mag(SNR);
% dx_dataset=dx_dataset+noisy;
% 
% sig=ddx_dataset;
% noisy = randn(size(sig,1),1).*std(sig)/db2mag(SNR);
% ddx_dataset=ddx_dataset+noisy;


fr_dataset=u_exc_dataset-m*ddx_dataset;
fr_test=u_exc_test-m*ddx_test;


% PLOT SYSTEM INPUT AND OUTPUT TIME SIGNALS
% FIGURE 5a
figure()
subplot()
set(gcf,'units','normalized','outerposition',[0 0 1 1])
subplot(321)
plot(tspan_dataset,x_dataset,'r-','LineWidth',1.5); hold on
% plot(tspan_test,x_test,'k-');hold on
ylabel('$x (m)$')
xlim([0 tspan_dataset(end)])
legend('Dataset','Location',[0.379900450088758 0.928451135389521 0.0841724556245262 0.0404958684139015])
subplot(322)
% plot(tspan_dataset,x_dataset,'r-'); hold on
plot(tspan_test,x_test,'k-','LineWidth',1.5);hold on
% ylabel('$x (m)$')
xlim([0 tspan_test(end)])
legend(['Ref. Sol. (Test)'],'Location',[0.767102457995406 0.928451139911182 0.137033539200528 0.0404958684139015])
subplot(323)
plot(tspan_dataset,dx_dataset,'r-','LineWidth',1.5); hold on
% plot(tspan_test,dx_test,'k-');hold on
xlim([0 tspan_dataset(end)])
ylabel('$\dot{x}$ (m/s)')
subplot(324)
% plot(tspan_dataset,dx_dataset,'r-'); hold on
plot(tspan_test,dx_test,'k-','LineWidth',1.5);hold on
xlim([0 tspan_test(end)])
% ylabel('$\dot{x}$ (m/s)')
subplot(325)
plot(tspan_dataset,fr_dataset,'r-','LineWidth',1.5); hold on
% plot(tspan_test,ddx_test,'k-');hold on
xlim([0 tspan_dataset(end)])
ylabel('$f_r$ (N)')
xlabel('t(s)')
subplot(326)
% plot(tspan_dataset,ddx_dataset,'r-'); hold on
plot(tspan_test,fr_test,'k-','LineWidth',1.5);hold on
xlim([0 tspan_test(end)])
% ylabel('$\ddot{x} (m/s^2)$')
xlabel('t(s)')

% Save figure
%exportgraphics(gcf,'figures/duffing_3D_time_signals.pdf','ContentType','vector')

% FIGURE 5b
figure(2)
set(gcf,'units','normalized','outerposition',[0 0 1 0.8])
plot3(x_dataset,dx_dataset,fr_dataset,'ro','MarkerFaceColor','r','MarkerSize',1); hold on
plot3(x_test,dx_test,fr_test,'k.','MarkerSize',10)
xlabel('$x$(m)')
ylabel('$\dot{x}$(m/s)')
zlabel('$f_r$ (N)')
[h,icons]=legend('Dataset','Ref. Sol. (Test)','Position',[0.780853309541509 0.825465884143214 0.162154095464538 0.151928194913458])
icons = findobj(icons,'Type','line');
set(icons(1:2),'MarkerSize',4);

% Save figure
exportgraphics(gcf,'figures/duffing_3D_dataset.pdf','ContentType','image','Resolution',600)

%% PLOT DENSITY OF POINTS

% FIGURE 7
figure(1)
set(gcf,'units','normalized','outerposition',[0 0 1 1])
h=scatterhist(x_dataset,dx_dataset,'MarkerSize',0.5,'Color','r','NBins',[200,150]); hold on
xlabel('$x$ (m)','FontSize',24)
ylabel('$\dot{x}$ (m/s)','FontSize',24)

binscatter(x_dataset,dx_dataset,[200 150],'FaceAlpha',0.8); hold on
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
colormap(gca,flip(cmp))
ylabel(cb,'Bin Counts','FontSize',24,'Rotation',270,'Interpreter','Latex')

%Save figure
exportgraphics(gcf,'figures/duffing_3D_density_plot.pdf','ContentType','image','Resolution',600)

%% DATA-DRIVEN SOLVER (CONSTANT HYPERPARAMETER)

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

% Run the data-driven solver
[Z_rec,opt_data,elapsedTime] = DDCNSys_pred(Z_star,u_test,Z0_test,datadriven_solver);

x_rec=Z_rec(:,1);
dx_rec=Z_rec(:,2);
ddx_rec=Z_rec(:,3);
fr_rec=Z_rec(:,4);

RMSE_x=rms(x_rec-x_test)/rms(x_test)*100;
RMSE_dx=rms(dx_rec-dx_test)/rms(dx_test)*100;
RMSE_fr=rms(fr_rec-fr_test)/rms(fr_test)*100;

% save(strcat('results\duffing_3D_conv_N',num2str(N_dataset),'.mat'))
% save(strcat('results\case_2_swept_sine_25_T2_A200_',num2str(N_dataset),'.mat'))
% fileName = ['results/' datestr(now, 'dd-mmm-yyyy_HHMMSS') '_duffing_oscillator_3D']
% save(fileName)


%% PLOT ANIMATION OF DATA DRIVEN SOLVER PROGRESS


Z_test=[x_test dx_test fr_test];

% Define video parameters
animation_config.file_name='animation/case_2_sine_25_T2_A200_2500001';
animation_config.Quality   = 50;  
animation_config.FrameRate = 5;
animation_config.Frame_steps=100;

% 
% myVideo           = VideoWriter('animation/test2','MPEG-4');
% myVideo.Quality   = 100;  
% myVideo.FrameRate = 1;
% myVideo.Frame_steps=50;
close all
plot_animation_driven_solver_formulation4(Z_rec,Z_star,Z_test,opt_data,tspan_test,u_test,datadriven_solver,animation_config)



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
legend('Ref. Sol. (Test)','Data-driven Pred.','Abs. Error','Position',[0.855581709990948 0.879909955777893 0.122161016963136 0.106869838878143])
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
ylabel('$f_r$ (N)')
xlim([0 tspan_test(end)])
subplot(414)
bar(tspan_test(2:jj+1),vertcat(opt_data(:).steps),'b'); hold on
% ylim([0,20])
xlabel('t(s)')
ylabel('$\#Iterations$')
xlim([0 tspan_test(end)])

% Save figure
exportgraphics(gcf,'figures/duffing_3D_pred_time_data_driven.pdf','ContentType','vector')

%% Plot restoring force
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
[h,icons]=legend('Dataset','Ref. Sol. (Test)','Data-driven Pred.','Position',[0.780853309541509 0.825465884143214 0.162154095464538 0.151928194913458])
icons = findobj(icons,'Type','line');
set(icons(1:2),'MarkerSize',4);

% Save figure
exportgraphics(gcf,'figures/duffing_3D_pred_rfs1_data_driven.pdf','ContentType','image','Resolution',600)
