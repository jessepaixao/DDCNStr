%% Data-driven Model-free Computing for Nonlinear Dynamics Systems 
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
close all
rng(123,'twister')

% Aqcquisition parameters
Fs=5000;                     % Sampling frequency

% % 2) Continuous Sweept-sine Exciation
% Excitation to generate training data
RMSu = 50;                        % RMS Amplitude
tspan_dataset  = [0:1/Fs:100]';      % Time vector
N_dataset = length(tspan_dataset);    % Number of samples
T_exc = tspan_dataset(end);         % Period of swept
f_ini = 1;                        % Initial frequency
f_end = 20;                       % Final frequency
u_exc=sin(2*pi*(f_ini*tspan_dataset+(f_end-f_ini)/(2*T_exc)*tspan_dataset.^2));
% u_exc=sin(2*pi*f_ini*tspan_train);
u_exc_dataset = u_exc/rms(u_exc(:,1))*RMSu;

% Excitation to generate test data
RMSu =50;                   % RMS Amplitude
% RMSu =60;                   % RMS Amplitude
tspan_test  = [0:1/Fs:5]';   % Time vector
N_test = length(tspan_test); % Number of sample
T_exc = tspan_test(end);     % Period of swept
f_ini = 5;                  % Initial frequency
f_end = 20;                  % Final frequency
f_exc = 12; 
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


fnl_dataset=u_exc_dataset-m*ddx_dataset-c*dx_dataset-k1*x_dataset;
fnl_test=u_exc_test-m*ddx_test-c*dx_test-k1*x_test;

% Noise
% SNR=50;                     % Signal-to-noise ratio
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

%% FIGURE 2

linewidth=1;

% FIGURE 2a
figure(1)
set(gcf,'units','normalized','outerposition',[0 0 1 0.9])
subplot(221)
plot(tspan_dataset,x_dataset,'r-','LineWidth',linewidth); hold on
ylabel('$x (m)$')
xlim([0 tspan_dataset(end)])
legend('Dataset','Position',[0.366570156643023 0.937442980691757 0.0975027490702616 0.0560465126037598])
subplot(222)
plot(tspan_test,x_test,'k-','LineWidth',linewidth);hold on
xlim([0 tspan_test(end)])
legend('Ref. Sol. (Test)','Position',[0.723825611231659 0.935506597835713 0.180320901211383 0.0560465126037594])
subplot(223)
plot(tspan_dataset,fnl_dataset,'r-','LineWidth',linewidth); hold on
xlim([0 tspan_dataset(end)])
ylabel('$f_{nl}$ (N)')
xlabel('t(s)')
subplot(224)
plot(tspan_test,fnl_test,'k-','LineWidth',linewidth);hold on
xlim([0 tspan_test(end)])
xlabel('t(s)')

% Save figure
exportgraphics(gcf,'figures/duffing_2D_time_signals.pdf','ContentType','vector',"Resolution",600)

% FIGURE 2b
figure(2)
set(gcf,'units','normalized','outerposition',[0 0 1 0.6])
plot(x_dataset,fnl_dataset,'ro','MarkerFaceColor','r','MarkerSize',6); hold on
plot(x_test,fnl_test,'k.','MarkerSize',10)
legend('Dataset','Ref. Sol.(Test)','Location','southeast')
xlabel('$x$ (m)')
ylabel('$f_{nl}$ (N)')

% Save figure
exportgraphics(gcf,'figures/duffing_2D_dataset.pdf','ContentType','vector',"Resolution",600)

%% DATA-DRIVEN SOLVER (CONSTANT HYPERPARAMETER)

Z_star=[x_dataset fnl_dataset];

u_test=u_exc_test;
Z0_test=[x_test(1) fnl_test(1)];

% Define parameters of Newmark
datadriven_solver.newmark_gamma=0.5;
datadriven_solver.newmark_beta=0.25;

datadriven_solver.Fs=Fs; % Sampling frequency of dataset

% Define constants of objective function
datadriven_solver.A1=2/max(x_dataset)^2;
datadriven_solver.A2=2/max(fnl_dataset)^2;
% datadriven_solver.A1=1;
% datadriven_solver.A2=1;
datadriven_solver.m=m;
datadriven_solver.c=c;
datadriven_solver.k1=k1;

datadriven_solver.conv_tol=1e-10;
datadriven_solver.max_steps=100;

% Run the data-driven solver
[Z_rec,opt_data,elapsedTime] = DDCNSys_2D_pred(Z_star,u_test,Z0_test,datadriven_solver);

x_rec=Z_rec(:,1);
dx_rec=Z_rec(:,2);
ddx_rec=Z_rec(:,3);
fn_rec=Z_rec(:,4);

RMSE_x=rms(x_rec-x_test)/rms(x_test)*100;
RMSE_fn=rms(fn_rec-fnl_test)/rms(fnl_test)*100;

% fileName = ['results/' datestr(now, 'dd-mmm-yyyy_HHMMSS') '_introductory_example_2D']
% save(fileName)

%% PLOT ANIMATION OF DATA DRIVEN SOLVER PROGRESS

close all
jj=10059;


% Set data-driven solver parameters
A1=datadriven_solver.A1;
A2=datadriven_solver.A2;
% A1=k1;
% A2=1;

% Define constants
dt=1/Fs;
g1=m/(Beta*dt^2)+c*Gamma/(Beta*dt)+k1;
g2=m/(Beta*dt^2)+c*Gamma/(Beta*dt);
g3=m/(Beta*dt)+(Gamma/Beta-1)*c;
g4=(1/(2*Beta)-1)*m+dt*(Gamma/(2*Beta)-1)*c;


fn_elp_min=min(fnl_dataset);
fn_elp_max=max(fnl_dataset);
x_cntr_min=(u_test(jj+1)+g2*x_rec(jj)+g3*dx_rec(jj)+g4*ddx_rec(jj)-fn_elp_min)/g1;
x_cntr_max=(u_test(jj+1)+g2*x_rec(jj)+g3*dx_rec(jj)+g4*ddx_rec(jj)-fn_elp_max)/g1;
x_cntr=linspace(x_cntr_min,x_cntr_max,5000);
fn_cntr=u_test(jj+1)+g2*x_rec(jj)+g3*dx_rec(jj)+g4*ddx_rec(jj)-g1*x_cntr;


linewidth=1;

% FIGURE 3b
fig=figure(1)
set(gcf,'units','normalized','outerposition',[0 0 1 0.8])
% Subfigure
tiledlayout(1,3);
nexttile(1,[1 2]);
plot(x_dataset,fnl_dataset,'ro','MarkerFaceColor','r','MarkerSize',4); hold on
plot(x_cntr,fn_cntr,'b-','LineWidth',1.5)
for kk=1:opt_data(jj+1).steps-1
    x_star_i=opt_data(jj+1).x_star_i(kk);
    fn_star_i=opt_data(jj+1).fn_star_i(kk);

    Lc2=(1/((g1^2/A1+1/A2)))*(u_test(jj+1)+g2*x_rec(jj)+g3*dx_rec(jj)+g4*ddx_rec(jj)-g1*x_star_i-fn_star_i)^2;
    a2=Lc2/A1;
    b2=Lc2/A2;

    % Parametric form of ellipse equation
    tt=linspace(0,2*pi,5000);
    x_elp=x_star_i+sqrt(a2)*cos(tt);
    fn_elp=fn_star_i+sqrt(b2)*sin(tt);


    plot(opt_data(jj+1).x_star_i(kk),opt_data(jj+1).fn_star_i(kk),'kx','LineWidth',linewidth);
    plot(x_elp,fn_elp,'k-','LineWidth',0.1)
    plot(opt_data(jj+1).x(kk+1),opt_data(jj+1).fn(kk+1),'mo','MarkerFaceColor','m','MarkerSize',4);
    plot(x_test(jj+1),fnl_test(jj+1),'o','Linewidth',1,'MarkerSize',6,'MarkerFaceColor','y','MarkerEdgeColor','k'); hold on
end
legend('Dataset','Constraint','$z^{*(k)}_j$','Ellipse','$z^{(k)}_{i}$','Ref. Sol.($t_i^{(k)}$)','Position',[0.0927602729875811 0.923259311739601 0.583537813834825 0.0764357134534849],...
    'Orientation','horizontal')
ylabel('$f_{nl}$ (N)')
xlabel('$x$ (m)')


% Subfigure
% subplot(1,3,[3])
nexttile;
plot(x_dataset,fnl_dataset,'ro','MarkerFaceColor','r','MarkerSize',4); hold on
plot(x_cntr,fn_cntr,'b-','LineWidth',1.5)
for kk=1:opt_data(jj+1).steps-1

    x_star_i=opt_data(jj+1).x_star_i(kk);
    fn_star_i=opt_data(jj+1).fn_star_i(kk);

    Lc2=(1/((g1^2/A1+1/A2)))*(u_test(jj+1)+g2*x_rec(jj)+g3*dx_rec(jj)+g4*ddx_rec(jj)-g1*x_star_i-fn_star_i)^2;
    a2=Lc2/A1;
    b2=Lc2/A2;

    % Parametric form of ellipse equation
    tt=linspace(0,2*pi,5000);
    x_elp=x_star_i+sqrt(a2)*cos(tt);
    fn_elp=fn_star_i+sqrt(b2)*sin(tt);

    plot(opt_data(jj+1).x_star_i(kk),opt_data(jj+1).fn_star_i(kk),'kx','LineWidth',linewidth);
    plot(x_elp,fn_elp,'k-','LineWidth',0.1)
    plot(opt_data(jj+1).x(kk+1),opt_data(jj+1).fn(kk+1),'mo','MarkerFaceColor','m','MarkerSize',4);
    plot(x_test(jj+1),fnl_test(jj+1),'o','Linewidth',1,'MarkerSize',6,'MarkerFaceColor','y','MarkerEdgeColor','k'); hold on
end

% xlabel('$x$ (m)')
title('Zoom')
xlim([0.11 0.13])
ylim([-100 300])

% Save figure
exportgraphics(gcf,'figures/duffing_2D_ellipses_t3.pdf','ContentType','vector',"Resolution",600)
%% PLOT TESTE DATASET X RECONSTRUCTED 
% close all

jj=length(x_rec)-1;


% FIGURE 3a
figure(1)
set(gcf,'units','normalized','outerposition',[0 0 1 1])
% Subplot 1
% subplot(311)
tiledlayout(3,1);
nexttile;
plot(tspan_test(1:jj+1),x_test(1:jj+1),'k-','MarkerSize',10,'Linewidth',1.5,'MarkerSize',12); hold on
plot(tspan_test(1:jj+1),x_rec(1:jj+1),'b--','MarkerSize',10,'Linewidth',1,'MarkerSize',12); hold on
plot(tspan_test(1:jj+1),(x_rec(1:jj+1)-x_test(1:jj+1)),'g-','Linewidth',1,'MarkerSize',12); hold on
ylabel('$x$(m)')
xlim([0 tspan_test(end)])
legend('Ref. Sol.(Test)','Data-driven Pred.','Abs. Error','Position',[0.814874414987007 0.860141179740593 0.176713387326348 0.134762401526625])
% Subplot 2
% subplot(312)
nexttile;
plot(tspan_test(1:jj+1),fnl_test(1:jj+1),'k-','MarkerSize',10,'Linewidth',1.5,'MarkerSize',12); hold on
plot(tspan_test(1:jj+1),fn_rec(1:jj+1),'b--','MarkerSize',10,'Linewidth',1,'MarkerSize',12); hold on
plot(tspan_test(1:jj+1),(fn_rec(1:jj+1)-fnl_test(1:jj+1)),'g-','Linewidth',1,'MarkerSize',12); hold on
xlabel('t(s)')
ylabel('$f_{nl}$ (N)')
xlim([0 tspan_test(end)])
% Subplot 3
% subplot(313)
nexttile;
bar(tspan_test(2:jj+1),vertcat(opt_data(:).steps),'b','EdgeColor','none'); hold on
xlabel('t(s)')
ylabel('$\#Iterations$')
xlim([0 tspan_test(end)])

% Save figure
exportgraphics(gcf,'figures/duffing_2D_pred_data_driven.pdf','ContentType','vector',"Resolution",600)
