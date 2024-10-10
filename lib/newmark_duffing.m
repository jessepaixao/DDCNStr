%% IMPLEMENTATION OF NEWMARK METHOD FOR SOLVING FORCED DUFFING OSCILLATOR
%  Description:the newmark implementation folllwed the algorithm presented
%  by Chopra 2015
%  Versions:
%   1.0 - 22/04/24
%   1.1 - 14/05/24
%  Author: Jessé Paixao
% 
function [x_i,dx_i,ddx_i] = newmark_duffing(m,c,k1,k2,k3,u_exc,t,xc0,Gamma,Beta,tol)

dt=t(2)-t(1);

% Initial Calculations
x0=xc0(1);
dx0=xc0(2);
kT0=k1;
fs0=k1*x0+k3*x0^3;
ddx0=(u_exc(1)-c*dx0-fs0)/m;

a1=1/(Beta*dt^2)*m+Gamma/(Beta*dt)*c;
a2=1/(Beta*dt)*m+(Gamma/Beta-1)*c;
a3=(1/(2*Beta)-1)*m+dt*(Gamma/(2*Beta)-1)*c;

x_i=zeros(length(u_exc),1);
x_i(1)=x0;
dx_i=zeros(length(u_exc),1);
dx_i(1)=dx0;
ddx_i=zeros(length(u_exc),1);
ddx_i(1)=ddx0;
fs_i=zeros(length(u_exc),1);
fs_i(1)=fs0;
kT_i=zeros(length(u_exc),1);
kT_i(1)=kT0+a1;
u_hat_i=zeros(length(u_exc),1);

% Calculations for each time instant
for ii=1:length(u_exc)-1

    x_i(ii+1)=x_i(ii);
    fs_i(ii+1)=fs_i(ii);
    kT_i(ii+1)=kT_i(ii);

    u_hat_i(ii+1)=u_exc(ii+1)+a1*x_i(ii)+a2*dx_i(ii)+a3*ddx_i(ii);

    R_hat=u_hat_i(ii+1)-fs_i(ii+1)-a1*x_i(ii+1);

    % Calculations for each iteration (Newton-Raphson)
    while abs(R_hat)>tol

        kT_hat_j=kT_i(ii+1)+a1;
        Dx=R_hat/kT_hat_j;
        x_ij=x_i(ii+1);
        x_i(ii+1)=x_i(ii+1)+Dx;

        % fs_i(ii+1)=fs_i(ii+1)+k1*(x_i(ii+1)-x_ij)+k2*(x_i(ii+1)^2-x_ij^2)+k3*(x_i(ii+1)^3-x_ij^3);
        % kT_i(ii+1)=k1+k2*2*x_i(ii+1)+k3*3*x_i(ii+1)^2;

        fs_i(ii+1)=fs_i(ii+1)+k1*(x_i(ii+1)-x_ij)+k3*(x_i(ii+1)^3-x_ij^3);
        kT_i(ii+1)=k1+k3*3*x_i(ii+1)^2;

        R_hat=u_hat_i(ii+1)-fs_i(ii+1)-a1*x_i(ii+1);


    end

    % Calculations for velocity and acceleration
    dx_i(ii+1)=Gamma/(Beta*dt)*(x_i(ii+1)-x_i(ii))+(1-Gamma/Beta)*dx_i(ii)+dt*(1-Gamma/(2*Beta))*ddx_i(ii);
    ddx_i(ii+1)=1/(Beta*dt^2)*(x_i(ii+1)-x_i(ii))-(1/(Beta*dt))*dx_i(ii)-(1/(2*Beta)-1)*ddx_i(ii);

end

end
