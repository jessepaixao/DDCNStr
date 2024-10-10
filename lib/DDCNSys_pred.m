function [Z_rec,opt_data,elapsedTime] = DDCNSys_pred(Z_star,u_test,Z0_test,datadriven_solver)

% Set data-driven solver parameters
A1=datadriven_solver.A1;
A2=datadriven_solver.A2;
A3=datadriven_solver.A3;
m=datadriven_solver.m;

% Define dataset
x_star=Z_star(:,1);
dx_star=Z_star(:,2);
fr_star=Z_star(:,3);


% Prepare knnsearch distance algorithm
% X=[x_star dx_star fr_star];
X=[x_star*sqrt(A1/2) dx_star*sqrt(A2/2) fr_star*sqrt(A3/2)];

% Setup nearest neighbor searcher algorithm
% Mdl1=createns(X,'NSMethod','kdtree','Distance','euclidean');
Mdl1=createns(X,'NSMethod','exhaustive','Distance','euclidean');
% Mdl1 = ExhaustiveSearcher(X);
% Mdl1=createns(X,'NSMethod','hnsw','Distance','euclidean');
% Mdl1=createns(X,'hnsw');


% Select initial point "nonlinear state" (x,dx,fn) for solver
% ii=randi(length(x_star));
fr_pred=Z0_test(1);
x_pred=Z0_test(2);
dx_pred=Z0_test(3);


% Define parameters of Newmark integration scheme
Gamma=datadriven_solver.newmark_gamma;
Beta=datadriven_solver.newmark_beta;
dt=1/datadriven_solver.Fs;

% Initialized reconstructed vectors
x_rec=zeros(length(u_test),1);
dx_rec=zeros(length(u_test),1);
ddx_rec=zeros(length(u_test),1);
fr_rec=zeros(length(u_test),1);
x_rec(1,1)=Z0_test(1);
dx_rec(1,1)=Z0_test(2);
fr_rec(1,1)=Z0_test(3);

tol=datadriven_solver.conv_tol;
max_steps=datadriven_solver.max_steps;

opt_data=struct;

N=length(u_test);

% if datadriven_solver.monitoring=='on'
tic
for jj=1:N-1
    % 
    clc
    fprintf('Iteration: %d / %d \n', [jj,N-1])
    % jj

    % Select initial a dataset point (z_star)
    % ii=1;
    ii=randi(length(x_star));
    x_star_i=x_star(ii);
    dx_star_i=dx_star(ii);
    fr_star_i=fr_star(ii);

    f_obj=1;
    conv=1;

    % Vector to store evolution of fn and x along optimization
    fr_vec=zeros(max_steps,1);
    x_vec=zeros(max_steps,1);
    dx_vec=zeros(max_steps,1);
    fr_star_vec=zeros(max_steps,1);
    x_star_vec=zeros(max_steps,1);
    dx_star_vec=zeros(max_steps,1);
    conv_vec=zeros(max_steps,1);
    f_obj_vec=zeros(max_steps,1);

    kk=1;
    fr_vec(kk)=fr_pred;
    x_vec(kk)=x_pred;
    dx_vec(kk)=dx_pred;
    fr_star_vec(kk)=fr_star_i;
    x_star_vec(kk)=x_star_i;
    dx_star_vec(kk)=dx_star_i;
    conv_vec(kk)=conv;
    f_obj_vec(kk)=A1/2*(x_pred-x_star_i)^2+A2/2*(dx_pred-dx_star_i)^2+A3/2*(fr_pred-fr_star_i)^2;

    
    while conv>tol && kk<max_steps

        kk=kk+1;

        % Solve first optimization problem (staggered scheme)
        C1=-Gamma/(Beta*dt)*x_rec(jj)+(1-Gamma/Beta)*dx_rec(jj)+dt*(1-Gamma/(2*Beta))*ddx_rec(jj);
        C2=m/(Beta*dt^2)*x_rec(jj)+m/(Beta*dt)*dx_rec(jj)+(1/(2*Beta)-1)*m*ddx_rec(jj);
        A=[A1+A2*(Gamma/(Beta*dt))^2 -A3*m/(Beta*dt^2);
            m/(Beta*dt^2) 1];
        b=[A1*x_star_i+A2*Gamma/(Beta*dt)*dx_star_i-A3*m/(Beta*dt^2)*fr_star_i-C1*A2*Gamma/(Beta*dt);
           u_test(jj+1)+C2];
        sol=A\b;

        % Store previous solution before updating
        x_pred_old=x_pred;
        dx_pred_old=dx_pred;
        fr_pred_old=fr_pred;
    
        % Update solution of the first optimization problem (z_t+1)
        x_pred=sol(1);
        dx_pred=Gamma/(Beta*dt)*(x_pred-x_rec(jj))+(1-Gamma/Beta)*dx_rec(jj)+dt*(1-Gamma/(2*Beta))*ddx_rec(jj);
        fr_pred=sol(2);
            
        % Solve second optimization problem (staggered scheme); 
        IdxKDT = knnsearch(Mdl1,[x_pred*sqrt(A1/2) dx_pred*sqrt(A2/2) fr_pred*sqrt(A3/2)]);
    

        % Convergence criterion
        % conv=abs(x_pred_old-x_pred)/x_norm_factor+abs(dx_pred_old-dx_pred)/dx_norm_factor+abs(fr_pred_old-fr_pred)/fr_norm_factor;
        f_obj_new=A1/2*(x_pred-x_star_i)^2+A2/2*(dx_pred-dx_star_i)^2+A3/2*(fr_pred-fr_star_i)^2;
        % conv=abs(x_pred_old-x_pred)+abs(dx_pred_old-dx_pred)+abs(fr_pred_old-fr_pred);
        conv=abs(f_obj_new-f_obj);

        % Update "nonlinear state"
        x_star_i=x_star(IdxKDT);
        dx_star_i=dx_star(IdxKDT);
        fr_star_i=fr_star(IdxKDT);

        f_obj=f_obj_new;
               
        
        % Update "nonlinear dynamics solution" 
        fr_vec(kk)=fr_pred;
        x_vec(kk)=x_pred;
        dx_vec(kk)=dx_pred;
        fr_star_vec(kk)=fr_star_i;
        x_star_vec(kk)=x_star_i;
        dx_star_vec(kk)=dx_star_i;
        conv_vec(kk)=conv;
        f_obj_vec(kk)=f_obj;

        

    end

    opt_data(jj+1).steps=kk;
    opt_data(jj+1).fr=fr_vec;
    opt_data(jj+1).x=x_vec;
    opt_data(jj+1).dx=dx_vec;
    opt_data(jj+1).fr_star_i=fr_star_vec;
    opt_data(jj+1).x_star_i=x_star_vec;
    opt_data(jj+1).dx_star_i=dx_star_vec;

    % Update "nonlinear dynamics solution" for the timestep jj+1
    x_rec(jj+1,1)=x_pred;
    dx_rec(jj+1,1)=dx_pred;
    ddx_rec(jj+1,1)=1/(Beta*dt^2)*(x_rec(jj+1)-x_rec(jj))-1/(Beta*dt)*dx_rec(jj)-(1/(2*Beta)-1)*ddx_rec(jj);
    fr_rec(jj+1,1)=fr_pred;

    
end
elapsedTime = toc;

Z_rec=[x_rec dx_rec ddx_rec fr_rec];

end