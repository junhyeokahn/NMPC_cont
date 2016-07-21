clear all; close all; clc;
warning('off','MATLAB:nargchk:deprecated');
%% Constants
n = 2;
m = 1;

%% Setup Geodesic Numerics

% W = @(x) [4.258279173109496,  -0.934234854771844;
%         -0.934234854771844,   3.766923169705589];
% dW = @(x){zeros(2), zeros(2)};
% 
% geodesic_N = 4;
% 
% [geo_Prob,geo_Ke,geo_we,T_e,T_dot_e,geo_Aeq] = ...
%         setup_geodesic_calc(n,geodesic_N,W,dW);

%Assemble geodesic struct for MPC
% geodesic_MPC = struct('geo_Prob',geo_Prob,'geodesic_N',geodesic_N,'T_e',T_e,'T_dot_e',T_dot_e,...
%                       'geo_Aeq',geo_Aeq);

%% Setup NMPC Problem

f  = @(x) [-1*x(1) + 2*x(2);
           -3*x(1) + 4*x(2) - 0.25*(x(2)^3)];
B = [0.5;-2];

df = @(x) [-1, 2;
           -3, 4-0.75*(x(2)^2)];
B_w = [0;1];

w_max = 0.1;

% M_ccm = W(0)\eye(2);
lambda =  1.742857142857143;
% d_bar = (w_max*sqrt(max(eig(M_ccm)))/lambda);

M_ccm = diag([39.0251, 486.0402]);
d_bar = 1;

N_mpc = 50;

P = [7.9997, -12.2019;
    -12.2019, 27.0777];
alpha = 10;

Tp = 1.5;
delta = 0.1;
dt = 0.005;

state_constr_low = [-4.94;-4.94]; 
ctrl_constr_low = -1.793*ones(m,1);
% ctrl_constr_low = -1.15*ones(m,1);

x_eq = [0;0];

[NMPC_Prob,L_e,L_e_full] = setup_NMPC(n,m,...
    f,B,df, state_constr_low,ctrl_constr_low,...
    N_mpc,Tp,delta,dt,...
    P,alpha,M_ccm,d_bar^2,...
    x_eq);

%% Test MPC Solve

test_state = [3.4;-2.4];
tic
[NMPC_state,NMPC_ctrl,converged_MPC] = compute_NMPC(NMPC_Prob,...
    test_state,test_state,state_constr_low,x_eq,...
    n,m,N_mpc,L_e_full);
toc
disp(converged_MPC);

% pause;
% Visualize
close all
figure(); 
hold on
%RCI set
Ellipse_plot(M_ccm*(1/d_bar^2), NMPC_state(1,:)',25,'k');
%Terminal set
Ellipse_plot(P*(1/alpha),NMPC_state(end,:)',25,'r');
%Implemnted segment of MPC traj
plot(NMPC_state(1:(delta/dt)+1,1),...
     NMPC_state(1:(delta/dt)+1,2),'r-','linewidth',2);
%Full MPC traj
plot(NMPC_state(:,1),NMPC_state(:,2),'b-','linewidth',1);
grid on
axis equal
xlabel('x_1'); ylabel('x_2');

pause;

%% Test Geodesic Numerics

% tic
% [X, X_dot,J_opt,converged_geo] = ...
%     compute_geodesic_tom(geo_Prob,n,geodesic_N,...
%             NMPC_state(1,:)',test_state,...
%             T_e,T_dot_e,geo_Aeq);
% toc;
% disp(converged_geo);

                
%% Setup Auxiliary controller
% aux_Prob = setup_opt_aux(m);
% 
% tic
% [ctrl_opt,converged_aux] = compute_opt_aux(aux_Prob,geo_Ke,X,X_dot,J_opt,...
%                             W,f,B,NMPC_ctrl(1,:)',lambda);
% toc;
% disp(converged_aux);
% pause;
%% Set up non-linear sim

ode_options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);

t_end = 10;
solve_t = (0:dt:t_end)';
T_steps = length(solve_t)-1;

t = cell(T_steps,1);
Geod = cell(T_steps,1);
Aux_ctrl = zeros(T_steps,m);

True_state = cell(T_steps,1);
True_ctrl = zeros(T_steps,m);

dt_MPC = delta;
solve_MPC = (0:dt_MPC:t_end)';
T_steps_MPC = length(solve_MPC)-1;

state_0 = test_state;
state_0_MPC = state_0;

MPC_state = cell(T_steps_MPC,1);
MPC_ctrl = cell(T_steps_MPC,1);

ctrl_solve_time = zeros(T_steps,3);
ctrl_solve_time(:,1) = NaN;

opt_solved = zeros(T_steps,3);
opt_solved(:,1) = NaN;

geo_energy = zeros(T_steps,2);
geo_energy(:,2) = NaN;

i_mpc = 0;

      
%% Simulate

for i = 1:T_steps
    fprintf('%d/%d \n',i,T_steps);
    %First Solve MPC
    if (mod(solve_t(i),delta)==0)
        geo_energy(i,1) = (state_0-state_0_MPC)'*M_ccm*(state_0-state_0_MPC);
%         [~, ~,J_opt,~] = compute_geodesic_tom(geo_Prob,...
%             n,geodesic_N,state_0_MPC,state_0,T_e,T_dot_e,geo_Aeq);
%         geo_dist(i,1) = J_opt;
        
        tic
        [NMPC_state,NMPC_ctrl,opt_solved(i,1)] = ...
         compute_NMPC(NMPC_Prob,state_0,state_0_MPC,...
                    state_constr_low,x_eq,n,m,N_mpc,L_e);
        ctrl_solve_time(i,1) = toc;
        
        i_mpc = i_mpc + 1;
        
        MPC_state{i_mpc} = NMPC_state;
        MPC_ctrl{i_mpc} = NMPC_ctrl;
        
        x_nom = MPC_state{i_mpc}(1,:);
        u_nom = MPC_ctrl{i_mpc}(1,:);
        
%         [~, ~,J_opt,~] = compute_geodesic_tom(geo_Prob,...
%             n,geodesic_N,x_nom',state_0,T_e,T_dot_e,geo_Aeq);
        
        geo_energy(i,2) = (state_0-x_nom')'*M_ccm*(state_0-x_nom');
%         geo_dist(i,2) = J_opt;
        
        %update starting state for next MPC problem
        state_0_MPC = MPC_state{i_mpc}(end,:)';
    else
        i_mpc_use = round((mod(solve_t(i),delta))/dt)+1;
        x_nom = MPC_state{i_mpc}(i_mpc_use,:);
        u_nom = MPC_ctrl{i_mpc}(i_mpc_use,:);
    end
    
    %Optimal Control
    tic
%     [X, X_dot,J_opt,opt_solved(i,2)] = compute_geodesic_tom(geo_Prob,...
%     n,geodesic_N,x_nom',state_0,T_e,T_dot_e,geo_Aeq);
%     ctrl_solve_time(i,2) = toc;
%     
%     Geod{i} = X';
    geo_energy(i,1) = (state_0-x_nom')'*M_ccm*(state_0-x_nom');
    
%     tic
%     [Aux_ctrl(i,:),opt_solved(i,3)] = compute_opt_aux(aux_Prob,geo_Ke,...
%         X,X_dot,J_opt,W,f,B,u_nom,lambda);
%     ctrl_solve_time(i,3) = toc;   
    Aux_ctrl(i,:) = ([-1.3696, 5.1273]*(state_0-x_nom'))';
    
    True_ctrl(i,:) = u_nom+Aux_ctrl(i,:);

    %Simulate Optimal
%     w_dist = 0.0*w_max + (w_max)*sin(2*pi*0.5*solve_t(i));
    w_dist = w_max;
    [d_t,d_state] = ode113(@(t,d_state)ode_sim(t,d_state,True_ctrl(i,:)',...
        f,B,B_w,w_dist),[solve_t(i),solve_t(i+1)],state_0,ode_options);
    t{i} = d_t;
    True_state{i} = d_state;
    state_0 = d_state(end,:)';    
end

%% Plot
close all

% State Trajectory
figure()
hold on
for i = 1:T_steps
    plot(t{i},True_state{i}(:,1),'r-','linewidth',2);
    plot(t{i},True_state{i}(:,2),'b-','linewidth',2);
end
grid on
xlabel('Time [s]'); ylabel('States'); legend('x_1','x_2');
set(findall(gcf,'type','text'),'FontSize',32);set(gca,'FontSize',32)

%Control Trajectory
figure()
for i_mpc = 1:T_steps_MPC
    stairs([(i_mpc-1)*delta:dt:i_mpc*delta]',...
             MPC_ctrl{i_mpc},'b-','linewidth',2); hold on
end
stairs(solve_t(1:end-1),True_ctrl,'r-','linewidth',2);
xlabel('Time [s]');
ylabel('u(t)'); 
set(findall(gcf,'type','text'),'FontSize',32);set(gca,'FontSize',32)
grid on

figure()
plot(solve_t(1:end-1),Aux_ctrl,'b-','linewidth',2); 
xlabel('Time [s]');
ylabel('Ancillary control k(x^{*},x)');
set(findall(gcf,'type','text'),'FontSize',32);set(gca,'FontSize',32)
grid on

%2D State Plot
figure()
for i_mpc = 1:T_steps_MPC
    plot(MPC_state{i_mpc}(:,1),...
          MPC_state{i_mpc}(:,2),'r-','linewidth',2);
    hold on; 
    Ellipse_plot(M_ccm*(1/(d_bar^2)), MPC_state{i_mpc}(1,:)',30,'g');
%     Ellipse_plot(M_ccm*(1/(d_bar^2)), MPC_state{i_mpc}(end,:)',30,'r');
end
for i = 1:T_steps
    plot(True_state{i}(:,1),...
          True_state{i}(:,2),'b-','linewidth',2);
%     plot(Geod{i}(:,1),Geod{i}(:,2),'k-','linewidth',2);
end
Ellipse_plot(P*(1/(alpha)), x_eq,30,'k');
xlabel('x_1'); ylabel('x_2');
set(findall(gcf,'type','text'),'FontSize',32);set(gca,'FontSize',32)
grid on; 

figure()
% plot(solve_t(1:end-1),ctrl_solve_time(:,1),'ro','markersize',10,'markerfacecolor','g');
hold on
% plot(solve_t(1:end-1),ctrl_solve_time(:,2),'bs','markersize',10,'markerfacecolor','m');
plot(solve_t(1:end-1),ctrl_solve_time(:,3),'rd','markersize',10,'markerfacecolor','k');
grid on
% legend('Aux');
xlabel('Time [s]');
ylabel('Solve time [s]'); %title('Solve time');
xlabel('Time [s]'); 
set(findall(gcf,'type','text'),'FontSize',32);set(gca,'FontSize',32)

figure()
plot(solve_t(1:end-1),opt_solved(:,1),'ro','markersize',20,'markerfacecolor','g');
hold on
plot(solve_t(1:end-1),opt_solved(:,2),'bs','markersize',15,'markerfacecolor','m');
plot(solve_t(1:end-1),opt_solved(:,3),'rd','markersize',10,'markerfacecolor','k');
grid on
legend('NMPC (1,2,3)','Geodesic (0,1,6)','Aux (0)');
xlabel('Time [s]');
title('Convergence');
set(findall(gcf,'type','text'),'FontSize',32);set(gca,'FontSize',32)
 
figure()
plot(solve_t(1:end-1),sqrt(geo_energy(:,1)),'b-','linewidth',2);
hold on
plot(solve_t(1:end-1),d_bar*ones(T_steps,1),'r-','linewidth',2);
plot(solve_t(1:end-1),sqrt(geo_energy(:,2)),'bo','markersize',10,'markerfacecolor','b','linewidth',2);
grid on
legend('d(x^{*},x)','RCI bound');
xlabel('Time [s]');
set(findall(gcf,'type','text'),'FontSize',32);set(gca,'FontSize',32)

%%

save('NMPC_allgower_alg_run.mat');



