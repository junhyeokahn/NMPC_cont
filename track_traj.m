clear all; close all; clc;
warning('off','MATLAB:nargchk:deprecated');
%% Constants
n = 4;
m = 1;

%% Setup Geodesic Numerics

load 'metric_FLR.mat';

W = @(x) W_mat(wrapToPi(x(1)));
dW = @(x) {dW_x1_mat(wrapToPi(x(1))), zeros(4), zeros(4), zeros(4)};

sigma_ThBw = 0.0491;
w_lower = 3.5171;

geodesic_N = 2;

[geo_Prob,geo_Ke,geo_we,T_e,T_dot_e,geo_Aeq] = ...
        setup_geodesic_calc(n,geodesic_N,W,dW);

% Assemble geodesic struct for MPC
geodesic_MPC = struct('geo_Prob',geo_Prob,'W',W,'geodesic_N',geodesic_N,'T_e',T_e,'T_dot_e',T_dot_e,...
                      'geo_Aeq',geo_Aeq);

%% Setup NMPC Problem

mass = 1;
l = 1;
I = (1/3)*m*(2*l)^2;
g = 9.81;
sigma = 100;
J = 1;
b = 1;

f  = @(x) [x(2);
           (mass*g*l/I)*sin(x(1)) - (sigma/I)*(x(1)-x(3));
           x(4);
           (sigma/J)*(x(1)-x(3)) - (b/J)*x(4)];
       
B = [zeros(3,1);
     1/J];

df = @(x) [0, 1, 0, 0;
           (mass*g*l/I)*cos(x(1))-(sigma/I), 0, (sigma/I), 0;
           zeros(1,3), 1;
           (sigma/J), 0, -(sigma/J), -(b/J)];

B_w = [0,1/I,0,0;
       0,0,0,1/J]';

% B_w = eye(4);

w_max = 0.15;

M_ccm = eye(4);
lambda =  2.5;
d_bar = (w_max*sigma_ThBw/lambda);

P = eye(4);
alpha = 1e-3;

Tp = 8;
delta = 0.1;
dt = 0.005;

N_mpc = 80;

state_constr_low = -[pi;3.83;pi;3.83];
ctrl_constr_low = -33*ones(m,1);

q_eq = 30*(pi/180);
th_eq = q_eq - (mass*g*l/sigma)*sin(q_eq);
x_eq = [q_eq;0;th_eq;0];
u_eq = -sigma*(q_eq - th_eq);

[NMPC_Prob,L_e,L_e_full] = setup_NMPC(n,m,...
    f,B,df, state_constr_low,ctrl_constr_low,...
    N_mpc,Tp,delta,dt,...
    P,alpha,geodesic_MPC,d_bar^2,...
    x_eq,u_eq);

%% Test MPC Solve

link_ang = -pi + 15*(pi/180);
mot_ang = link_ang + 5*(pi/180);

test_state = [link_ang;
             -30*(pi/180);
             mot_ang;
               0];
tic
[NMPC_state,NMPC_ctrl,converged_MPC] = compute_NMPC(NMPC_Prob,...
    test_state,test_state,state_constr_low,x_eq,...
    n,m,N_mpc,L_e_full);
toc
disp(converged_MPC);

% Visualize
close all
figure(); 
hold on
%RCI set
Ellipse_plot(M_ccm(1:2,1:2)*((1/w_lower)/d_bar^2), NMPC_state(1,1:2)',25,'k');
%Terminal set
Ellipse_plot(P(1:2,1:2)*(1/alpha),x_eq(1:2),25,'r');

%Full MPC traj
plot(NMPC_state(:,1),NMPC_state(:,2),'r-','linewidth',2);
plot(NMPC_state(:,3),NMPC_state(:,4),'k-','linewidth',2);

grid on
axis equal
xlabel('$q, \theta$','interpreter','latex'); ylabel('$\dot{q}, \dot{\theta}$','interpreter','latex');
set(findall(gcf,'type','text'),'FontSize',32);set(gca,'FontSize',32)

% pause;

%% Test Geodesic Numerics

test_rand = randn(4,1); test_rand = test_rand/norm(test_rand);
test_state = NMPC_state(1,:)' + (d_bar*sqrt(w_lower))*test_rand;
% test_state = x_eq; %just to see cool curved lines

tic
[X, X_dot,J_opt,converged_geo] = ...
    compute_geodesic_tom(geo_Prob,n,geodesic_N,...
            NMPC_state(1,:)',test_state,...
            T_e,T_dot_e,geo_Aeq);
toc;
disp(converged_geo);
disp('Geo dist: '); disp(sqrt(J_opt));

figure(1)
plot(X(1,:),X(2,:),'b-','linewidth',2);
plot(X(3,:),X(4,:),'b-','linewidth',2);

% pause;
                
%% Setup Auxiliary controller
aux_Prob = setup_opt_aux(m);

tic
[ctrl_opt,converged_aux] = compute_opt_aux(aux_Prob,geo_Ke,X,X_dot,J_opt,...
                            W,f,B,NMPC_ctrl(1,:)',lambda);
toc;
disp(converged_aux);
disp('opt_control:'); disp(ctrl_opt);

pause;

%% Set up non-linear sim


ode_options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);

t_end = Tp;
solve_t = (0:dt:t_end)';
T_steps = length(solve_t)-1;

t = cell(T_steps,1);
Geod = cell(T_steps,1);
Aux_ctrl = zeros(T_steps,m);

x_act = zeros(T_steps+1,n);
x_act(1,:) = test_state';
state_0 = test_state;

True_ctrl = zeros(T_steps,m);

ctrl_solve_time = zeros(T_steps,3);
ctrl_solve_time(:,1) = NaN;

opt_solved = NaN(T_steps,3);

geo_energy = zeros(T_steps,2);
geo_energy(:,2) = NaN;

      
%% Simulate

for i = 1:T_steps
%     fprintf('%d/%d \n',i,T_steps);
    
    x_nom = NMPC_state(i,:)';
    u_nom = NMPC_ctrl(i,:)';
    
    %Optimal Control
    tic
    [X, X_dot,J_opt,opt_solved(i,2)] = compute_geodesic_tom(geo_Prob,...
    n,geodesic_N,x_nom,state_0,T_e,T_dot_e,geo_Aeq);
    ctrl_solve_time(i,2) = toc;
    
    Geod{i} = X';
    geo_energy(i,1) = J_opt;
    
    tic
    [Aux_ctrl(i,:),opt_solved(i,3)] = compute_opt_aux(aux_Prob,geo_Ke,...
        X,X_dot,J_opt,W,f,B,u_nom,lambda);
    ctrl_solve_time(i,3) = toc;   
    
    True_ctrl(i,:) = u_nom+Aux_ctrl(i,:);

    %Simulate Optimal
%     w_dist = w_max*[sin(2*pi*solve_t(i))*cosd(30);
%                     cos(2*pi*solve_t(i))*sind(30)];

    dist_dir = (X_dot(:,geo_Ke+1))'*(W(state_0)\eye(4))*B_w;
    w_dist = w_max*(dist_dir'/norm(dist_dir));
    
    [d_t,d_state] = ode113(@(t,d_state)ode_sim(t,d_state,True_ctrl(i,:)',...
        f,B,B_w,w_dist),[solve_t(i),solve_t(i+1)],state_0,ode_options);
%     t{i} = d_t;
%     True_state{i} = d_state;
    state_0 = d_state(end,:)';
%     state_0 = state_0 + (f(state_0) + B*True_ctrl(i,:)' + B_w*w_dist)*dt;

    x_act(i+1,:) = state_0';
end

%% Plot
close all

% State Trajectory
figure()
hold on
plot(solve_t, x_act,'linewidth',2);
grid on
xlabel('Time [s]'); ylabel('States'); 
h_leg = legend('$q$','$\dot{q}$','$\theta$','$\dot{\theta}$');
set(h_leg,'interpreter','latex');
set(findall(gcf,'type','text'),'FontSize',32);set(gca,'FontSize',32)

%Control Trajectory
figure()
hold on
plot(solve_t,NMPC_ctrl,'b-','linewidth',2);
plot(solve_t(1:end-1),True_ctrl,'r-','linewidth',2);
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
hold on
plot(NMPC_state(:,1),NMPC_state(:,2),'r--','linewidth',2);
plot(NMPC_state(:,3),NMPC_state(:,4),'k--','linewidth',2);
plot(x_act(:,1),x_act(:,2),'r-','linewidth',2);
plot(x_act(:,3),x_act(:,4),'k-','linewidth',2);

Ellipse_plot(P(1:2,1:2)*(1/(alpha)), x_eq(1:2),30,'k');
xlabel('$q, \theta$','interpreter','latex'); 
ylabel('$\dot{q}, \dot{\theta}$','interpreter','latex');
set(findall(gcf,'type','text'),'FontSize',32);set(gca,'FontSize',32)
grid on; 
axis equal

%Solve Time
figure()
hold on
plot(solve_t(1:end-1),ctrl_solve_time(:,2),'bs','markersize',10,'markerfacecolor','m');
plot(solve_t(1:end-1),ctrl_solve_time(:,3),'rd','markersize',10,'markerfacecolor','k');
grid on
legend('Geo','Aux');
xlabel('Time [s]');
ylabel('Solve time [s]'); title('Solve time');
set(findall(gcf,'type','text'),'FontSize',32);set(gca,'FontSize',32)

%Solve success
figure()
hold on
plot(solve_t(1:end-1),opt_solved(:,2),'bs','markersize',15,'markerfacecolor','m');
plot(solve_t(1:end-1),opt_solved(:,3),'rd','markersize',10,'markerfacecolor','k');
grid on
legend('Geodesic (0,1,6)','Aux (0)');
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

save('FLR_sim.mat');