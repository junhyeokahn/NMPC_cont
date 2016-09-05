clear all; close all; clc;
warning('off','MATLAB:nargchk:deprecated');
%% Constants
n = 6;
m = 2;

%% Obstacle info

obs_loc = [-2;-2];
obs_rad = 0.2;
obs = struct('loc',obs_loc,'r',obs_rad);

%% Setup Metric

load 'metric_PVTOL.mat';

W = @(x) W_mat(x);
dW = @(x) {zeros(6), zeros(6), dW_p_mat(x),...
           dW_vy_mat(x),zeros(6), zeros(6)};

sigma_ThBw = 0.1193;
w_upper = 898.581;
w_lower = 9.1692;
lambda =  0.56;

%% Dynamics

mass = 0.486;
J = 0.00383;
g = 9.81;
len = 0.25;

f  = @(x) [x(4)*cos(x(3)) - x(5)*sin(x(3));
           x(4)*sin(x(3)) + x(5)*cos(x(3));
           x(6);
           x(6)*x(5)-g*sin(x(3));
           -x(6)*x(4)-g*cos(x(3));
           0];
       
B = [zeros(1,4),1/mass, len/J;
     zeros(1,4),1/mass,-len/J]';

df = @(x) [0,0,-x(4)*sin(x(3))-x(5)*cos(x(3)),cos(x(3)),-sin(x(3)),0;
           0,0, x(4)*cos(x(3))-x(5)*sin(x(3)),sin(x(3)), cos(x(3)),0;
           zeros(1,5),1;
           0,0,-g*cos(x(3)), 0, x(6), x(5);
           0,0, g*sin(x(3)), -x(6), 0, -x(4);
           zeros(1,6)];

B_w = [zeros(1,3),1,0,0;
       zeros(1,3),0,1,0]';

%% Bounds

w_max = 0.02;

M_ccm = eye(n);
d_bar = (w_max*sigma_ThBw/lambda);
euc_bound = d_bar*sqrt(w_upper);
ctrl_bound = 9.1*w_max;

obs.infl = euc_bound + len;

P = eye(n);
alpha = 1e-3;

state_constr_low = -[10;10;pi/4;3;1;pi/3.5]+euc_bound;
ctrl_constr = [0.1*mass*g+ctrl_bound, 3*mass*g-ctrl_bound;
               0.1*mass*g+ctrl_bound, 3*mass*g-ctrl_bound];

%% Setup Geodesic numerics

geodesic_N = 2;

[geo_Prob,geo_Ke,geo_we,T_e,T_dot_e,geo_Aeq] = ...
        setup_geodesic_calc(n,geodesic_N,W,dW,state_constr_low-euc_bound);

% Assemble geodesic struct for MPC
geodesic_MPC = struct('geo_Prob',geo_Prob,'W',W,'geodesic_N',geodesic_N,'T_e',T_e,'T_dot_e',T_dot_e,...
                      'geo_Aeq',geo_Aeq);

%% Setup NMPC numerics

Tp = 6;
delta = 1.0;
dt = 0.002;

x_eq = [0;0;0;0;0;0];
u_eq = [0.5*mass*g; 0.5*mass*g];

N_mpc = 100;

[NMPC_Prob,L_e,L_e_full,MPC_st] = setup_NMPC(n,m,...
    f,B,df, state_constr_low,ctrl_constr,...
    N_mpc,Tp,delta,dt,...
    P,alpha,geodesic_MPC,d_bar^2,...
    x_eq,u_eq,obs);

%% Test MPC Solve

test_state = [-3;
              -3;
               0;
               0.5;
               0;
               0]; 
      
warm = struct('sol',0,'shift',0,...
              's_t',MPC_st,'Tp',Tp,...
              'state',[],'ctrl',[]);
tic
[NMPC_state,NMPC_ctrl,warm.state,warm.ctrl,converged_MPC] = compute_NMPC(NMPC_Prob,...
    test_state,test_state,state_constr_low,x_eq,u_eq,...
    n,m,N_mpc,L_e_full,warm);
toc
disp(converged_MPC);

warm.sol = 1;

% Visualize
close all
figure(); 
hold on
%RCI set
Ellipse_plot(M_ccm(1:2,1:2)*(1/euc_bound^2), NMPC_state(1,1:2)',25,'k');
%Terminal set
Ellipse_plot(P(1:2,1:2)*(1/alpha),x_eq(1:2),25,'r');

%Full MPC traj
plot(NMPC_state(:,1),NMPC_state(:,2),'r-','linewidth',2);
quiver(NMPC_state(1:(0.1/dt):end-1,1),NMPC_state(1:(0.1/dt):end-1,2),...
       -0.1*sin(NMPC_state(1:(0.1/dt):end-1,3)),0.1*cos(NMPC_state(1:(0.1/dt):end-1,3)));

grid on
axis equal
xlabel('X'); ylabel('Z','interpreter','latex');
set(findall(gcf,'type','text'),'FontSize',32);set(gca,'FontSize',32)

pause;

%% Test Geodesic Numerics

% test_rand = randn(n,1); test_rand = test_rand/norm(test_rand);
% test_state = NMPC_state(1,:)' + (d_bar*sqrt(w_lower))*test_rand;
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

pause;
                
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

dt_MPC = delta;
solve_MPC = (0:dt_MPC:t_end)';
T_steps_MPC = length(solve_MPC)-1;

MPC_state = cell(T_steps_MPC,1);
MPC_ctrl = cell(T_steps_MPC,1);

x_act = zeros(T_steps+1,n);
x_act(1,:) = test_state';

Geod = cell(T_steps,1);

Aux_ctrl = zeros(T_steps,m);
True_ctrl = zeros(T_steps,m);

ctrl_solve_time = zeros(T_steps,3);
ctrl_solve_time(:,1) = NaN;

opt_solved = NaN(T_steps,3);

geo_energy = zeros(T_steps,2);
geo_energy(:,2) = NaN;

state_0 = test_state;
state_0_MPC = NMPC_state(1,:)';

i_mpc = 0;

      
%% Simulate

for i = 1:T_steps
    fprintf('%d/%d \n',i,T_steps);
    %First Solve MPC
    if (mod(solve_t(i),delta)==0)
        [~, ~,J_opt,~] = compute_geodesic_tom(geo_Prob,...
            n,geodesic_N,state_0_MPC,state_0,T_e,T_dot_e,geo_Aeq);
        geo_energy(i,1) = J_opt;
        
        tic
        [NMPC_state,NMPC_ctrl,warm.state,warm.ctrl,opt_solved(i,1)] = ...
         compute_NMPC(NMPC_Prob,...
            state_0,state_0_MPC,state_constr_low,x_eq,u_eq,...
            n,m,N_mpc,L_e,warm);
        ctrl_solve_time(i,1) = toc;
        
        warm.shift = delta;
        
        i_mpc = i_mpc + 1;
        
        MPC_state{i_mpc} = NMPC_state;
        MPC_ctrl{i_mpc} = NMPC_ctrl;
        
        x_nom = MPC_state{i_mpc}(1,:);
        u_nom = MPC_ctrl{i_mpc}(1:2,:);
        
        [~, ~,J_opt,~] = compute_geodesic_tom(geo_Prob,n,geodesic_N,...
            x_nom',state_0,T_e,T_dot_e,geo_Aeq);
        
        geo_energy(i,2) = J_opt;
        
        %update starting state for next MPC problem
        state_0_MPC = MPC_state{i_mpc}(end,:)';
    else
        i_mpc_use = round((mod(solve_t(i),delta))/dt)+1;
        x_nom = MPC_state{i_mpc}(i_mpc_use,:);
        u_nom = MPC_ctrl{i_mpc}(i_mpc_use:i_mpc_use+1,:);
    end
    
    %Optimal Control
    tic
    [X, X_dot,J_opt,opt_solved(i,2)] = compute_geodesic_tom(geo_Prob,...
        n,geodesic_N,x_nom',state_0,T_e,T_dot_e,geo_Aeq);
    ctrl_solve_time(i,2) = toc;
    
    Geod{i} = X';
    geo_energy(i,1) = J_opt;
    
    tic
    [Aux_ctrl(i,:),opt_solved(i,3)] = compute_opt_aux(aux_Prob,geo_Ke,X,X_dot,J_opt,...
                            W,f,B,u_nom(1,:)',lambda);
    ctrl_solve_time(i,3) = toc;
    
    True_ctrl(i,:) = u_nom(1,:)+Aux_ctrl(i,:);
    
    %Simulate Optimal
    %     w_dist = w_max*[sin(2*pi*solve_t(i))*cosd(30);
    %                     cos(2*pi*solve_t(i))*sind(30)];
    
    %     dist_dir = (X_dot(:,geo_Ke+1))'*(W(state_0)\eye(4))*B_w;
    %     w_dist = w_max*(dist_dir'/norm(dist_dir));
    w_dist = zeros(2,1);
    
    [d_t,d_state] = ode45(@(t,d_state)ode_sim(t,d_state,[solve_t(i);solve_t(i+1)],u_nom,Aux_ctrl(i,:),...
        f,B,B_w,w_dist),[solve_t(i),solve_t(i+1)],state_0,ode_options);

    state_0 = d_state(end,:)';
    x_act(i+1,:) = state_0';
end

%% Plot
close all

% State Trajectory
figure()
hold on
plot(solve_t, x_act(:,3:6),'linewidth',2);
grid on
xlabel('Time [s]'); ylabel('States'); 
h_leg = legend('$\phi$','$v_x$','$v_z$','$\dot{\phi}$');
set(h_leg,'interpreter','latex');
set(findall(gcf,'type','text'),'FontSize',32);set(gca,'FontSize',32)

%Control Trajectory
figure()
hold on
plot(solve_t(1:end-1),True_ctrl(:,1)-Aux_ctrl(:,1),'b:','linewidth',2);
plot(solve_t(1:end-1),True_ctrl(:,2)-Aux_ctrl(:,2),'r:','linewidth',2);
plot(solve_t(1:end-1),True_ctrl(:,1),'b-','linewidth',2);
plot(solve_t(1:end-1),True_ctrl(:,2),'r-','linewidth',2);
xlabel('Time [s]');
ylabel('u(t)'); 
set(findall(gcf,'type','text'),'FontSize',32);set(gca,'FontSize',32)
grid on

figure()
plot(solve_t(1:end-1),Aux_ctrl,'b-','linewidth',2); 
xlabel('Time [s]');
ylabel('k(x^{*},x)');
set(findall(gcf,'type','text'),'FontSize',32);set(gca,'FontSize',32)
grid on

%2D State Plot
figure()
hold on
for i_mpc = 1:T_steps_MPC
    plot(MPC_state{i_mpc}(:,1),MPC_state{i_mpc}(:,2),'r--','linewidth',2);
end
plot(x_act(:,1),x_act(:,2),'r-','linewidth',2);

Ellipse_plot(P(1:2,1:2)*(1/(alpha)), x_eq(1:2),30,'k');
xlabel('$X$','interpreter','latex'); 
ylabel('$Z$','interpreter','latex');
set(findall(gcf,'type','text'),'FontSize',32);set(gca,'FontSize',32)
grid on; 
axis equal

%Solve Time
figure()
hold on
plot(solve_t(1:end-1),ctrl_solve_time(:,1),'ro','markersize',10,'markerfacecolor','g');
plot(solve_t(1:end-1),ctrl_solve_time(:,2),'bs','markersize',10,'markerfacecolor','m');
plot(solve_t(1:end-1),ctrl_solve_time(:,3),'rd','markersize',10,'markerfacecolor','k');
grid on
legend('MPC','Geo','Aux');
xlabel('Time [s]');
ylabel('Solve time [s]'); title('Solve time');
set(findall(gcf,'type','text'),'FontSize',32);set(gca,'FontSize',32)

%Solve success
figure()
hold on
plot(solve_t(1:end-1),opt_solved(:,1),'ro','markersize',15,'markerfacecolor','g');
plot(solve_t(1:end-1),opt_solved(:,2),'bs','markersize',15,'markerfacecolor','m');
plot(solve_t(1:end-1),opt_solved(:,3),'rd','markersize',10,'markerfacecolor','k');
grid on
legend('MPC (1,2,3)','Geodesic (0,1,6)','Aux (0)');
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

% save('NMPC_allgower_alg_run.mat');



