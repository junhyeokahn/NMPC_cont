clear all; close all; clc;
warning('off','MATLAB:nargchk:deprecated');
%% Constants
n = 6;
m = 2;

%% Obstacle info

obs_loc = [[3;-4],...
           [0.7;-3],...
           [-1;-0.5],...
           [2.5;-0.5],...
           [-4;1],...
           [1.0;1.7],...
           [2.5;3.8],...
           [-2.5;-4.5],...
           [-2;4]];
obs_loc_mpc = [[0.7;-3],...
           [-1;-0.5],...
           [2.5;-0.5],...
           [1.0;1.7],...
           [2.5;3.8],...
           [-2.5;-4.5]];
obs_rad = [1,0.9,0.8,1.2,1,0.9,0.5,1,0.6];
obs_rad_mpc = [0.9,0.8,1.2,0.9,0.5,1];
obs = struct('n_obs',length(obs_rad),'pos',obs_loc,'r',obs_rad);
obs_mpc = struct('n_obs',length(obs_rad_mpc),'pos',obs_loc_mpc,'r',obs_rad_mpc);

%% Setup Metric

load 'metric_PVTOL.mat';
% W_mat = @(x) W_mat(x);

% dW = @(x) {zeros(6), zeros(6), dW_p_mat(x),...
%            dW_vy_mat(x),zeros(6), zeros(6)};
dW = struct('p',dW_p_mat,'vy',dW_vy_mat);

sigma_ThBw = 0.3296;
lambda =  0.8283;
n_W = [3,4];

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

w_max = 0.1;

M_ccm = W_upper\eye(n);
d_bar = (w_max*sigma_ThBw/lambda);
ctrl_bound = 5.95*w_max;
euc_bound = d_bar*sqrt(diag(W_upper));

In = eye(n);
M_ccm_pos = (1/d_bar^2)*((In(1:2,:)*W_upper*In(1:2,:)')\eye(2));
[U_pos,S_pos,V_pos] = svd(M_ccm_pos);
    
%Rescale ellipsoids by obstacle + robot radius
M_obs = zeros(2,2,obs.n_obs);
for i = 1:obs.n_obs
    S_new = (sqrt(S_pos\eye(2)) + (obs_rad(i)+len)*eye(2))^2\eye(2);
    M_obs(:,:,i) = U_pos*S_new*V_pos';
end
obs.M_obs = M_obs;

M_obs_mpc = zeros(2,2,obs_mpc.n_obs);
for i = 1:obs_mpc.n_obs
    S_new = (sqrt(S_pos\eye(2)) + (obs_mpc.r(i)+len)*eye(2))^2\eye(2);
    M_obs_mpc(:,:,i) = U_pos*S_new*V_pos';
end
obs_mpc.M_obs = M_obs_mpc;

P = 2.5*eye(n);
alpha = 1e-3;

state_constr_low = -[5.5;5.5;pi/4;2;1;pi/3]+euc_bound;
ctrl_constr = [0.1*mass*g+ctrl_bound, 2*mass*g-ctrl_bound;
               0.1*mass*g+ctrl_bound, 2*mass*g-ctrl_bound];

%% Setup Geodesic numerics

geodesic_N = 2;

setup_geodesic_MPC(n,geodesic_N,W_mat,dW,n_W); %initializes geodesic_MPC struct
global geodesic_MPC;

[geo_Prob,geo_Ke,T_e,T_dot_e,geo_Aeq] = ...
        setup_geodesic_calc(n,geodesic_N,W_mat,dW,n_W);
    
geo_warm = struct('sol',0,'result',[]);    

%% Setup MP numerics

Tp = 28;
dt = 0.001;
N_mp = 120;

x_eq = [4.5;4.5;0;0;0;0];
u_eq = [0.5*mass*g; 0.5*mass*g];

T_mpc = 3;
dt_sim = 1/500;
delta = 1;
N_mpc = 14;

% Setup motion planning problem
[MP_Prob,L_e_mp,~] = setup_MP(n,m,...
    f,B,df, state_constr_low,ctrl_constr,...
    N_mp,Tp,dt,...
    P,alpha,d_bar^2,...
    x_eq,obs,'MP');

load MP_WARM.mat;

%% Test MP Solve

test_state = [-4.4;
              -5;
               0;
               1.3;
               0;
               0]; 
% test_state = mp_warm.state(1,:)';
      
tic
[MP_state,MP_ctrl,converged_MP,mp_warm] = compute_MP(MP_Prob,...
    test_state,test_state,state_constr_low,ctrl_constr,x_eq,u_eq,...
    n,m,N_mp,L_e_mp,mp_warm);
toc
disp('MP:'); disp(converged_MP);

mp_warm.sol = 1;

% Visualize
close all
figure(); 
hold on
%Obstacles
for i_ob = 1:obs.n_obs
    Ellipse_plot(eye(2)*(1/(obs.r(i_ob)+len)^2),obs.pos(:,i_ob), 25,'r',1);
end
%RCI set
for i = 1:(T_mpc/dt):length(MP_state)
    Ellipse_plot(M_ccm_pos,MP_state(i,1:2)',25,'k');
end
%Terminal set
Ellipse_plot(P(1:2,1:2)*(1/alpha),x_eq(1:2),25,'r');

%Full MP traj
plot(MP_state(:,1),MP_state(:,2),'b-','linewidth',1);
quiver(MP_state(1:(0.5/dt):end-1,1),MP_state(1:(0.5/dt):end-1,2),...
       (MP_ctrl(1:(0.5/dt):end-1,1)+MP_ctrl(1:(0.5/dt):end-1,2)).*sin(MP_state(1:(0.5/dt):end-1,3)),...
       (MP_ctrl(1:(0.5/dt):end-1,1)+MP_ctrl(1:(0.5/dt):end-1,2)).*cos(MP_state(1:(0.5/dt):end-1,3)));

plot(test_state(1),test_state(2),'go','markersize',10,'markerfacecolor','g');
grid on
axis equal
xlabel('X'); ylabel('Z','interpreter','latex');
set(findall(gcf,'type','text'),'FontSize',32);set(gca,'FontSize',32)

%% Test Geodesic Numerics

% test_rand = randn(n,1); test_rand = test_rand/norm(test_rand);
% test_state = NMPC_state(1,:)' + (d_bar*sqrt(w_lower))*test_rand;
% test_state = x_eq; %just to see cool curved lines

tic
[~, ~,J_opt,converged_geo,geo_result,geo_Prob] = ...
    compute_geodesic_tom(geo_Prob,n,geodesic_N,...
            MP_state(1,:)',test_state,...
            T_e,T_dot_e,geo_Aeq,geo_warm);
toc;
disp('Geo dist: ');disp(converged_geo);
disp(sqrt(J_opt));
geo_Prob.CHECK = 1;
geo_warm.sol = 1;
geo_warm.result = geo_result;

tic
[X, X_dot,J_opt,converged_geo,geo_result_MPC,geo_Prob_MPC] = ...
    compute_geodesic_tom(geodesic_MPC.geo_Prob,n,geodesic_N,...
            MP_state(1,:)',test_state,...
            T_e,T_dot_e,geo_Aeq,geodesic_MPC.warm);
toc;
disp('MPC Geo dist: '); disp(converged_geo);
disp(sqrt(J_opt));
geo_Prob_MPC.CHECK = 1;
geodesic_MPC.geo_Prob = geo_Prob_MPC;
geodesic_MPC.warm.sol = 1;
geodesic_MPC.warm.result = geo_result_MPC;

figure(1)
plot(X(1,:),X(2,:),'b-','linewidth',2);

% pause;

%% Setup MPC numerics

[MPC_Prob,L_e,L_e_mpc,MPC_st] = setup_NMPC(n,m,...
    f,B,df, state_constr_low,ctrl_constr,...
    N_mpc,T_mpc,delta,dt,...
    P,alpha,d_bar^2,M_ccm,...
    x_eq,obs_mpc,'MPC');

load MPC_WARM.mat;
% global NMPC_GEOD;
% NMPC_GEOD = zeros(200,2);
% global GEOD_ITER; 
% GEOD_ITER = 1;

% mpc_warm = struct('Tp',T_mpc,'shift',0,'sol',0,'solve_t',0,...
%                   's_t',MPC_st,'state',[],'ctrl',[],'result',[]);

%% Test MPC solve
tic
[MPC_state,~,converged_MPC,mpc_warm,MPC_Prob] = compute_NMPC(MPC_Prob,...
    test_state,MP_state(1,:)',state_constr_low,ctrl_constr,MP_state,MP_ctrl,...
    n,m,N_mpc,L_e_mpc,mpc_warm,dt);
toc
disp('MPC:');disp(converged_MPC);

MPC_Prob.CHECK = 1;
mpc_warm.sol = 1;

figure(1)
plot(MPC_state(:,1),MPC_state(:,2),'r-','linewidth',2);

% pause;
                
%% Setup Auxiliary controller
aux_Prob = setup_opt_aux(m);

tic
[ctrl_opt,converged_aux] = compute_opt_aux(aux_Prob,geo_Ke,X,X_dot,J_opt,...
                            W_mat,f,B,MP_ctrl(1,:)',lambda);
toc;
disp('opt_control:');disp(converged_aux);
disp(ctrl_opt);

% pause;
% keyboard;

%% Set up non-linear sim

ode_options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);

% dt_sim = 1/400;
t_end = 5*delta;
solve_t = (0:dt_sim:t_end)';
T_steps = length(solve_t)-1;

dt_MPC = delta;
solve_MPC = (0:dt_MPC:t_end)';
T_steps_MPC = length(solve_MPC)-1;

MPC_state = cell(T_steps_MPC,1);
MPC_ctrl = cell(T_steps_MPC,1);

accel_nom = zeros(T_steps,2);
w_dist = zeros(T_steps,2);

x_act = zeros(T_steps+1,n);
x_act(1,:) = test_state';

Geod = cell(T_steps,1);

Aux_ctrl = zeros(T_steps,m);
True_ctrl = zeros((t_end/dt)+1,m);
Nom_ctrl = zeros((t_end/dt)+1,m);

ctrl_solve_time = zeros(T_steps,3);
ctrl_solve_time(:,1) = NaN;

opt_solved = NaN(T_steps,3);

geo_energy = zeros(T_steps,2);
geo_energy(:,2) = NaN;

state_0 = test_state;
state_0_MPC = MP_state(1,:)';

i_mpc = 0;

      
%% Simulate

for i = 1:T_steps

    %First Solve MPC
    if (mod(solve_t(i),delta)==0)
        
        fprintf('%d/%d:',i,T_steps);
        
        [~, ~,J_opt,~,~,geo_Prob] = compute_geodesic_tom(geo_Prob,...
            n,geodesic_N,state_0_MPC,state_0,T_e,T_dot_e,geo_Aeq,geo_warm);
        geo_energy(i,1) = J_opt;

%         figure(1)
%         plot(x_act(1:i,1),x_act(1:i,2),'k','linewidth',2);
        
        tic
        [MPC_x,MPC_u,opt_solved(i,1),mpc_warm,MPC_Prob] = ...
         compute_NMPC(MPC_Prob,state_0,state_0_MPC,state_constr_low,ctrl_constr,MP_state,MP_ctrl,...
            n,m,N_mpc,L_e_mpc,mpc_warm,dt);
        ctrl_solve_time(i,1) = toc;
        
        fprintf('%d, %.2f \n', opt_solved(i,1),ctrl_solve_time(i,1));
        
        mpc_warm.solve_t = solve_t(i);
        mpc_warm.shift = delta;
        mpc_warm.sol = 1;
        
        i_mpc = i_mpc + 1;
        
        MPC_state{i_mpc} = MPC_x;
        MPC_ctrl{i_mpc} = MPC_u;
        
        x_nom = MPC_state{i_mpc}(1,:);
        u_nom = MPC_ctrl{i_mpc}(1:round(dt_sim/dt)+1,:);
        
        [~, ~,J_opt,~,geo_result,geo_Prob] = compute_geodesic_tom(geo_Prob,n,geodesic_N,...
            x_nom',state_0,T_e,T_dot_e,geo_Aeq,geo_warm);
        
        geo_energy(i,2) = J_opt;
        geo_warm.result = geo_result;
        
        %update starting state for next MPC problem
        state_0_MPC = MPC_state{i_mpc}(round(delta/dt)+1,:)';
    else
        i_mpc_use = round((mod(solve_t(i),delta))/dt)+1;
        x_nom = MPC_state{i_mpc}(i_mpc_use,:);
        u_nom = MPC_ctrl{i_mpc}(i_mpc_use:i_mpc_use+round(dt_sim/dt),:);
    end
    
    %Nominal acceleration
    accel_nom(i,:) = -g*[sin(x_nom(3)), cos(x_nom(3))] + ...
                    [0, [1/mass, 1/mass]*u_nom(1,:)'];
    
    %Optimal Control
    tic
    [X, X_dot,J_opt,opt_solved(i,2),geo_result,geo_Prob] = compute_geodesic_tom(geo_Prob,...
        n,geodesic_N,x_nom',state_0,T_e,T_dot_e,geo_Aeq,geo_warm);
    ctrl_solve_time(i,2) = toc;
    
    Geod{i} = X';
    geo_energy(i,1) = J_opt;
    geo_warm.result = geo_result;
    
    tic
    [Aux_ctrl(i,:),opt_solved(i,3)] = compute_opt_aux(aux_Prob,geo_Ke,X,X_dot,J_opt,...
                            W_mat,f,B,u_nom(1,:)',lambda);
    ctrl_solve_time(i,3) = toc;
    
    True_ctrl(1+(i-1)*(dt_sim/dt):1+i*(dt_sim/dt),:) = u_nom+kron(ones((dt_sim/dt)+1,1),Aux_ctrl(i,:));
    Nom_ctrl(1+(i-1)*(dt_sim/dt):1+i*(dt_sim/dt),:) = u_nom;
    
    %Simulate Optimal
%     w_dist(i,:) = w_max*[cos(x_act(i,3));
%                         -sin(x_act(i,3))]';
    
    dist_dir = (X_dot(:,geo_Ke+1))'*(W_mat(state_0)\eye(n))*B_w;
    w_dist(i,:) = w_max*(dist_dir/norm(dist_dir));
    %     w_dist = zeros(2,1);
    
    [d_t,d_state] = ode113(@(t,d_state)ode_sim(t,d_state,[solve_t(i):dt:solve_t(i+1)]',u_nom,Aux_ctrl(i,:),...
        f,B,B_w,w_dist(i,:)'),[solve_t(i),solve_t(i+1)],state_0,ode_options);

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
subplot(2,1,1)
hold on
plot(0:dt:t_end,Nom_ctrl(:,1),'b-','linewidth',2);
plot(0:dt:t_end,True_ctrl(:,1),'r-','linewidth',2);
xlabel('Time [s]');
ylabel('u(t)'); 
set(findall(gcf,'type','text'),'FontSize',32);set(gca,'FontSize',32)
grid on

subplot(2,1,2)
hold on
plot(0:dt:t_end,Nom_ctrl(:,2),'b-','linewidth',2);
plot(0:dt:t_end,True_ctrl(:,2),'r-','linewidth',2);
xlabel('Time [s]');
ylabel('u(t)'); 
set(findall(gcf,'type','text'),'FontSize',32);set(gca,'FontSize',32)
grid on

figure()
plot(solve_t(1:end-1),Aux_ctrl,'linewidth',2); 
xlabel('Time [s]');
ylabel('k(x^{*},x)');
set(findall(gcf,'type','text'),'FontSize',32);set(gca,'FontSize',32)
grid on

%2D State Plot
figure()
hold on
plot(MP_state(:,1),MP_state(:,2),'b--','linewidth',1.5);
for i_mpc = 1:T_steps_MPC
    plot(MPC_state{i_mpc}(1:(delta/dt),1),MPC_state{i_mpc}(1:(delta/dt),2),'r--','linewidth',2);
    Ellipse_plot(M_ccm_pos,MPC_state{i_mpc}(1,1:2)',30,'k');
end
plot(x_act(:,1),x_act(:,2),'k-','linewidth',2);

Ellipse_plot(P(1:2,1:2)*(1/(alpha)), x_eq(1:2),30,'k');
xlabel('$X$','interpreter','latex'); 
ylabel('$Z$','interpreter','latex');
set(findall(gcf,'type','text'),'FontSize',32);set(gca,'FontSize',32)
grid on; 
% axis equal

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

%% Compute nominal accelerations

figure()  
plot(solve_t(1:end-1), accel_nom,'linewidth',2); 
hold on
plot(solve_t(1:end-1), w_dist,'linewidth',2);
grid on; xlabel('Time [s]'); ylabel('Acceleration [m/s^2]');
hl = legend('$\bar{a}_y$','$\bar{a}_z$','$w_y$','$w_z$');
set(hl,'interpreter','latex');
set(findall(gcf,'type','text'),'FontSize',32);set(gca,'FontSize',32)

save('quad_sim.mat');

%% Evaluate cost

%Control cost weighting
R = eye(m);

J_cost = zeros((t_end/dt)+1,1);
for i = 1:(t_end/dt)+1
    J_cost(i) = True_ctrl(i,:)*R*True_ctrl(i,:)';
end

disp('NET COST:'); disp(trapz([0:dt:t_end],J_cost));

%Evaluate nominal cost
J_nom_mp = zeros((t_end/dt)+1,1);
for i = 1:(t_end/dt)+1
    J_nom_mp(i) = MP_ctrl(i,:)*R*MP_ctrl(i,:)';
end
disp('NOMINAL MP COST:'); disp(trapz([0:dt:t_end],J_nom_mp));

figure()
plot([0:dt:t_end],cumtrapz([0:dt:t_end],J_nom_mp),'r-','linewidth',2);
grid on
xlabel('Time [s]'); grid on;
hold on

J_nom_mpc = zeros((t_end/dt)+1,1);
for i = 1:(t_end/dt)+1
    J_nom_mpc(i) = Nom_ctrl(i,:)*R*Nom_ctrl(i,:)';
end
disp('NOMINAL MPC COST:'); disp(trapz([0:dt:t_end],J_nom_mpc));
plot([0:dt:t_end],cumtrapz([0:dt:t_end],J_nom_mpc),'b-','linewidth',2);

