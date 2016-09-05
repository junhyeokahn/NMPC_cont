clear all; close all; clc;
warning('off','MATLAB:nargchk:deprecated');
%% Constants
n = 4;
m = 1;

%% Setup Geodesic Numerics

load 'metric_FLR.mat';

W = @(x) W_mat(wrapToPi(x(1)));
dW = @(x) {dW_x1_mat(wrapToPi(x(1))), zeros(4), zeros(4), zeros(4)};

w_lower = 8.0232;

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
     1];

df = @(x) [0, 1, 0, 0;
           (mass*g*l/I)*cos(x(1))-(sigma/I), 0, (sigma/I), 0;
           zeros(1,3), 1;
           (sigma/J), 0, -(sigma/J), -(b/J)];

B_w = [0,1,0,0;
       0,0,0,1]';

w_max = 0.1;

M_ccm = eye(4);
lambda =  3.35;
d_bar = (w_max*sqrt(1/w_lower)/lambda);

P = eye(4);
alpha = 1e-3;

Tp = 8;
delta = 0.1;
dt = 0.01;

N_mpc = 80;

state_constr_low = -[pi;0.77*pi;pi;0.77*pi];
ctrl_constr_low = -33*ones(m,1);

x_eq = [0;0;0;0];

[NMPC_Prob,L_e,L_e_full] = setup_NMPC(n,m,...
    f,B,df, state_constr_low,ctrl_constr_low,...
    N_mpc,Tp,delta,dt,...
    P,alpha,geodesic_MPC,d_bar^2,...
    x_eq);

%% Setup Auxiliary controller
aux_Prob = setup_opt_aux(m);

%% Setup MC

N_MC = 100;

ode_options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);

t_end = Tp;
solve_t = (0:dt:t_end)';
T_steps = length(solve_t)-1;

state_nom = cell(1,N_MC);
ctrl_nom = cell(1,N_MC);

state_act = cell(1,N_MC);
ctrl_act = cell(1,N_MC);

ctrl_aux = cell(1,N_MC);

ctrl_solve_time = cell(1,N_MC);
ctrl_solved = cell(1,N_MC);

Geod = cell(1,N_MC);
geo_energy = cell(1,N_MC);

disturbance = cell(1,N_MC);

for it = 1:N_MC
    fprintf('it: %d',it);
    if (mod(it,10)==0)
        save(sprintf('FLR_sim_MC_%d.mat',it));
    end
    
    found_traj = 0;
    while(~found_traj)
        link_ang = wrapToPi(pi/2 + pi*rand(1));
        mot_ang = abs(wrapToPi(link_ang + (-15+30*rand(1))*(pi/180)))*sign(link_ang);
        
        test_state = [link_ang;
            (-25 + 50*rand(1))*(pi/180);
            mot_ang;
            (-15 + 30*rand(1))*(pi/180)];
        
        [NMPC_state,NMPC_ctrl,converged_MPC] = compute_NMPC(NMPC_Prob,...
            test_state,test_state,state_constr_low,x_eq,...
            n,m,N_mpc,L_e_full);
        if (converged_MPC ==1) || (converged_MPC == 2) || (converged_MPC == 3)
            found_traj = 1;
            fprintf(' found nom \n');
        end
        
    end
    
    state_nom{it} = NMPC_state;
    ctrl_nom{it} = NMPC_ctrl;
    
    x_act = zeros(T_steps+1,n);
    
    start_type = rand(1);
    if (start_type <= 0.5)
        start_vec = randn(4,1); start_vec = start_vec/norm(start_vec);
        test_state = NMPC_state(1,:)' + (d_bar*sqrt(w_lower))*start_vec;
    end
    x_act(1,:) = test_state';
    
    u_act = zeros(T_steps,m);
    u_aux = zeros(T_steps,m);  
        
    solve_time = zeros(T_steps,3);
    solve_time(:,1) = NaN;   
    
    solved = NaN(T_steps,3);
        
    geod = cell(T_steps,1);
    energy = zeros(T_steps,1);
        
    dist_type = rand(1);
    if (dist_type <= 0.5)
        vec_dist = randn(1,2); vec_dist = vec_dist/norm(vec_dist);
        disturbance{it} = kron(ones(T_steps+1,1),w_max*vec_dist);        
    else
        vec_dist = randn(1,2); vec_dist = vec_dist/norm(vec_dist);
        phase_dist = atan2(vec_dist(2),vec_dist(1));
        disturbance{it}(:,1) = w_max*sin(2*pi*solve_t)*cos(phase_dist);
        disturbance{it}(:,2) = w_max*cos(2*pi*solve_t)*sin(phase_dist);
    end
    
    state_0 = test_state;
    
    %% Simulate

    for i = 1:T_steps
        %     fprintf('%d/%d \n',i,T_steps);
        
        x_nom = NMPC_state(i,:)';
        u_nom = NMPC_ctrl(i,:)';
        
        %Optimal Control
        tic
        [X, X_dot,J_opt,solved(i,2)] = compute_geodesic_tom(geo_Prob,...
            n,geodesic_N,x_nom,state_0,T_e,T_dot_e,geo_Aeq);
        solve_time(i,2) = toc;
        
        geod{i} = X';
        energy(i) = J_opt;
        
        tic
        [u_aux(i,:),solved(i,3)] = compute_opt_aux(aux_Prob,geo_Ke,...
            X,X_dot,J_opt,W,f,B,u_nom,lambda);
        solve_time(i,3) = toc;
        
        u_act(i,:) = u_nom'+u_aux(i,:);
        
        %Simulate Optimal
        w_dist = disturbance{it}(i,:)';
        
        [d_t,d_state] = ode113(@(t,d_state)ode_sim(t,d_state,u_act(i,:)',...
            f,B,B_w,w_dist),[solve_t(i),solve_t(i+1)],state_0,ode_options);
        state_0 = d_state(end,:)';
        x_act(i+1,:) = state_0';
    end
    
    %% Record
    state_act{it} = x_act;
    ctrl_act{it} = u_act;
    ctrl_aux{it} = u_aux;
    
    ctrl_solve_time{it} = solve_time;
    ctrl_solved{it} = solved;   
    
    Geod{it} = geod;
    geo_energy{it} = energy;
    
    
end

%% Plot
close all

%Initial states
figure()
hold on
for it = 1:N_MC
    plot(state_act{it}(1,1),state_act{it}(1,3),'rx','markersize',15);
end
xlabel('q'); ylabel('\theta');
title('Initial conditions');
grid on
set(findall(gcf,'type','text'),'FontSize',32);set(gca,'FontSize',32)

%Control Trajectory
figure()
for it = 1:N_MC
    plot(solve_t(1:end-1),ctrl_act{it},'linewidth',2);
    hold all
end
xlabel('Time [s]');
ylabel('u(t)'); 
set(findall(gcf,'type','text'),'FontSize',32);set(gca,'FontSize',32)
grid on

figure()
for it = 1:N_MC
    plot(solve_t(1:end-1),ctrl_aux{it},'linewidth',2); 
    hold all
end
plot(solve_t(1:end-1),7.15*ones(T_steps,1),'r-','linewidth',2);
plot(solve_t(1:end-1),-7.15*ones(T_steps,1),'r-','linewidth',2);
xlabel('Time [s]');
ylabel('Ancillary control k(x^{*},x)');
grid on
set(findall(gcf,'type','text'),'FontSize',32);set(gca,'FontSize',32)

%Solve Time
figure()
hold on
for it = 1:N_MC
    plot(solve_t(1:end-1),ctrl_solve_time{it}(:,2),'b-','linewidth',2);
    plot(solve_t(1:end-1),ctrl_solve_time{it}(:,3),'r-','linewidth',2);
end
grid on
legend('Geo','Aux');
xlabel('Time [s]');
ylabel('Solve time [s]'); title('Solve time');
set(findall(gcf,'type','text'),'FontSize',32);set(gca,'FontSize',32)

%Solve success
figure()
hold on
for it = 1:N_MC
    plot(solve_t(1:end-1),ctrl_solved{it}(:,2),'b-','linewidth',2);
    plot(solve_t(1:end-1),ctrl_solved{it}(:,3),'r-','linewidth',2);
end
grid on
legend('Geodesic (0,1,6)','Aux (0)');
xlabel('Time [s]');
title('Convergence');
set(findall(gcf,'type','text'),'FontSize',32);set(gca,'FontSize',32)
 
%Energy
figure()
for it = 1:N_MC
    plot(solve_t(1:end-1),sqrt(geo_energy{it}),'linewidth',2);
    hold all
end
plot(solve_t(1:end-1),d_bar*ones(T_steps,1),'r-','linewidth',2);
grid on
xlabel('Time [s]'); title('Energy');
set(findall(gcf,'type','text'),'FontSize',32);set(gca,'FontSize',32)

%Disturbance
figure()
hold on
for it = 1:N_MC
    plot(solve_t,disturbance{it}(:,1),'b-','linewidth',2);
    plot(solve_t,disturbance{it}(:,2),'r-','linewidth',2);
end
grid on
xlabel('Time [s]'); ylabel('Nm'); title('Disturbance');
legend('Link','Motor');
set(findall(gcf,'type','text'),'FontSize',32);set(gca,'FontSize',32)


