clear all; close all; clc;
warning('off','MATLAB:nargchk:deprecated');
         
%% Load non-linear system

load_FLR_config;

%% Setup Geodesic numerics

geodesic_N = 3;

setup_geodesic_MPC(n,geodesic_N,W_fnc,dW_fnc,n_W); %initializes geodesic_MPC struct
global geodesic_MPC;

[geo_Prob,geo_Ke,T_e,T_dot_e,geo_Aeq] = ...
        setup_geodesic_calc(n,geodesic_N,W_fnc,dW_fnc,n_W);
    
geo_solver = 'snopt';    
    
geo_warm = struct('sol',0,'result',[]);    

%% Setup MP numerics

Tp = 7;
dt = 0.001;
N_mp = 70;

dt_sim = 0.002;

% Setup motion planning problem
[MP_Prob,L_e_mp,MP_st] = setup_MP(n,m,...
    f,B,df,state_constr,ctrl_constr,...
    N_mp,Tp,dt,...
    P,alpha,(0.98*d_bar)^2,...
    x_eq,obs,'MP');

load MP_WARM_FLR.mat;
% mp_warm = struct('Tp',Tp,'shift',0,'sol',0,...
%                   's_t',MP_st,'state',[],'ctrl',[],'result',[]);

%% Test MP Solve
      
tic
[MP_state,MP_ctrl,converged_MP,mp_warm] = compute_MP(MP_Prob,...
    test_state,test_state,state_constr_low,ctrl_constr,x_eq,u_eq,...
    n,m,N_mp,L_e_mp,mp_warm);
toc
disp('MP:'); disp(converged_MP);

mp_warm.sol = 1;
save('MP_WARM_PVTOL.mat','mp_warm');

%% Visualize

visualize_FLR;

%% Test Geodesic Numerics

tic
[X, X_dot,J_opt,converged_geo,geo_result,geo_Prob] = ...
    compute_geodesic_tom(geo_Prob,n,geodesic_N,...
            MP_state(1,:)',test_state,...
            T_e,T_dot_e,geo_Aeq,geo_warm,geo_solver);
toc;
disp('Geo dist: ');disp(converged_geo);
disp(sqrt(J_opt));
geo_Prob.CHECK = 1;
geo_warm.sol = 1;
geo_warm.result = geo_result;

tic
[~, ~,J_opt,converged_geo,geo_result_MPC,geo_Prob_MPC] = ...
    compute_geodesic_tom(geodesic_MPC.geo_Prob,n,geodesic_N,...
            MP_state(1,:)',test_state,...
            T_e,T_dot_e,geo_Aeq,geodesic_MPC.warm,'npsol');
toc;
disp('MPC Geo dist: '); disp(converged_geo);
disp(sqrt(J_opt));
geo_Prob_MPC.CHECK = 1;
geodesic_MPC.geo_Prob = geo_Prob_MPC;
geodesic_MPC.warm.sol = 1;
geodesic_MPC.warm.result = geo_result_MPC;
         
%% Setup Auxiliary controller

tic
[ctrl_opt,converged_aux] = compute_opt_aux(geo_Ke,X,X_dot,J_opt,...
                            W_fnc,f,B,MP_ctrl(1,:)',lambda);
toc;
disp('opt_control:');disp(converged_aux);
disp(ctrl_opt);

%% Set up non-linear sim

ode_options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);

t_end = Tp;
solve_t = (0:dt_sim:t_end)';
T_steps = length(solve_t)-1;

w_dist = zeros(T_steps,2);

x_act = zeros(T_steps+1,n);

Geod = cell(T_steps,1);

Aux_ctrl = zeros(T_steps,m);
True_ctrl = zeros((t_end/dt)+1,m);

ctrl_solve_time = zeros(T_steps,3);
ctrl_solve_time(:,1) = NaN;

opt_solved = NaN(T_steps,3);

geo_energy = zeros(T_steps,2);
geo_energy(:,2) = NaN;

x_act(1,:) = test_state';
state_0 = test_state;

%% Simulate
disp('Ready to Simulate');
keyboard;

FL_ctrl = 0; %compare with Feedback linearization controller

for i = 1:T_steps
    
    x_nom = MP_state(1+(i-1)*(dt_sim/dt),:);
    u_nom = MP_ctrl(1+(i-1)*(dt_sim/dt):1+i*(dt_sim/dt),:);
    
    %Optimal Control
    tic
    [X, X_dot,J_opt,opt_solved(i,2),geo_result,geo_Prob] = compute_geodesic_tom(geo_Prob,...
        n,geodesic_N,x_nom',state_0,T_e,T_dot_e,geo_Aeq,geo_warm,geo_solver);
    ctrl_solve_time(i,2) = toc;
    
    Geod{i} = X';
    geo_energy(i,1) = J_opt;
    geo_warm.result = geo_result;
    geo_warm.sol = 1;
    
    if (~FL_ctrl)
        tic
        [Aux_ctrl(i,:),opt_solved(i,3)] = compute_opt_aux(geo_Ke,X,X_dot,J_opt,...
            W_fnc,f,B,u_nom(1,:)',lambda);
        ctrl_solve_time(i,3) = toc;
    else
        xi_nom = phi(x_nom');
        xi_act = phi(state_0);
        
        v_aux = -K_FL*(xi_act - xi_nom);
        u_net = (1/B_tilde)*(v_aux + (f_tilde(x_nom) + B_tilde*u_nom(1,:)') - (f_tilde(state_0)));
        Aux_ctrl(i,:) = u_net' - u_nom(1,:);
    end
    
    True_ctrl(1+(i-1)*(dt_sim/dt):1+i*(dt_sim/dt),:) = u_nom+kron(ones((dt_sim/dt)+1,1),Aux_ctrl(i,:));
    
    %Disturbance model
    dist_dir = (X_dot(:,geo_Ke+1))'*(W_mat(state_0)\eye(n))*B_w;
    w_dist(i,:) = w_max*(dist_dir/norm(dist_dir));
    w_dist(i,:) = -(w_max/sqrt(2))*[1,1];
    
    [d_t,d_state] = ode113(@(t,d_state)ode_sim(t,d_state,[solve_t(i):dt:solve_t(i+1)]',u_nom,Aux_ctrl(i,:),...
        f,B,B_w,w_dist(i,:)'),[solve_t(i),solve_t(i+1)],state_0,ode_options);
    
    state_0 = d_state(end,:)';
    x_act(i+1,:) = state_0';
end

%% Plots
close all;
plot_FLR;

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

%% Movie sim

plot_FLR_movie(solve_t,x_act,round(0.05/dt_sim));
