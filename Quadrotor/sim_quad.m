clear; close all; clc;
warning('off','MATLAB:nargchk:deprecated');

%% Load config file

load_quad_config;

%% Generate trajectory

dt = 1/500;
[t,state_nom, ctrl_dyn_nom,thrust_nom, accel_nom] = generate_quad_traj_10(dt,Jq,mq,g,poly_file);
Tp = t(end);

% th_nom: thrust
% ctrl_dyn_nom: th_dot, rd, pd, yd
% state_nom: x,y,z, vx, vy,vz, r, p, y

%Generate timeseries data for simulation
thrust_dot = ctrl_dyn_nom(:,1)/mq;
MP_state = [state_nom(:,1:6),thrust_nom/mq,state_nom(:,7:8)]; %xc_nom : [pos,vel,th,r,p]
MP_ctrl = ctrl_dyn_nom; %[uc_nom, yaw_dot_nom]
MP_yaw = state_nom(:,9);

%% Visualize

visualize_quad;
keyboard;

%% Setup geodesic solver

geodesic_N = 2;
[geo_Prob,geo_Ke,T_e,T_dot_e,geo_Aeq] = ...
        setup_geodesic_calc(n,geodesic_N,W_fnc,dW_fnc,n_W);
    
geo_solver = 'npsol';    
    
geo_warm = struct('sol',0,'result',[]);

%% Setup simulation

%Contraction controller rate
dt_sim = 1/500;

%Initial conditions
thrust_init = g;
x_init = [state_nom(1,:)'+[0.01*randn(6,1);zeros(3,1)];
         zeros(3,1);
          thrust_init];

xc_init = [x_init(1:8); thrust_init];

%% Test Geodesic numerics

tic
[Xc, Xc_dot,J_opt,converged_geo,geo_result,geo_Prob] = ...
    compute_geodesic_tom(geo_Prob,n,geodesic_N,...
            MP_state(1,:)',xc_init,...
            T_e,T_dot_e,geo_Aeq,geo_warm,geo_solver);
toc;
disp('Geo dist: ');disp(converged_geo);
disp(sqrt(J_opt));
geo_Prob.CHECK = 1;
geo_warm.sol = 1;
geo_warm.result = geo_result;

%% Setup Auxiliary controller

tic
[ctrl_opt,converged_aux] = compute_opt_aux(geo_Ke,Xc,Xc_dot,J_opt,...
                            W_fnc,f_ctrl,B_ctrl,MP_ctrl(1,1:3)',lambda);
toc;
disp('opt_control:');disp(converged_aux);
disp(ctrl_opt);

%% Setup nonlinear sim

ode_options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);

t_end = Tp;
solve_t = (0:dt_sim:t_end)';
T_steps = length(solve_t)-1;

%Store disturbance
w_dist = zeros(T_steps,3);

%Store actual state
x_act = zeros(T_steps+1,13);

%Store geodesics
Geod = cell(T_steps,1);

%Store control history
Aux_ctrl = zeros(T_steps,4);
True_ctrl = zeros((t_end/dt)+1,4);
Nom_ctrl = zeros((t_end/dt)+1,4);

%Computation times
ctrl_solve_time = NaN(T_steps,2);

%Solve success
opt_solved = NaN(T_steps,2);

%Geodesic distances
geo_energy = NaN(T_steps,1);

%Initialize
x_act(1,:) = x_init';
state = x_init;
state_xc = xc_init;
yaw = x_init(9);

%% Simulate

disp('Ready to Simulate');
keyboard;

for i = 1:T_steps
    
    xc_nom = MP_state(1+(i-1)*(dt_sim/dt),:);
    yaw_nom = MP_yaw(1+(i-1)*(dt_sim/dt));
    uc_nom = MP_ctrl(1+(i-1)*(dt_sim/dt):1+i*(dt_sim/dt),:);
    
    %Optimal Control
    tic
    [Xc, Xc_dot,J_opt,opt_solved(i,1),geo_result,geo_Prob] = compute_geodesic_tom(geo_Prob,...
        n,geodesic_N,xc_nom',state_xc,T_e,T_dot_e,geo_Aeq,geo_warm,geo_solver);
    ctrl_solve_time(i,1) = toc;
    
    Geod{i} = Xc';
    geo_energy(i,1) = J_opt;
    geo_warm.result = geo_result;
    geo_warm.sol = 1;
    
    tic
    [Aux_ctrl(i,1:3),opt_solved(i,2)] = compute_opt_aux(geo_Ke,Xc,Xc_dot,J_opt,...
        W_fnc,f_ctrl,B_ctrl,uc_nom(1,1:3)',lambda);
    Aux_ctrl(i,4) = 2*lambda*(yaw_nom-yaw);
    ctrl_solve_time(i,2) = toc;
    
    True_ctrl(1+(i-1)*(dt_sim/dt):1+i*(dt_sim/dt),:) = uc_nom+kron(ones((dt_sim/dt)+1,1),Aux_ctrl(i,:));
    Nom_ctrl(1+(i-1)*(dt_sim/dt):1+i*(dt_sim/dt),:) = uc_nom;
    
    %Disturbance model
    w_dist(i,:) = w_max*(1/sqrt(3))*[1;1;1]';
    
    [d_t,d_state] = ode113(@(t,d_state)quad_ode(t,d_state,[solve_t(i):dt:solve_t(i+1)]',uc_nom,Aux_ctrl(i,:),...
        f,B,B_w,w_dist(i,:)'),[solve_t(i),solve_t(i+1)],state,ode_options);
    
    state = d_state(end,:)';
    x_act(i+1,:) = state';
    state_xc = [state(1:6);state(13);state(7:8)];
    yaw = wrapToPi(state(9));
end


%% Plot

close all
plot_quad;


