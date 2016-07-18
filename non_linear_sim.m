clear all; close all; clc;

%% Constants
n = 2;
m = 1;

%% Setup Geodesic Numerics

W = @(x) [4.258279173109496,  -0.934234854771844;
        -0.934234854771844,   3.766923169705589];
dW = @(x){zeros(2), zeros(2)};

geodesic_N = 4;

[geo_Prob,geo_Ke,geo_we,T_e,T_dot_e,geo_Aeq] = ...
        setup_geodesic_calc(n,geodesic_N,W,dW);

%% Setup Auxiliary controller
Q = kron(geo_we',eye(m)); 
Q_bar = Q'*Q;
aux_Prob = setup_opt_aux(m);

%% Test Auxiliary control computation
start_p = [0;0];
end_p = [3;3];

% Test TOMLAB
tic
[X, X_dot,J_opt,converged_geo] = ...
    compute_geodesic_tom(geo_Prob,n,geodesic_N,...
            start_p,end_p,...
            T_e,T_dot_e,geo_Aeq);
toc;
disp(converged_geo);

% Validate solution
Eps =0;
for k = 1:geo_Ke+1
    x_k = X(:,k);
    W_geo = W(x_k);
    M_xdot = W_geo\X_dot(:,k);
    e_k =  (X_dot(:,k)'*M_xdot);
    Eps = Eps + (1/2)*geo_we(k)*(abs(e_k-J_opt)/J_opt);
end

disp('Exit flag'); disp(exitflag);
disp('Eps'); disp(Eps);

%Numerics for optimal differential controller

f  = @(x) [-1*x(1) + 2*x(2);
           -3*x(1) + 4*x(2) - 0.25*(x(2)^3)];
B = [0.5;-2];
B_w = [0;1];
w_dist = 0;

df = @(x) [-1, 2;
           -3, 4-0.75*(x(2)^2)];

tic
[~,converged_aux] = compute_opt_aux(aux_Prob,geo_Ke,X,X_dot,J_opt,...
                            W,f,B,start_p,lambda);
toc;
disp(converged_aux);

%% Simulate non-linear system

ode_options = odeset('RelTol', 1e-4, 'AbsTol', 1e-7);

dt = 0.002;
solve_t = (0:dt:8)';
T_steps = length(solve_t)-1;

t_opt = cell(T_steps,1);
t_s = cell(T_steps,1);

state_0_opt = end_p;
state_0_s = end_p;

state_opt = cell(T_steps,1);
state_s = cell(T_steps,1);

ctrl_opt = zeros(T_steps,m);
ctrl_s = zeros(T_steps,m);

J_opt = zeros(T_steps,1);
J_s = zeros(T_steps,1);

ctrl_solve_time = zeros(T_steps,2);

solved = ones(T_steps,2);

% rho = @(x) 51.88116196-0.001458799296458*x(1)-...
%        0.001463848144924*x(2)+2.241303600302532*x(2)^2+...
%        7.261640169359435*x(1)^2+1.585568386570665*x(1)*x(2)+...
%        2.655722707839840*x(1)^4+0.611430475833738*x(1)^3*x(2)+...
%        4.105229604485951*x(1)^2*x(2)^2+...
%        2.152384589971836*x(1)*x(2)^3+...
%        6.938732192947612*x(2)^4;

% rho =@(x) 5825.652472-1.994747096772980*x(1)-0.169841718823589*x(2)+...
%           3.817406236202479e+03*x(1)^2+4.225004246942979e+03*x(2)^2+...
%           7.486700415836236e+03*x(1)^4+6.462177245609084e+03*x(1)^2*x(2)^2+...
%           9.390587718981402e+03*x(2)^4+18.253220705351527*x(1)^2*x(2)+...
%           17.584174648858632*x(2)^3-4.746329713593079e+02*x(1)*x(2)-...
%           21.122699043530364*x(1)^3-32.298861002297762*x(1)*x(2)^2-...
%           1.930025253840884e+03*x(1)^3*x(2)-2.349756284415636e+03*x(1)*x(2)^3;

rho = @(x) 51.13139186-0.001647772177198*x(1)-0.001321483617423*x(2)+...
    0.933737285292152*x(2)^2+7.759857665781637*x(1)^2+0.948789019015335*x(1)*x(2)+...
    2.612862946905970*x(1)^4+0.411362758122320*x(1)^3*x(2)+4.092350682546736*x(1)^2*x(2)^2+...
    3.739815009936324*x(1)*x(2)^3+11.498438444671228*x(2)^4;
 %% 
for i = 1:T_steps
    %Optimal Control
    tic
    [X, X_dot,J_opt(i),solved(i,1)] = ...
    compute_geodesic_tom(geo_Prob,n,geodesic_N,...
            start_p,state_0_opt,...
            T_e,T_dot_e,geo_Aeq);
    ctrl_solve_time(i,1) = toc;
    
    tic
    [ctrl_opt(i,:),solved(i,2)] = compute_opt_aux(aux_Prob,geo_Ke,X,X_dot,J_opt,...
                            W,f,B,start_p,lambda);
    ctrl_solve_time(i,2) = toc;
    
    %Standard control
    [X_s, X_dot_s,J_s(i),exitflag] = compute_geodesic_tom(geo_Prob,n,geodesic_N,...
            start_p,state_0_s,...
            T_e,T_dot_e,geo_Aeq);
    ctrl_s(i,:) = aux_control(rho,X_s,X_dot_s,geo_Ke,geo_we,...
        B,W,m);
    
    %Simulate Optimal
    [d_t,d_state] = ode113(@(t,d_state)ode_sim(t,d_state,ctrl_opt(i,:)',f,B,B_w,w_dist),...
        [solve_t(i),solve_t(i+1)],state_0_opt,ode_options);
    t_opt{i} = d_t;
    state_opt{i} = d_state;
    state_0_opt = d_state(end,:)';
    
    %Simulate Standard
    [d_t_s,d_state_s] = ode113(@(t,d_state_s)ode_sim(t,d_state_s,ctrl_s(i,:)',f,B,B_w,w_dist),...
        [solve_t(i),solve_t(i+1)],state_0_s,ode_options);
    t_s{i} = d_t_s;
    state_s{i} = d_state_s;
    state_0_s = d_state_s(end,:)';
end

%% plot
close all;
figure()
subplot(2,1,1)
hold on
for i = 1:T_steps
    plot(t_s{i},state_s{i}(:,1),'r-','linewidth',2);
    plot(t_s{i},state_s{i}(:,2),'b-','linewidth',2);
%     plot(t_opt{i},state_opt{i}(:,3),'k-','linewidth',2);
end
grid on; 
xlabel('Time [s]'); 
set(findall(gcf,'type','text'),'FontSize',32);set(gca,'FontSize',32)

subplot(2,1,2)
hold on
for i = 1:T_steps
    plot(t_opt{i},state_opt{i}(:,1),'r-','linewidth',2);
    plot(t_opt{i},state_opt{i}(:,2),'b-','linewidth',2);
%     plot(t_opt{i},state_opt{i}(:,3),'k-','linewidth',2);
end
grid on; 
xlabel('Time [s]');
set(findall(gcf,'type','text'),'FontSize',32);set(gca,'FontSize',32)

figure()
stairs(solve_t,ctrl_s,'r-','linewidth',2); hold on
stairs(solve_t,ctrl_opt,'b-','linewidth',2);
xlabel('Time [s]');
ylabel('u(t)'); legend('\rho-multiplier','Optimized control');
grid on
set(findall(gcf,'type','text'),'FontSize',32);set(gca,'FontSize',32)

figure()
plot(solve_t(1:end-1),ctrl_solve_time(:,1),'go','markersize',10,'markerfacecolor','g');
hold on
plot(solve_t(1:end-1),ctrl_solve_time(:,2),'ro','markersize',10,'markerfacecolor','r');
grid on
legend('Geodesic computation time','Aux computation');
xlabel('Time [s]');
ylabel('Solve time [s]');
set(findall(gcf,'type','text'),'FontSize',32);set(gca,'FontSize',32)

figure()
plot(solve_t(1:end-1),solved(:,1),'go','markersize',10,'markerfacecolor','g');
hold on
plot(solve_t(1:end-1),solved(:,2),'ro','markersize',10,'markerfacecolor','r');
grid on
legend('Geodesic Solved','Aux Solved');
xlabel('Time [s]');
set(findall(gcf,'type','text'),'FontSize',32);set(gca,'FontSize',32)


figure()
plot(solve_t(1:end-1),J_s,'r-','linewidth',2); hold on
plot(solve_t,J_opt(1)*exp(-2*lambda*solve_t),'k--','linewidth',4); 
plot(solve_t(1:end-1),J_opt,'b-','linewidth',2);
grid on
xlabel('Time [s]');
ylabel('Geodesic energy');
xlabel('Time [s]'); 
legend('\rho-multiplier','Contraction bound','Optimized control');
set(findall(gcf,'type','text'),'FontSize',32);set(gca,'FontSize',32)


%%

save('rho_diff_comp.mat');