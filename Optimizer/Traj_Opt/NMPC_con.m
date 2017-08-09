% function c = NMPC_con(xu,Prob)
function c = NMPC_con(xu,Prob,n,N,P,D,f,B,B_full,Tp)
%Dynamics, and terminal

obs = Prob.user.obs;

no = obs.n_obs;
global geodesic_MPC;

global US_A;
% global X_EQ;
% global X_ACT;
% global GEOD_ITER;
% global NMPC_GEOD;

c = zeros(n*(N+1)+2+no*(N+1),1);

%% Dynamics constraints

c(1:n*(N+1)) = (2/Tp)*D*xu(1:n*(N+1)) -...
    (NMPC_dyn(f,xu(1:n*(N+1)),n,N) + B_full*xu(n*(N+1)+1:end-1));

%% Initial RPI constraint

% tic
[~, X_dot,J_opt,~,geo_result,~] = compute_geodesic_tom(geodesic_MPC.geo_Prob,...
    n,geodesic_MPC.geodesic_N,xu(1:n),Prob.user.x_act,geodesic_MPC.T_e,geodesic_MPC.T_dot_e,geodesic_MPC.geo_Aeq,geodesic_MPC.warm,geodesic_MPC.solver);
% NMPC_GEOD(GEOD_ITER,2) = toc;
US_A = X_dot(:,1);
geodesic_MPC.warm.sol = 1; geodesic_MPC.warm.result = geo_result;

c(n*(N+1)+1) = J_opt;
% NMPC_GEOD(GEOD_ITER,1) = J_opt;
% GEOD_ITER = GEOD_ITER +1;

%% Terminal constraints
t_star = xu(end);
x_eq = interp1(Prob.user.t_nom,Prob.user.x_nom,t_star);

c(n*(N+1)+2) = (xu(n*N+1:n*(N+1))-x_eq')'*P*(xu(n*N+1:n*(N+1))-x_eq');

%% Obstacle constraints

% c_obs = zeros(no*(N+1),1);
for k = 1:N+1
    for i = 1:no
        o_pos = obs.pos(:,i);
        Mo = obs.M{k}(:,:,i);
        x_k = xu(1+(k-1)*n:2+(k-1)*n);
        c(n*(N+1)+2+(i-1)*(N+1)+k,1) = (o_pos-x_k)'*Mo*(o_pos-x_k);
    end
end

end

