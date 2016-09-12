function c = NMPC_con(xu,Prob)
%Dynamics, and terminal

n = Prob.user.n;
N = Prob.user.N;

no = Prob.user.obs.n_obs;

geo = Prob.user.geo_MPC;

global US_A;
global US_B;

c = zeros(n*(N+1)+2+no*(N+1),1);

%% Dynamics constraints

c(1:n*(N+1)) = (2/Prob.user.Tp)*Prob.user.D*xu(1:n*(N+1)) -...
    (NMPC_dyn(Prob.user.f,xu(1:n*(N+1)),n,N) + Prob.user.B_full*xu(n*(N+1)+1:end));

%% Initial RPI constraint

[~, X_dot,J_opt,~,geo_result] = compute_geodesic_tom(geo.geo_Prob,...
    n,geo.geodesic_N,xu(1:n),Prob.user.x_act,geo.T_e,geo.T_dot_e,geo.geo_Aeq,geo.warm);
US_A = X_dot(:,1);
geo.warm.sol = 1; geo.result = geo_result;

c(n*(N+1)+1) = J_opt;

% c(end-1) = (Prob.user.x_act- xu(1:n))'*(Prob.user.M)*(Prob.user.x_act-xu(1:n));

%% Terminal constraint
c(n*(N+1)+2) = (xu(n*N+1:n*(N+1))-Prob.user.x_eq)'*Prob.user.P*(xu(n*N+1:n*(N+1))-Prob.user.x_eq);

%% Obstacle constraints

c_obs = zeros(no*(N+1),1);
for i = 1:no
    o_pos = Prob.user.obs.pos(:,i);
    Mo = Prob.user.obs.M_obs(:,:,i);
    for k = 1:N+1
        x_k = xu(1+(k-1)*n:2+(k-1)*n);
%         c_obs((i-1)*(N+1)+k,1) = exp(-(o_pos-x_k)'*Mo*(o_pos-x_k));
        c_obs((i-1)*(N+1)+k,1) = (o_pos-x_k)'*Mo*(o_pos-x_k);
    end
end
US_B = c_obs;

c(n*(N+1)+3:end) = c_obs;
end

