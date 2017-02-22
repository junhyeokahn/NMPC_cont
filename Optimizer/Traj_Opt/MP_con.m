function c = MP_con(xu,Prob)
%Dynamics, and terminal

n = Prob.user.n;
N = Prob.user.N;

no = Prob.user.obs.n_obs;

global US_A;
global geodesic_MPC;

c = zeros(n*(N+1)+2+no*(N+1),1);

%% Dynamics constraints

c(1:n*(N+1)) = (2/Prob.user.Tp)*Prob.user.D*xu(1:n*(N+1)) -...
    (NMPC_dyn(Prob.user.f,xu(1:n*(N+1)),n,N) + Prob.user.B_full*xu(n*(N+1)+1:end));

%% Initial RPI constraint

[~, X_dot,J_opt,~,geo_result,~] = compute_geodesic_tom(geodesic_MPC.geo_Prob,...
    n,geodesic_MPC.geodesic_N,xu(1:n),Prob.user.x_act,geodesic_MPC.T_e,geodesic_MPC.T_dot_e,geodesic_MPC.geo_Aeq,geodesic_MPC.warm,geodesic_MPC.solver);
US_A = X_dot(:,1);
geodesic_MPC.warm.sol = 1; geodesic_MPC.warm.result = geo_result;

c(n*(N+1)+1) = J_opt;

%% Terminal constraint
c(n*(N+1)+2) = (xu(n*N+1:n*(N+1))-Prob.user.x_eq)'*Prob.user.P*(xu(n*N+1:n*(N+1))-Prob.user.x_eq);

%% Obstacle constraints

% c_obs = zeros(no*(N+1),1);
for i = 1:no
    o_pos = Prob.user.obs.pos(:,i);
    Mo = Prob.user.obs.M_obs(:,:,i);
    for k = 1:N+1
        x_k = xu(1+(k-1)*n:2+(k-1)*n);
%         c_obs((i-1)*(N+1)+k,1) = (o_pos-x_k)'*Mo*(o_pos-x_k);
        c(n*(N+1)+2+(i-1)*(N+1)+k,1) = (o_pos-x_k)'*Mo*(o_pos-x_k);
    end
end

end

