function conJ = NMPC_conJ(xu,Prob,n,N,P,D,df,~,Tp)
%Dynamics, and terminal
% n = Prob.user.n;
% N = Prob.user.N;
% m = Prob.user.m;

obs = Prob.user.obs;
no = obs.n_obs;

global US_A;
global geodesic_MPC;
global NMPC_CONJ;

% conJ = zeros(n*(N+1)+2+no*(N+1),(n+m)*(N+1));

%% Dynamics constraints

NMPC_CONJ(1:n*(N+1),1:n*(N+1)) = (2/Tp)*D - ...
            df_all(df,xu(1:n*(N+1)),n,N);

% conJ(1:n*(N+1),n*(N+1)+1:end) = -Prob.user.B_full;

%% Initial RPI constraint

w_poly = geodesic_MPC.W.w_poly_fnc(xu(1:n));
M = (geodesic_MPC.W.W_eval(w_poly))\eye(n);
NMPC_CONJ(n*(N+1)+1,1:n) = -2*US_A'*M;


%% Terminal constraint
NMPC_CONJ(n*(N+1)+2,n*N+1:n*(N+1)) = 2*(P*(xu(n*N+1:n*(N+1))-Prob.user.x_eq))';

%% Obstacle constraints

% cJ_obs = zeros(no*(N+1),(n+m)*(N+1));

for k = 1:N+1
    for i = 1:no
        o_pos = obs.pos(:,i);
        Mo = obs.M{k}(:,:,i);
        x_k = xu(1+(k-1)*n:2+(k-1)*n);
        NMPC_CONJ(n*(N+1)+2+(i-1)*(N+1)+k,1+(k-1)*n:2+(k-1)*n) = -2*(o_pos-x_k)'*Mo;
    end
end

conJ = NMPC_CONJ;


end

