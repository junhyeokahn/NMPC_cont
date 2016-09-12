function conJ = NMPC_conJ(xu,Prob)%,D,n,m,N,df,B_full,Tp,P)
%Dynamics, and terminal
n = Prob.user.n;
N = Prob.user.N;
m = Prob.user.m;

no = Prob.user.obs.n_obs;

global US_A;
global US_B;

conJ = zeros(n*(N+1)+2+no*(N+1),(n+m)*(N+1));

%% Dynamics constraints

conJ(1:n*(N+1),1:n*(N+1)) = (2/Prob.user.Tp)*Prob.user.D - ...
            df_all(Prob.user.df,xu(1:n*(N+1)),n,N);

conJ(1:n*(N+1),n*(N+1)+1:end) = -Prob.user.B_full;

%% Initial RPI constraint

% conJ(end-1,1:n) = NaN;

M = (Prob.user.geo_MPC.W(xu(1:n)))\eye(n);
conJ(n*(N+1)+1,1:n) = -2*US_A'*M;


%% Terminal constraint
conJ(n*(N+1)+2,n*N+1:n*(N+1)) = 2*(Prob.user.P*(xu(n*N+1:n*(N+1))-Prob.user.x_eq))';

%% Obstacle constraints

cJ_obs = zeros(no*(N+1),(n+m)*(N+1));

for i = 1:no
    o_pos = Prob.user.obs.pos(:,i);
    Mo = Prob.user.obs.M_obs(:,:,i);
    for k = 1:N+1
        x_k = xu(1+(k-1)*n:2+(k-1)*n);
%         cJ_obs((i-1)*(N+1)+k,1+(k-1)*n:2+(k-1)*n) = US_B((i-1)*(N+1)+k)*(2)*(o_pos-x_k)'*Mo;
        cJ_obs((i-1)*(N+1)+k,1+(k-1)*n:2+(k-1)*n) = -2*(o_pos-x_k)'*Mo;
    end
end

conJ(n*(N+1)+3:end,:) = cJ_obs;


end

