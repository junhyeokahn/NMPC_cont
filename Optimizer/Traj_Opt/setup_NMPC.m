function [NMPC_Prob,L_e,L_e_full,s_t] = ...
    setup_NMPC(MP_x,MP_u,MP_t,n,m,...
               f,B,df,state_con,u_con,...
               N,Tp,delta,dt,...
               P,alpha,RPI_bound,...
               x_eq,obs,R,Name)

%% Constants

%State bounds
x_L = state_con(:,1);
x_U = state_con(:,2);

%Control bounds
u_L = u_con(:,1);
u_U = u_con(:,2);

%State cost weighting
Q = 0.5*diag([1;1;0;0;0;0]);

%Number of collocation points
K = N;

%CGL nodes
[s_t,w] = clencurt(K); %t_t: [-1, 1] : <-> : [0, Tp]
s = fliplr(s_t); %t: [1, -1]

%% Final solution interpolation matrix

tau_e = 0:dt:delta;
s_e = (2*tau_e - Tp)/Tp; %[-1, s_delta]

tau_full = 0:dt:Tp; 
s_e_full = (2*tau_full - Tp)/Tp; %[-1, 1]

%Lagrange polynomial evaluation at the interpolation points
L_e = compute_Lagrange(length(s_e)-1,N,s_e,s_t);
L_e_full = compute_Lagrange(length(s_e_full)-1,N,s_e_full,s_t);

%% Get Differentiation matrix

D = ChebyshevDiffMatrix(N,s); %arranged for forward time
D = sparse(kron(D,eye(n)));

%% Variables

%State node values: [x_t0,...,x_tN]
%Control node values: [u_t0,...,u_tN]
%Rejoin MP time: T_end

% n_vars = (N+1)*(n+m) + 1;

%% Define problem

% u_eq = zeros(m,1);

x_eq_all = kron(ones(N+1,1),x_eq);
u_eq_all = kron(ones(N+1,1),zeros(m,1));

xu_eq = [x_eq_all;u_eq_all];

Q_bar = kron(diag(w),Q); R_bar = kron(diag(w),R);
% Q_tilde = Q_bar + kron(diag([zeros(N,1);(2/Tp)]),P);
Q_tilde = Q_bar;

F = sparse(blkdiag(Q_tilde,R_bar));
% F_pattern = sparse(F~=0);
    
B_full = sparse(kron(eye(N+1),B));

xu_L = [kron(ones(N+1,1),x_L);
        kron(ones(N+1,1),u_L);
        0];
xu_U = [kron(ones(N+1,1),x_U);
        kron(ones(N+1,1),u_U);
        MP_t(end)]; 
       
MPC_cost = @(xu) (Tp/2)*(xu(1:end-1)-xu_eq)'*F*(xu(1:end-1)-xu_eq)-xu(end);
MPC_grad = @(xu) [Tp*F*(xu(1:end-1)-xu_eq);-1];
MPC_hess = @(xu) blkdiag(Tp*F,0);

MPC_con = @(xu,Prob) NMPC_con(xu,Prob,n,N,P,D,f,B,B_full,Tp);
MPC_conJ = @(xu,Prob) NMPC_conJ(xu,Prob,n,N,P,D,f,df,B,B_full,Tp);

c_L = [zeros(n*(N+1)+2,1); ones(obs.n_obs*(N+1),1)];
c_U = [zeros(n*(N+1),1);RPI_bound;alpha;Inf*ones(obs.n_obs*(N+1),1)];

xu0 = [zeros((n+m)*(N+1),1);1];

global NMPC_CONJ;
NMPC_CONJ = zeros(n*(N+1)+2+obs.n_obs*(N+1),(n+m)*(N+1)+1);
NMPC_CONJ(1:n*(N+1),n*(N+1)+1:end-1) = -B_full;

% Name = 'NMPC';
NMPC_Prob = conAssign(MPC_cost,MPC_grad,MPC_hess,[],...
            xu_L,xu_U,Name, xu0,...
            [], 0, [],[],[],...
            MPC_con,MPC_conJ,[],[],...
            c_L,c_U,...
            [],[],[],[]);
        
        
NMPC_Prob.SOL.optPar(10) = 1e-4;
NMPC_Prob.SOL.optPar(12) = 1e-4;

%% Create MPC obstacle structure
obs.M = cell(N+1,1);
for k = 1:N+1
    obs.M{k} = zeros(2,2,obs.n_obs);
end
 
%% Declare user variables

NMPC_Prob.user.x_act = zeros(n,1);
NMPC_Prob.user.D = D;
NMPC_Prob.user.n = n;
NMPC_Prob.user.m = m;
NMPC_Prob.user.N = N;
NMPC_Prob.user.f = f;
NMPC_Prob.user.df = df;
NMPC_Prob.user.B_full = B_full;
NMPC_Prob.user.Tp = Tp;
NMPC_Prob.user.P = P;
NMPC_Prob.user.obs = obs;
NMPC_Prob.user.x_nom = MP_x;
NMPC_Prob.user.u_nom = MP_u;
NMPC_Prob.user.t_nom = MP_t;

end