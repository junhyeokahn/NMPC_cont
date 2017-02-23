function [NMPC_Prob,L_e,L_e_full,s_t] = ...
    setup_NMPC(n,m,...
               f,B,df,state_con,u_con,...
               N,Tp,delta,dt,...
               P,alpha,RPI_bound,...
               x_eq,obs,Name)

%% Constants

%State bounds
x_L = state_con(:,1);
x_U = state_con(:,2);

%Control bounds
u_L = u_con(:,1);
u_U = u_con(:,2);

%State cost weighting
Q = diag(zeros(n,1));

%Control cost weighting
R = eye(m);

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

% n_vars = (N+1)*(n+m);

%% Define problem

% u_eq = zeros(m,1);

x_eq_all = kron(ones(N+1,1),zeros(n,1));
u_eq_all = kron(ones(N+1,1),zeros(m,1));

xu_eq = [x_eq_all;u_eq_all];

Q_bar = kron(diag(w),Q); R_bar = kron(diag(w),R);
% Q_tilde = Q_bar + kron(diag([zeros(N,1);(2/Tp)]),P);
Q_tilde = Q_bar;

F = sparse(blkdiag(Q_tilde,R_bar));
% F_pattern = sparse(F~=0);
    
B_full = sparse(kron(eye(N+1),B));

xu_L = [kron(ones(N+1,1),x_L);
        kron(ones(N+1,1),u_L)];
xu_U = [kron(ones(N+1,1),x_U);
        kron(ones(N+1,1),u_U)]; 
       
MPC_cost = @(xu) (Tp/2)*(xu-xu_eq)'*F*(xu-xu_eq);% + Obs_cost(xu,n,N,obs) ;
MPC_grad = @(xu) Tp*F*(xu-xu_eq);% + Obs_grad(xu,n,m,N,obs);
MPC_hess = @(xu) Tp*F;

MPC_con = @(xu,Prob) NMPC_con(xu,Prob,n,N,P,D,f,B_full,Tp,obs);
MPC_conJ = @(xu,Prob) NMPC_conJ(xu,Prob,n,N,P,D,df,B_full,Tp,obs);

c_L = [zeros(n*(N+1)+2,1); ones(obs.n_obs*(N+1),1)];
c_U = [zeros(n*(N+1),1);RPI_bound;alpha;Inf*ones(obs.n_obs*(N+1),1)];

xu0 = zeros((n+m)*(N+1),1);

global NMPC_CONJ;
NMPC_CONJ = zeros(n*(N+1)+2+obs.n_obs*(N+1),(n+m)*(N+1));
NMPC_CONJ(1:n*(N+1),n*(N+1)+1:end) = -B_full;

% Name = 'NMPC';
NMPC_Prob = conAssign(MPC_cost,MPC_grad,MPC_hess,[],...
            xu_L,xu_U,Name, xu0,...
            [], 0, [],[],[],...
            MPC_con,MPC_conJ,[],[],...
            c_L,c_U,...
            [],[],[],[]);
        
        
NMPC_Prob.SOL.optPar(10) = 1e-4;
NMPC_Prob.SOL.optPar(12) = 1e-4;
        
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
NMPC_Prob.user.x_eq = x_eq;
NMPC_Prob.user.obs = obs;

end