function [NMPC_Prob,Aeq,L_e] = setup_MNMPC(n,m,f,B,df,...
            x_L, u_L,...
            N,P,alpha,Tp,dt,delta)

%% Constants

%State bounds
x_U = -x_L;

%Control bounds
u_U = -u_L;

%State cost weighting
Q = 0.5*eye(n);

%Control cost weighting
R = 1*eye(m);

%Number of collocation points
K = N;

%CGL nodes
[s_t,w] = clencurt(K); %t_t: [-1, 1] : <-> : [0, Tp]
s = fliplr(s_t); %t: [1, -1]

%% Final solution interpolation matrix

tau_e = 0:dt:delta;
s_e = (2*tau_e - Tp)/Tp; %[-1, s_delta]
% 
% tau_full = 0:dt:Tp;
% t_e_full = (2*tau_full - Tp)/Tp;
% 
%Lagrange polynomial evaluation at the interpolation points
L_e = compute_Lagrange(length(s_e)-1,N,s_e,s_t);

%% Get Differentiation matrix

D = ChebyshevDiffMatrix(N,s);
D = kron(D,eye(n));

%% Variables

%State node values: [x_t0,...,x_tN]
%Control node values: [u_t0,...,u_tN]

% n_vars = (N+1)*(n+m);

%% Define problem

Q_bar = kron(diag([w(1:end-1);0]),Q); R_bar = kron(diag(w),R);
Q_tilde = Q_bar + kron(diag([zeros(N,1);(2/Tp)]),P);

F = blkdiag(Q_tilde,R_bar);
    
B_full = kron(eye(N+1),B);

xu_L = [kron(ones(N+1,1),x_L);
        kron(ones(N+1,1),u_L)];
xu_U = [kron(ones(N+1,1),x_U);
        kron(ones(N+1,1),u_U)]; 
    
Aeq = [eye(n), zeros(n,N*n+(N+1)*m)];
b_L = zeros(n,1);
b_U = zeros(n,1);
    
MPC_cost = @(xu) (Tp/2)*xu'*F*xu;
MPC_grad = @(xu) Tp*F*xu;

c_L = zeros(n*(N+1)+1,1);
c_U = [zeros(n*(N+1),1);alpha];

xu0 = zeros((n+m)*(N+1),1);

Name = 'NMPC';
NMPC_Prob = conAssign(MPC_cost,MPC_grad,[],[],...
            xu_L,xu_U,Name, xu0,...
            [], 0, Aeq, b_L,b_U,...
            'NMPC_conMN','NMPC_conJMN',[],[],...
            c_L,c_U,...
            [],[],[],[]);
        
NMPC_Prob.user.D = D;
NMPC_Prob.user.n = n;
NMPC_Prob.user.m = m;
NMPC_Prob.user.N = N;
NMPC_Prob.user.f = f;
NMPC_Prob.user.df = df;
NMPC_Prob.user.B_full = B_full;
NMPC_Prob.user.Tp = Tp;
NMPC_Prob.user.P = P;
end