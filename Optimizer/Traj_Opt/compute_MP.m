function [MP_state,MP_ctrl,converged,warm] = ...
    compute_MP(Prob,act_p,mpc_state,state_constr,ctrl_constr,x_eq,u_eq,...
                 n,m,N,L_e,warm)

%adjust initial nominal state guess
for i = 1:n
    if mpc_state(i) > state_constr(i,2)
        mpc_state(i) = state_constr(i,2)-0.01*(state_constr(i,2)-state_const(i,1));
    elseif mpc_state(i) < state_constr(i,1)
        mpc_state(i) = state_constr(i,1)+0.01*(state_constr(i,2)-state_const(i,1));
    end
end

%Solution guess
if (warm.sol==0) %got nothing
    In = eye(n);
    x0 = zeros((N+1)*n,1);
    for i = 1:n
        x0 = x0 + kron(linspace(mpc_state(i),x_eq(i),N+1), In(i,:))';
    end
    
    Im = eye(m);
    u0 = zeros((N+1)*m,1);
    for j = 1:m
        u0 = u0 + kron(u_eq(j)*ones(N+1,1), Im(:,j));
    end
    
    Prob = modify_x_0(Prob,[x0;u0]);
elseif (warm.sol==0.5) %have some guess of homotopy
    Prob = modify_x_0(Prob,warm.result.x_k);
end

%Update constraint information
Prob.user.x_act = act_p;

%Recall warm solution
if (warm.sol==1) %have actual full solution
    Prob = WarmDefSOL('snopt',Prob,warm.result);
end

Prob = ProbCheck(Prob,'snopt');

Result = snoptTL(Prob);

converged = Result.Inform; %GOOD: {1,2,3}

%Compute trajectories
MP_state = zeros(size(L_e,2),n);
x_nom = zeros(N+1,n);
for i = 1:n
    c = Result.x_k(i:n:n*(N+1)-(n-i))';
    MP_state(:,i) = (c*L_e)';
    x_nom(:,i) = c';
end


MP_ctrl = zeros(size(L_e,2),m);
u_nom = zeros(N+1,m);
for j = 1:m
    c = Result.x_k(n*(N+1)+j:m:end-(m-j))';
    MP_ctrl(:,j) = (c*L_e)';
    u_nom(:,j) = c';
end

warm.result = Result;
warm.state = x_nom;
warm.ctrl = u_nom;

end