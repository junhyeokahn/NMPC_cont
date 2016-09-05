function [NMPC_state,NMPC_ctrl,x_nom,u_nom,converged] = ...
    compute_NMPC(Prob,act_p,mpc_state,state_constr,x_eq,u_eq,...
                 n,m,N,L_e,warm)

%adjust initial nominal state guess
for i = 1:n
    if abs(mpc_state(i)) > abs(state_constr(i));
        mpc_state(i) = 0.98*sign(mpc_state(i))*abs(state_constr(i));
    end
end

%Solution guess
if (~warm.sol)
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
else
    %Find: desired (real) time points to evaluate previous soln
    tau = (1/2)*(warm.Tp*warm.s_t+warm.Tp) + warm.shift;
    tau = tau(tau<=warm.Tp);
    N_prev = length(tau);
    
    %Find: evaluation (scaled) time points
    s_e_des = (2*tau - warm.Tp)/warm.Tp;
    
    %Compute Lagrange
    L_e_prev = compute_Lagrange(length(s_e_des)-1,...
           N,s_e_des,warm.s_t); %scaled nodes are the same
    
    x0 = zeros((N+1)*n,1);
    u0 = zeros((N+1)*m,1);
    for k = 1:N_prev
        x_prev = warm.state'*L_e_prev(:,k);
        for i = 1:n
            if abs(x_prev(i)) > abs(state_constr(i));
                x_prev(i) = 0.98*sign(x_prev(i))*abs(state_constr(i));
            end
        end
        x0(1+(k-1)*n:k*n) = x_prev;
        u0(1+(k-1)*m:k*m) = warm.ctrl'*L_e_prev(:,k);
    end
           
    In = eye(n);
    for i = 1:n
        x0 = x0 + kron([zeros(1,N_prev),...
          linspace(x_prev(i),x_eq(i),N+1-N_prev)], In(i,:))';
    end
    
    Im = eye(m);
    for j = 1:m
        u0 = u0 + kron([zeros(N_prev,1);
                     u_eq(j)*ones(N+1-N_prev,1)], Im(:,j));
    end    
end
Prob = modify_x_0(Prob,[x0;u0]);
Prob.user.x_act = act_p;

Prob = ProbCheck(Prob, 'snopt');
Result = snoptTL(Prob);

converged = Result.Inform; %GOOD: {1,2,3}

%Compute trajectories
NMPC_state = zeros(size(L_e,2),n);
x_nom = zeros(N+1,n);
for i = 1:n
    c = Result.x_k(i:n:n*(N+1)-(n-i))';
    NMPC_state(:,i) = (c*L_e)';
    x_nom(:,i) = c';
end

NMPC_ctrl = zeros(size(L_e,2),m);
u_nom = zeros(N+1,m);
for j = 1:m
    c = Result.x_k(n*(N+1)+j:m:end-(m-j))';
    NMPC_ctrl(:,j) = (c*L_e)';
    u_nom(:,j) = c';
end

end