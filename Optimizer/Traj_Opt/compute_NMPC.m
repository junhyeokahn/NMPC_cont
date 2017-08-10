function [NMPC_state,NMPC_ctrl,T_end,converged,warm,Prob] = ...
    compute_NMPC(Prob,act_p,state_constr,ctrl_constr,MP_state,MP_ctrl,...
    n,m,N,L_e,warm,dt,E_s)

%% Solution guess
if (~warm.sol) %don't have guess for first MPC iteration
    
    x_term = MP_state((warm.Tp/dt)+1,:)';
    
    tau = (1/2)*(warm.Tp*warm.s_t+warm.Tp) ;
    state_prev = interp1([0:dt:warm.Tp]',MP_state(1:(warm.Tp/dt)+1,:),tau);
    ctrl_prev = interp1([0:dt:warm.Tp]',MP_ctrl(1:(warm.Tp/dt)+1,:),tau);
    
    x0 = zeros((N+1)*n,1);
    u0 = zeros((N+1)*m,1);
    for k = 1:N+1
        x_prev = state_prev(k,:)';
        for i = 1:n
            if x_prev(i) >= state_constr(i,2)
                x_prev(i) = state_constr(i,2)-0.01*(state_constr(i,2)-state_constr(i,1));
            elseif x_prev(i) < state_constr(i,1)
                x_prev(i) = state_constr(i,1)+0.01*(state_constr(i,2)-state_constr(i,1));
            end
        end
        
        x0(1+(k-1)*n:k*n) = x_prev;
        
        u_prev = ctrl_prev(k,:)';
        for j = 1:m
            if (u_prev(j) <= ctrl_constr(j,1))
                u_prev(j) = ctrl_constr(j,1);
            elseif (u_prev(j) >= ctrl_constr(j,2))
                u_prev(j) = ctrl_constr(j,2);
            end
        end
        u0(1+(k-1)*m:k*m) = u_prev;
    end
    
    Prob = modify_x_0(Prob,[x0;u0;warm.Tp]);
    
end

%% Update constraint information
Prob.user.x_act = act_p;

Prob = modify_c_U(Prob,E_s,n*(N+1)+1);
Prob = modify_x_L(Prob,min(warm.solve_t + warm.shift + warm.Tp,Prob.user.t_nom(end)),(n+m)*(N+1)+1);

%% Update scaled obstacles
obs_mpc = Prob.user.obs;
tau = (1/2)*(warm.Tp*warm.s_t+warm.Tp);
E_time_bound = (sqrt(E_s)*exp(-warm.lambda*tau) + warm.d_bar*(1-exp(-warm.lambda*tau))).^2;
for k = 1:N+1
    for i = 1:obs_mpc.n_obs
        S_new = (sqrt(E_time_bound(k)*(obs_mpc.S\eye(2))) + (obs_mpc.r(i))*eye(2))^2\eye(2);
        obs_mpc.M{k}(:,:,i) = obs_mpc.U*S_new*obs_mpc.V';
    end
end
Prob.user.obs = obs_mpc;

%% Recall warm solution
if (warm.sol)
      Prob = WarmDefSOL('snopt',Prob,warm.result);
end

if ~Prob.CHECK
    Prob = ProbCheck(Prob,'snopt');
end

%% Solve
Result = snoptTL(Prob);

converged = Result.Inform; %GOOD: {1,2,3}

%% Compute trajectories
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
    c = Result.x_k(n*(N+1)+j:m:end-1-(m-j))';
    NMPC_ctrl(:,j) = (c*L_e)';
    u_nom(:,j) = c';
end

T_end = Result.x_k(end);

warm.result = Result;
warm.state = x_nom;
warm.ctrl = u_nom;


end