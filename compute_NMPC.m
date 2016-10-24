function [NMPC_state,NMPC_ctrl,converged,warm,Prob] = ...
    compute_NMPC(Prob,act_p,~,state_constr,ctrl_constr,MP_state,MP_ctrl,...
    n,m,N,L_e,warm,dt)

%Solution guess
if (~warm.sol)
    
    x_term = MP_state((warm.Tp/dt)+1,:)';
    
    tau = (1/2)*(warm.Tp*warm.s_t+warm.Tp) ;
    state_prev = interp1([0:dt:warm.Tp]',MP_state(1:(warm.Tp/dt)+1,:),tau);
    ctrl_prev = interp1([0:dt:warm.Tp]',MP_ctrl(1:(warm.Tp/dt)+1,:),tau);
    
    x0 = zeros((N+1)*n,1);
    u0 = zeros((N+1)*m,1);
    for k = 1:N+1
        x_prev = state_prev(k,:)';
        for i = 1:n
            if abs(x_prev(i)) >= abs(state_constr(i));
                x_prev(i) = 0.98*sign(x_prev(i))*abs(state_constr(i));
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
    
    Prob = modify_x_0(Prob,[x0;u0]);
    
else
    %Find: desired (real) time points to evaluate previous soln
%     tau = (1/2)*(warm.Tp*warm.s_t+warm.Tp) ;
%     tau_int = tau+warm.shift;
%     tau_int = tau_int(tau_int <= warm.Tp);
%     state_prev = interp1(tau,warm.state,tau_int);
%     ctrl_prev = interp1(tau,warm.ctrl,tau_int);
%     N_prev = size(state_prev,1);
%     
%     x0 = zeros(N_prev*n,1);
%     u0 = zeros(N_prev*m,1);
%     for k = 1:N_prev
%         x_prev = state_prev(k,:)';
%         for i = 1:n
%             if abs(x_prev(i)) >= abs(state_constr(i));
%                 x_prev(i) = 0.98*sign(x_prev(i))*abs(state_constr(i));
%             end
%         end
%         x0(1+(k-1)*n:k*n) = x_prev;
%         
%         u_prev = ctrl_prev(k,:)';
%         for j = 1:m
%             if (u_prev(j) <= ctrl_constr(j,1))
%                 u_prev(j) = ctrl_constr(j,1);
%             elseif (u_prev(j) >= ctrl_constr(j,2))
%                 u_prev(j) = ctrl_constr(j,2);
%             end
%         end
%         u0(1+(k-1)*m:k*m) = u_prev;
%     end
%     %     x_term = warm.result.Prob.user.x_eq; %only used for first iteration
%     
%     %now forward shift
%     if (N+1-N_prev > 0)
%         tau = warm.solve_t + warm.shift + tau'; %[t_{i+1}, t_{i+1}+T]
%         tau = tau(tau > warm.solve_t + warm.Tp); %[t_i+T, t_{i+1}+T]
%         t_span = warm.solve_t + warm.Tp:dt:warm.solve_t+warm.shift+warm.Tp;
%         
%         i_start = round((warm.solve_t + warm.Tp)/dt) + 1;
%         i_end =  round((warm.solve_t + warm.shift + warm.Tp)/dt) + 1;
%         
%         if (N+1-N_prev > 1)
%             x_tail = interp1(t_span',MP_state(i_start:i_end,:),tau);
%             u_tail = interp1(t_span',MP_ctrl(i_start:i_end,:),tau);
%         else
%             x_tail = MP_state(i_start,:);
%             u_tail = MP_ctrl(i_start,:);
%         end
%         
%         for k = 1:N+1-N_prev
%             for i = 1:n
%                 if abs(x_tail(k,i)) >= abs(state_constr(i));
%                     x_tail(k,i) = 0.98*sign(x_tail(k,i))*abs(state_constr(i));
%                 end
%             end
%             for j = 1:m
%                 if (u_tail(k,j) <= ctrl_constr(j,1))
%                     u_tail(k,j) = ctrl_constr(j,1);
%                 elseif (u_tail(k,j) >= ctrl_constr(j,2))
%                     u_tail(k,j) = ctrl_constr(j,2);
%                 end
%             end
%         end
%         
%         %         x_term = x_tail(end,:)';
%         
%         x0 = [x0;reshape(x_tail',(N+1-N_prev)*n,1)];
%         u0 = [u0;reshape(u_tail',(N+1-N_prev)*m,1)];
%     end
    i_end =  round((warm.solve_t + warm.shift + warm.Tp)/dt) + 1;
    if (i_end > length(MP_state))
        i_end = length(MP_state);
    end
    x_term = MP_state(i_end,:)';
end

%Update constraint information
Prob.user.x_act = act_p;
Prob.user.x_eq = x_term;

%Recall warm solution
if (warm.sol)
      Prob = WarmDefSOL('snopt',Prob,warm.result);
end

if ~Prob.CHECK
    Prob = ProbCheck(Prob,'snopt');
end

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

warm.result = Result;
warm.state = x_nom;
warm.ctrl = u_nom;


end