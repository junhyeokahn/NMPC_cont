function [u_aux,converged, E_rate] = compute_quad_aux(Prob,...
    X_q,Xq_dot,X_xi,E,M,f,B,u_nom,u_prev,eps_u,lambda)

%%

k = 2;

%A*aux <= u_b
A = 2*Xq_dot(:,k)'*(M(X_xi(:,k))*B);

l_b = -Inf;
u_b = -2*lambda*E + 2*Xq_dot(:,1)'*(M(X_xi(:,1))*(f(X_xi(:,1)) + B*u_nom)) -...
             2*Xq_dot(:,k)'*(M(X_xi(:,k))*(f(X_xi(:,k)) + B*u_nom));

% Prob = replace_A(Prob,A,l_b,u_b);
% Prob = modify_c(Prob,2*eps_u*(u_nom-u_prev));
% 
% Prob = ProbCheck(Prob, 'qpopt');
% Result = qpoptTL(Prob);
% u_aux = Result.x_k;
% converged = Result.Inform; %GOOD: 0

% a + b'aux <= 0
a = clean(-u_b,1e-7);
b = clean(A',1e-7);
if (norm(b)==0) || (a <= 0)
    u_aux = zeros(length(b),1);
else
    u_aux = -(a/(b'*b))*b;
end
converged = 0;

E_rate = 2*Xq_dot(:,k)'*(M(X_q(:,k))*(f(X_xi(:,k)) + B*u_nom + B*u_aux)) - ...
         2*Xq_dot(:,1)'*(M(X_q(:,1))*(f(X_xi(:,1)) + B*u_nom));
         

end