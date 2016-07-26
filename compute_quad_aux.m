function [u_aux,converged] = compute_quad_aux(Prob,...
    X_q,Xq_dot,X_xi,E,M,f,B,u_nom,u_prev,eps_u,lambda)

%%

k = 2;
A = 2*Xq_dot(:,k)'*(M(X_q(:,k))*B);

l_b = -Inf;
u_b = -2*lambda*E + 2*Xq_dot(:,1)'*(M(X_q(:,1))*(f(X_xi(:,1)) + B*u_nom)) -...
             2*Xq_dot(:,k)'*(M(X_q(:,k))*(f(X_xi(:,k)) + B*u_nom));

Prob = replace_A(Prob,A,l_b,u_b);
Prob = modify_c(Prob,2*eps_u*(u_nom-u_prev));

Prob = ProbCheck(Prob, 'qpopt');
Result = qpoptTL(Prob);
u_aux = Result.x_k;
converged = Result.Inform; %GOOD: 0

% a = -u_b;
% b = A';
% if (norm(b)==0) || (a <= 0)
%     u_aux = zeros(length(b),1);
% else
%     u_aux = -(a/(b'*b))*b;
% end
% converged = 0;

end