function [u_aux,converged] = compute_opt_aux_FL(Prob,K_e,X,X_dot,E,...
                            M,f,B,u_nom,u_prev,eps_u,lambda)

%%

k = K_e + 1;
A = 2*X_dot(:,k)'*(M(X(:,k))*B);

l_b = -Inf;
u_b = -2*lambda*E + 2*X_dot(:,1)'*(M(X(:,1))*(f(X(:,1)) + B*u_nom)) -...
             2*X_dot(:,k)'*(M(X(:,k))*(f(X(:,k)) + B*u_nom));

% Prob = replace_A(Prob,A,l_b,u_b);
% Prob = modify_c(Prob,2*eps_u*(u_nom-u_prev));
% 
% Prob = ProbCheck(Prob, 'qpopt');
% Result = qpoptTL(Prob);
% u_aux = Result.x_k;
% converged = Result.Inform; %GOOD: 0

a = clean(-u_b,1e-3);
b = clean(A',1e-3);
if (norm(b)==0) || (a <= 0)
    u_aux = zeros(length(b),1);
else
    u_aux = -(a/(b'*b))*b;
end
converged = 0;

end