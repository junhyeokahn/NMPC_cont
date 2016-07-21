function [u_aux,converged] = compute_opt_aux(Prob,K_e,X,X_dot,E,...
                            W,f,B,u_nom,lambda)

%%

k = K_e + 1;
A = 2*X_dot(:,k)'*(W(X(:,k))\B);

l_b = -Inf;
u_b = -2*lambda*E + 2*X_dot(:,1)'*(W(X(:,1))\(f(X(:,1)) + B*u_nom)) -...
             2*X_dot(:,k)'*(W(X(:,k))\(f(X(:,k)) + B*u_nom));

Prob = replace_A(Prob,A,l_b,u_b);

Prob = ProbCheck(Prob, 'qpopt');
Result = qpoptTL(Prob);
u_aux = Result.x_k;

%%
converged = Result.Inform; %GOOD: 0

end



