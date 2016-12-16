function [u_aux,converged] = compute_opt_aux(Prob,K_e,X,X_dot,E,...
                            W,f,B,u_nom,lambda)

%%

w_poly = W.w_poly_fnc(X);

k = K_e + 1;
A = 2*X_dot(:,k)'*(W.W_eval(w_poly(:,k))\B);

% l_b = -Inf;
u_b = -2*lambda*E + 2*X_dot(:,1)'*(W.W_eval(w_poly(:,1))\(f(X(:,1)) + B*u_nom)) -...
                    2*X_dot(:,k)'*(W.W_eval(w_poly(:,k))\(f(X(:,k)) + B*u_nom));

% Prob = replace_A(Prob,A,l_b,u_b);
% Prob = modify_c(Prob,2*eps_u*(u_nom-u_prev));
% 
% Prob = ProbCheck(Prob, 'qpopt');
% Result = qpoptTL(Prob);
% u_aux = (Result.x_k)';

a = -u_b;
b = A';
if (norm(b)==0) || (a <= 0)
    u_aux = zeros(1,length(b));
else
    u_aux = -(a/(b'*b))*b';
end


%%
% converged = Result.Inform; %GOOD: 0

converged = 0;
end



