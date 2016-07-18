function [u_aux,converged] = compute_opt_aux(Prob,K_e,X,X_dot,E,...
                            W,f,B,u_nom,lambda)


% edWfe = zeros(K_e+1,1);
% eWdfe = zeros(K_e+1,1);
% eWe = E*ones(K_e+1,1);
% eB = zeros(K_e+1,(K_e+1)*m);
% 
% for k = 1:K_e+1
%    
%     W_geo = W(X(:,k));
%     Eta = W_geo\X_dot(:,k);
%     
%     dWf = zeros(n,n);
%     dW_geo = dW(X(:,k));
%     f_geo = f(X(:,k));
%     for i = 1:n
%         dWf = dWf + dW_geo{i}*f_geo(i);
%     end
%     edWfe(k) = Eta'*dWf*Eta;
%     We = W_geo*Eta;
%     eWdfe(k) = (We')*df(X(:,k))'*Eta;
%     eB(k,1+(k-1)*m:k*m) = Eta'*B;
%     eWe(k) = Eta'*We;
% end

% Prob = replace_A(Prob,2*eB,-Inf*ones((K_e+1)*m,1),-2*lambda*eWe + edWfe - 2*eWdfe);
% 
% Prob = ProbCheck(Prob, 'qpopt');
% Result = qpoptTL(Prob);
% 
% delta_u = Result.x_k;
% 
% for k = 1:K_e+1
%     W_geo = W(X(:,k));
%     Eta = W_geo\X_dot(:,k);
%     if norm(Eta'*B)<1e-3
% %         disp('found');
%         delta_u(1+(k-1)*m:k*m) = zeros(m,1);
%     end
% end
% u_aux = (1/2)*Q*delta_u;


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



