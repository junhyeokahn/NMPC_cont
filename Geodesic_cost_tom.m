function J = Geodesic_cost_tom(vars,w,n,K,W_fnc,Phi,Phi_dot)

J = 0;

global GEO_X;
global GEO_MXDOT;

for k = 1:K+1 %0 ---> K
    x_k = Phi(:,:,k)*vars;
    x_dot_k = Phi_dot(:,:,k)*vars;
    
    W = W_fnc(x_k);
    M = W\eye(n);
    
    M_xdot = M*x_dot_k;
    J = J + (1/2)*w(k)*(x_dot_k'*M_xdot); 
    
    GEO_X(:,k) = x_k;
    GEO_MXDOT(:,k) = M_xdot;
end


return