function GradObj = Geodesic_grad(~,w,K,N,n,Ti,~,dW,~,Phi_dot,~)

global GEO_X;
global GEO_MXDOT;

GradObj = zeros(n*(N+1),1);

for k = 1:K+1 %0 ---> K
%     x_k = Phi(:,:,k)*vars;
%     x_dot_k = Phi_dot(:,:,k)*vars;
%     
%     W = W_fnc(x_k);
%     M = W\eye(n);
%     
%     M_xdot = M*x_dot_k;

    M_xdot = GEO_MXDOT(:,k);
        
%     W_dx = dW(GEO_X(:,k));
    
    GradObj = GradObj + w(k)*Phi_dot(:,:,k)'*M_xdot;
    
    GradObj = GradObj+...
           -(1/2)*w(k)*(M_xdot'*dW.p(GEO_X(:,k))*M_xdot)*Ti(:,k,1) + ...
           -(1/2)*w(k)*(M_xdot'*dW.vy(GEO_X(:,k))*M_xdot)*Ti(:,k,2);
    
%     for j = 1:length(n_W)
%         i = n_W(j);
%         GradObj = GradObj+...
%            -(1/2)*w(k)*(M_xdot'*W_dx{i}*M_xdot)*Ti(:,k,j);
%     end
end


return