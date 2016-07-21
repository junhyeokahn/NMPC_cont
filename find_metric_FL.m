function [d_bar,euc_bound,cntrl_bound,lambda_opt] = find_metric_FL(m,r_vec,lambda_range)

%m: number of inputs = number of outputs
%r_vec: order of derivative when input appears

%size of Ac: sum(m(i)*r_vec(i))
Ac = [];
Bc = [];
for i = 1:m
    Ac_i = [zeros(r_vec(i)-1,1), eye(r_vec(i)-1);
            zeros(1,r_vec(i))];
    Ac = blkdiag(Ac,Ac_i);
    
    Bc_i = [zeros(r_vec(i)-1,1);1];
    Bc = blkdiag(Bc,Bc_i);
end
B_perp = null(Bc');

n = size(Ac,1);

euc_bound = NaN(length(lambda_range),1);
d_bar = euc_bound;
cntrl_bound = d_bar;

for i = 1:length(lambda_range)
    lambda = lambda_range(i);
    cvx_begin sdp quiet
    variable W_ccm(n,n) symmetric
    variables w_lower w_upper
    minimize (w_upper - w_lower)
    subject to
    W_ccm >= w_lower*eye(n);
    W_ccm <= w_upper*eye(n);
    w_lower >= 0.01;

    B_perp'*(Ac*W_ccm + W_ccm*Ac' + 2*lambda*W_ccm)*B_perp <= 0;
    cvx_end
    if strcmp(cvx_status,'Solved')==1
        condn = w_upper/w_lower;
        euc_bound(i) = sqrt(condn)/lambda;
        d_bar(i) = sqrt(1/w_lower)/lambda;
        L = chol(W_ccm);
        F = Ac*W_ccm + W_ccm*Ac' + 2*lambda*W_ccm;
        v = (eig((Bc'*inv(L))'*(Bc'*inv(L))));
        cntrl_bound(i) = 0.5*d_bar(i)*max(eig((inv(L))'*F*(inv(L))))/...
                          sqrt(min(v(v>0)));
    end
end

[~, pos] = min(euc_bound);
lambda_opt = lambda_range(pos);

end

