clear all; close all; clc;
% yalmip('clear');s

warning('off','MATLAB:lang:badlyScopedReturnValue');

%% Constants

n = 9;
% lambda = 1;
g = 9.81;
ccm_eps = 1e-1;

r_lim = pi/3;
p_lim = pi/3;
th_lim_low = 0.5*g;
th_lim_high = 2*g;

%% Pullback method

lambda = 0.8;

B_perp = [eye(6); zeros(3,6)];
Ac = [zeros(3), eye(3), zeros(3);
      zeros(3), zeros(3), eye(3);
      zeros(3,9)];
  
W_scale = blkdiag(diag([0.2,0.2,0.2]),diag([1,1,1]),...
                    1*eye(3));  
cond_l = 1; cond_u = 300;
eps = 0.5;
while(cond_u - cond_l > eps)
    condn = 0.5*(cond_l+cond_u);
    fprintf('condn: %f:', condn);
    
    cvx_begin sdp quiet
    variable W_xi(9,9) symmetric
    variables w_lower w_upper
    minimize(0)
    subject to
    w_lower >= 1;
    condn*w_lower >= trace(W_scale*W_xi);
    W_xi >= w_lower*eye(9);
    W_xi <= w_upper*eye(9);
    B_perp'*(Ac*W_xi + W_xi*Ac')*B_perp <= -2*lambda*(B_perp'*(W_xi)*B_perp);
    cvx_end
    
    if strcmp(cvx_status,'Solved')
        fprintf('feasible \n');
        cond_u = condn;
        W_best = W_xi;
    else
        fprintf('infeasible \n');
        cond_l = condn;
    end
end
W_xi = W_best;
M_xi = W_xi\eye(n);
M_xi = 0.5*(M_xi + M_xi');

b_T = @(x) [sin(x(9)); -cos(x(9))*sin(x(8)); cos(x(9))*cos(x(8))];

db_T_q =@(x) [0, cos(x(9));
    -cos(x(8))*cos(x(9)), sin(x(8))*sin(x(9));
    -sin(x(8))*cos(x(9)),-cos(x(8))*sin(x(9))];

Phi = @(x) blkdiag(eye(3),eye(3),-[b_T(x), db_T_q(x)*x(7)]);

M_pull = @(x) Phi(x)'*M_xi*Phi(x);

keyboard;

%% Test pullback metric

Bw = [zeros(3);
    eye(3);
    zeros(3)];
ctrl_N = 10;
r_range = linspace(-r_lim, r_lim, ctrl_N);
p_range = linspace(-p_lim, p_lim, ctrl_N);
th_range = linspace(th_lim_low, th_lim_high, 2*ctrl_N);

eig_CCM = zeros(ctrl_N, ctrl_N, ctrl_N);
eig_M = zeros(ctrl_N,ctrl_N, ctrl_N,2);
sigma_ThBw = zeros(ctrl_N,ctrl_N,ctrl_N);

for i4 = 1:length(r_range)
    for i5 = 1:length(p_range)
        for i6 = 1:length(th_range)
            
            x = [zeros(6,1);
                 th_range(i6); r_range(i4); p_range(i5)];
            
            M = M_pull(x);
            Theta = chol(M);
            Theta_Bw = Theta*Bw;
            sigma_ThBw(i4,i5,i6) = max(sqrt(eig(Theta_Bw'*Theta_Bw)));
            
            eig_M(i4,i5,i6,1) = min(eig(M));
            eig_M(i4,i5,i6,2) = max(eig(M));
        end
    end
end
d_bar = 0.1*max(sigma_ThBw(:))/lambda;
disp('d_bar'); disp(d_bar);
disp('M:'); disp(min(min(min(eig_M(:,:,:,1)))));
disp(max(max(max(eig_M(:,:,:,2)))));

%% Compute Pullback bound

[m_lower, M_lower_pull] = compute_QUAD_bound_pull(M_xi,r_lim,p_lim,th_lim_low,th_lim_high);
W_upper_pull = M_lower_pull\eye(n);

disp('euc_bounds (x)');
disp(d_bar*sqrt(diag(W_upper_pull)));

disp('euc_bounds (xi)');
disp(d_bar*sqrt(diag(W_xi)));


%% Save

save('metric_QUAD_pullback.mat','d_bar','M_xi','M_lower_pull','lambda');