clear all; close all; clc;
yalmip('clear');

warning('off','MATLAB:lang:badlyScopedReturnValue');

%% Constants

n = 5;
g = 9.81;
ccm_eps = 0.01;

%Dynamics constants

p_lim = pi/4;
pd_lim = pi/3;
vy_lim = 2;
vz_lim = 1.0;

%% Pick a solution

lambda = 1; 
condn = 132.8;
% lambda = 0.81;
% condn = 128.55;
return_metric = 1;

[sos_prob, w_lower, w_upper] = find_metric_PVTOL_SPOT_red(n,g,p_lim,pd_lim,vy_lim,vz_lim,...
                                condn,lambda,ccm_eps,return_metric);

%% Check CCM conditions red
load('metric_PVTOL_red_vectorized.mat');
disp('Checking CCM conditions and Computing control bound...');
lambda = 0.99*lambda;

dW_vz_fnc = @(x) zeros(n);
dW_pd_fnc = @(x) zeros(n);

m = 0.486;
J = 0.00383;
len = 0.25;

Bw = @(x)[zeros(1,2),cos(x(5)),-sin(x(5)),0]';
  
B_perp = @(x)[eye(2), zeros(2,1);
              zeros(3,2), [1;0;-x(4)]];

ctrl_N = 12;
p_range = linspace(-p_lim, p_lim, ctrl_N);
vy_range = linspace(-vy_lim, vy_lim, ctrl_N);
vz_range = linspace(-vz_lim, vz_lim, ctrl_N);

% sin_x = @(x) 0.7264*(x/(pi/4));
% cos_x = @(x) 0.8516 - 0.1464*(2*(x/(pi/4))^2 -1);

df_mat = @(x) [0,0,cos(x(5)),-sin(x(5)),-x(3)*sin(x(5))-x(4)*cos(x(5));
               0,0,sin(x(5)),cos(x(5)),x(3)*cos(x(5))-x(4)*sin(x(5));
               0,0,0,0,-g*cos(x(5));
               0,0,0,0, g*sin(x(5));
               zeros(1,5)];
f_mat = @(x) -g*sin(x(5)); %vy_dot
              
eig_CCM = zeros(ctrl_N, ctrl_N, ctrl_N);
eig_W = zeros(ctrl_N,ctrl_N,2);
sigma_ThBw = zeros(ctrl_N,ctrl_N,ctrl_N);

for i = 1:length(p_range)
    for j = 1:length(vy_range)
        for k = 1:length(vz_range)
            x = [0;0;vy_range(j);vz_range(k);p_range(i)];
            
            W = W_eval(w_poly_fnc(x));
            M = W\eye(n);
            Theta = chol(M);
            Theta_Bw = Theta*Bw(x);
            sigma_ThBw(i,j,k) = max(sqrt(eig(Theta_Bw'*Theta_Bw)));
            
            f = f_mat(x);
            df = df_mat(x);
            F = -W_eval(dw_poly_vy_fnc(x))*f + ...
                df*W + W*df' + 2*lambda*W;
            
            R_CCM = -B_perp(x)'*F*B_perp(x);
            
            eig_CCM(i,j,k) = min(eig(R_CCM));
            eig_W(i,j,1) = min(eig(W));
            eig_W(i,j,2) = max(eig(W));
        end
    end
end
d_bar = max(sigma_ThBw(:))/lambda;
disp('d_bar'); disp(d_bar);
% disp('euc_bound'); disp(d_bar*sqrt(w_upper));
disp('W:'); disp(min(min(eig_W(:,:,1))));
disp(max(max(eig_W(:,:,2))));
disp('CCM:'); disp(min(eig_CCM(:)));
disp(d_bar*sqrt(diag(W_upper)));







