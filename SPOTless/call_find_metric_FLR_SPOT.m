clear all; close all; clc;
yalmip('clear');

warning('off','MATLAB:lang:badlyScopedReturnValue');

%% Function approximations

% sin_x =  @(x) 1.1336*(x/(pi/2)) -0.1381*(4*(x/(pi/2))^3 - 3*(x/(pi/2)));
% cos_x = @(x) 0.472 - 0.4994*(2*(x/(pi/2))^2 -1) + ...
%              0.028*(8*(x/(pi/2))^4 - 8*(x/(pi/2))^2 + 1);

%% Constants

n = 4;
% lambda = 1;
g = 9.81;
ccm_eps = 1e-1;

%Dynamics constants
m = 1;
l = 1; %L = 2*l
I = (1/3)*m*(2*l)^2;
sigma = 100;
J = 1;
b = 1;

x1_lim = pi; %link angle
x2_lim = 5; %link rate
x3_lim = pi; %rotor angle
% x4_lim = pi; %rotor rate

%% Define problem

% lambda_range = linspace(0.1,4.0,13);
% lambda_range = (1/100)*round(lambda_range*100);
% lambda_range = 2.5;
% euc_bounds = NaN(length(lambda_range),1);
% d_bars = NaN(length(lambda_range),1);
% cond_bound = NaN(length(lambda_range),1);
% 
% eps = 1;
% condn_prev = 350;
% return_metric = 0;
% 
% for ll = 1:length(lambda_range)
%     lambda = lambda_range(ll);
%     
%     fprintf('**********\n');
%     fprintf('lambda: %f\n', lambda);
%     solved = 0;
%     
%     %Determine upper bound
%     cond_l = condn_prev;
%     cond_u = 1.2*condn_prev;
%     while (~solved) 
%         fprintf(' cond_u: %.2f: ', cond_u);
%         [sos_prob,~,~,~,~] = find_metric_FLR(n,m,g,l,I,J,b,sigma,x1_lim,x2_lim,x3_lim,...
%                                 cond_u,lambda,ccm_eps,return_metric);
%         if (sos_prob == 0)
%             solved = 1;
%             fprintf('feasible \n');
%         else
%             %shift up condition number range
%             fprintf('\n');
%             cond_l = cond_u;
%             cond_u = 1.2*cond_u;
%         end
%     end
%     if (solved)
%         euc_bounds(ll) = sqrt(cond_u)/lambda;
%         fprintf(' cond_l: %.2f, cond_u: %.2f\n', cond_l, cond_u);
%     else
%         continue;
%     end
%     
%     %Now do bisection search
%     while(cond_u - cond_l >= eps)
%         condn = (cond_l+cond_u)/2;
%         fprintf(' cond: %.2f', condn);
%         
%         [sos_prob, w_lower, w_upper,~,~] = find_metric_FLR(n,m,g,l,I,J,b,sigma,x1_lim,x2_lim,x3_lim,...
%                                 condn,lambda,ccm_eps,return_metric);
%         
%         if (sos_prob == 0)
%             fprintf(' feasible\n');
%             
%             euc_bounds(ll) = sqrt(double(w_upper/w_lower))/lambda;
%             d_bars(ll) = sqrt(double(1/w_lower))/lambda;           
% 
%             cond_u = condn;
%         else
%             fprintf(' infeasible\n');
%             cond_l = condn;
%         end
%     end
%     condn_prev = cond_u;
%     cond_bound(ll) = cond_u;
%     disp('Euc_bound:'); disp(euc_bounds(ll));
%     disp('d_bar:'); disp(d_bars(ll));
%     fprintf('**********\n');
%     
% end
% 
% pause;

%% Pick a solution

% load metric_FLR_comp_final.mat;
% lambda = lambda_range(11); 
% condn = cond_bound(11);
lambda = 2.5; 
condn = 379;
return_metric = 1;

[sos_prob, w_lower, w_upper,w_poly_fnc, dw_poly_x1_fnc, dw_poly_x2_fnc, W_eval, W_upper] = find_metric_FLR_SPOT(n,m,g,l,I,J,b,sigma,x1_lim,x2_lim,x3_lim,...
                                condn,lambda,ccm_eps,return_metric);

save('metric_FLR_vectorized.mat','W_eval','w_poly_fnc','dw_poly_x1_fnc','W_upper');

%% Compute aux control bound

disp('Checking CCM conditions and Computing control bound...');

B = [zeros(3,1);1/J];
Bw = [0,(1/I),0,0;
      0,0,0,(1/J)]';
B_perp = [eye(3); zeros(1,3)];

% d_bar = sqrt(double(1/w_lower))/lambda; %normalized

ctrl_N = 15;
x1_range = linspace(-x1_lim, x1_lim, ctrl_N);
x2_range = linspace(-x2_lim, x2_lim, ctrl_N);
x3_range = linspace(-x3_lim, x3_lim, ctrl_N);

df_mat = @(x) [0, 1, 0, 0;
    (m*g*l/I)*cos(x(1))-(sigma/I), 0, (sigma/I), 0;
    zeros(1,3), 1;
    (sigma/J), 0, -(sigma/J), -(b/J)];
f_mat = @(x) [x(2);
    (m*g*l/I)*sin(x(1)) - (sigma/I)*(x(1)-x(3));
    x(4);
    (sigma/J)*(x(1)-x(3)) - (b/J)*x(4)];

delta_u = zeros(ctrl_N,ctrl_N,ctrl_N);
eig_CCM = zeros(ctrl_N, ctrl_N, ctrl_N);
eig_W = zeros(ctrl_N,ctrl_N,2);
sigma_ThBw = zeros(ctrl_N,ctrl_N);

for i = 1:length(x2_range)
    for j = 1:length(x1_range)
        for k = 1:length(x3_range)
            x = [x1_range(j); x2_range(i); x3_range(k); 0];
            
            W = W_eval(w_poly_fnc(x));
            M = W\eye(4);
            Theta = chol(M);
            Theta_Bw = Theta*Bw;
            sigma_ThBw(i,j) = max(sqrt(eig(Theta_Bw'*Theta_Bw)));
            
            L = chol(W);
            
            f = f_mat(x);
            df = df_mat(x);
            F = -W_eval(dw_poly_x1_fnc(x))*f(1) + ...
                df*W + W*df' + 2*lambda*W;
            
            delta_u_den = eig((inv(L))'*(B*B')*inv(L));
            delta_u(i,j,k) = 0.5*max(eig((inv(L))'*F*inv(L)))/...
                    sqrt(min(delta_u_den(delta_u_den>0)));
            
            R_CCM = -B_perp'*F*B_perp;
            
            eig_CCM(i,j,k) = min(eig(R_CCM));
            eig_W(i,j,1) = min(eig(W));
            eig_W(i,j,2) = max(eig(W));
        end
    end
end
d_bar = max(sigma_ThBw(:))/lambda;
disp('d_bar'); disp(d_bar);
disp('euc_bound'); disp(d_bar*sqrt(w_upper));
disp('Control:'); disp(max(d_bar*delta_u(:)));
disp('W:'); disp(min(min(eig_W(:,:,1))));
disp(max(max(eig_W(:,:,2))));
disp('CCM:'); disp(min(eig_CCM(:)));


%% CVX approx method

% f = [x2;
%         (m*g*l/I)*sin_x1 - (sigma/I)*(x1-x3);
%         x4;
%         (sigma/J)*(x1-x3) - (b/J)*x4];

% phi = @(x) [x(1);
%             x(2);
%            (m*g*l/I)*sin(x(1)) - (sigma/I)*(x(1)-x(3));
%            (m*g*l/I)*x(2)*cos(x(1)) - (sigma/I)*(x(2)-x(4))];
%
% phi_d = @(x) [1, zeros(1,3);
%               0, 1, 0, 0;
%             (m*g*l/I)*cos(x(1))-(sigma/I), 0, (sigma/I), 0;
%              -(m*g*l/I)*x(2)*sin(x(1)), (m*g*l/I)*cos(x(1))-(sigma/I), 0, (sigma/I)];

%
% clear all; close all; clc;
%
%
%
% clear W
%
% Ac = [zeros(3,1), eye(3); zeros(1,4)];
% Bc = [zeros(3,1); 1];
% B_perp = [eye(3);
%     zeros(1,3)];
%
% % lambda_range = (1/100)*round(100*linspace(0.5,2.5,11));
% lambda_range = 2;
% euc_bounds = NaN(length(lambda_range),1);
% d_bars = NaN(length(lambda_range),1);
%
% for ll = 1:length(lambda_range)
%
%     lambda = lambda_range(ll);
%     disp('lambda:'); disp(lambda);
%
%     cvx_begin sdp
%     variable W(4,4) symmetric
%     variables w_lower w_upper
%     minimize(w_upper - w_lower)
%     subject to
%
%     W >= w_lower*eye(4);
%     W <= w_upper*eye(4);
%     w_lower >= 1;
%
%     B_perp'*(Ac*W + W*Ac' + 2*lambda*W)*B_perp <= 0;
%
%     cvx_end
%
%     if strcmp(cvx_status, 'Solved') == 1
%         solved = 1;
%         fprintf(' feasible\n');
%
%         euc_bounds(ll) = sqrt(double(w_upper/w_lower))/lambda;
%         disp('Euc_bound:'); disp(euc_bounds(ll));
%
%         d_bars(ll) = sqrt(double(1/w_lower))/lambda;
%         disp('d_bar:'); disp(d_bars(ll));
%     else
%         fprintf(' infeasible\n');
%     end
% %     clear w_lower w_upper W
%
% end






