clear all; close all; clc;
yalmip('clear');

%% Constants

x1_lim = 5;
x2_lim = 5;
ccm_eps = 0;

%% Define problem & combine line + bisection searches

% lambda_range = 1.742857142857143; %final optimal
% lambda_range = linspace(0.8,2,15); %line search range
% euc_bounds = NaN(length(lambda_range),1);
% d_bars = NaN(length(lambda_range),1);
% cond_bound = NaN(length(lambda_range),1);
% 
% eps = 1e-3;
% condn_prev = 1.63;
% return_metric = 0;

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
%         [sos_prob,~,~] = find_metric_Allgower_SPOT(x1_lim, x2_lim,cond_u, lambda, ccm_eps,return_metric);
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
%         fprintf(' cond_l: %.4f, cond_u: %.4f\n', cond_l, cond_u);
%     else
%         continue;
%     end
%     
%     %Now do bisection search
%     while(cond_u - cond_l >= eps)
%         condn = (cond_l+cond_u)/2;
%         fprintf(' cond: %.4f', condn);
%         
%         [sos_prob, w_lower, w_upper] = find_metric_Allgower_SPOT(x1_lim, x2_lim,condn, lambda, ccm_eps,return_metric);
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

% figure()
% plot(lambda_range, euc_bounds,'ro','markerfacecolor','g','markersize',20);
% grid on
% xlabel('\lambda');
% ylabel('$\|x^{*}-x\|/\bar{w}$','Interpreter','Latex');
% % title('Robustness optimization');
% set(findall(gcf,'type','text'),'FontSize',32);set(gca,'FontSize',32)

% pause;

%% Pick a solution

lambda = 1.742857142857143; 
condn = 1.6325;
return_metric = 1;

[sos_prob, w_lower, w_upper] = find_metric_Allgower_SPOT(x1_lim, x2_lim,condn, lambda, ccm_eps,return_metric);

%% Compute control bounds using optimal metric (W = const matrix)
load('metric_Allgower.mat');

M = W_sol\eye(2);
alpha = max(eig(M));
w = 0.1;
d_bar = sqrt(alpha)*w/lambda;

%Compare RCI sets
P_rci = diag([39.0251, 486.0402]);

figure()
Ellipse_plot(M*(1/d_bar^2),[0;0],20,'k');
Ellipse_plot(P_rci,[0;0],20,'r');

L = chol(W_sol);
df_mat = @(x) [-1, 2;
     -3, 4-0.75*x(2)^2];
B = [0.5;-2];
F = @(x) df_mat(x)*W_sol + W_sol*df_mat(x)' + 2*lambda*W_sol;

x2_range = linspace(-x2_lim,x2_lim,30);
delta_u = zeros(length(x2_range),1);

for i = 1:length(x2_range)
    F_sol = F([0;x2_range(i)]);
    
    delta_u(i) = 0.5*d_bar*max(eig((inv(L))'*F_sol*inv(L)))/...
                          sqrt(max(eig((inv(L))'*(B*B')*inv(L))));
end

figure()
plot(x2_range,delta_u,'b-','linewidth',2);
grid on
xlabel('x2'); ylabel('$\bar{\delta}_u$','interpreter','latex');
set(findall(gcf,'type','text'),'FontSize',28);set(gca,'FontSize',28)


