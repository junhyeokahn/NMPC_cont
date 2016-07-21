clear all; close all; clc;
warning('off','MATLAB:nargchk:deprecated');
%% Constants
n = 2;
m = 1;

%% Setup NMPC Problem

f  = @(x) [-1*x(1) + 2*x(2);
           -3*x(1) + 4*x(2) - 0.25*(x(2)^3)];
B = [0.5;-2];

df = @(x) [-1, 2;
           -3, 4-0.75*(x(2)^2)];
B_w = [0;1];

w_max = 0.1;

W = @(x) [4.258279173109496,  -0.934234854771844;
        -0.934234854771844,   3.766923169705589];

M_ccm = W(0)\eye(2);
lambda =  1.742857142857143;
d_bar = (w_max*sqrt(max(eig(M_ccm)))/lambda);

N_mpc = 50;

P = [7.9997, -12.2019;
    -12.2019, 27.0777];
alpha = 10;

Tp = 1.5;
delta = 0.1;
dt = 0.005;

state_constr_low = [-4.94;-4.94]; 
ctrl_constr_low_CCM = -1.793*ones(m,1);

x_eq = [0;0];

[NMPC_Prob_CCM,~,~] = setup_NMPC(n,m,...
    f,B,df, state_constr_low,ctrl_constr_low_CCM,...
    N_mpc,Tp,delta,dt,...
    P,alpha,M_ccm,d_bar^2,...
    x_eq);

M_alg = diag([39.0251, 486.0402]);
d_bar_alg = 1;
ctrl_constr_low_alg = -1.15*ones(m,1);

[NMPC_Prob_alg,L_e,L_e_full] = setup_NMPC(n,m,...
    f,B,df, state_constr_low,ctrl_constr_low_alg,...
    N_mpc,Tp,delta,dt,...
    P,alpha,M_alg,d_bar_alg^2,...
    x_eq);

%% Get feasibility domains
x1_range = linspace(-4.5,4.5,100);
x2_range = x1_range;

TOTAL_COUNT = length(x1_range)*length(x2_range);

[X1,X2] = meshgrid(x1_range, x2_range);

feas_CCM = zeros(size(X1));
feas_allg = feas_CCM;
count = 1;

for i = 1:length(x2_range)
    for j = 1:length(x1_range)
        
        test_state = [x1_range(j);x2_range(i)];
        fprintf('%d/%d \n',count,TOTAL_COUNT);
        
        [~,~,converged_CCM] = compute_NMPC(NMPC_Prob_CCM,...
            test_state,test_state,state_constr_low,x_eq,...
            n,m,N_mpc,L_e_full);
        
        if (converged_CCM==1) || (converged_CCM==2) || (converged_CCM==3)
            feas_CCM(i,j) = 1;
        end
        
        [~,~,converged_allg] = compute_NMPC(NMPC_Prob_alg,...
            test_state,test_state,state_constr_low,x_eq,...
            n,m,N_mpc,L_e_full);
        
        if (converged_allg==1) || (converged_allg==2) || (converged_allg==3)
            feas_allg(i,j) = 1;
        end
        
        count = count+ 1;
    end
end

%% Plot

figure()
subplot(1,2,1)
surface(X1,X2,feas_CCM);
%shading('interp');
grid on
xlabel('x_1'); ylabel('x_2');
set(findall(gcf,'type','text'),'FontSize',32);set(gca,'FontSize',32)

subplot(1,2,2)
surface(X1,X2,feas_allg);
%shading('interp');
grid on
xlabel('x_1'); ylabel('x_2');
set(findall(gcf,'type','text'),'FontSize',32);set(gca,'FontSize',32)



