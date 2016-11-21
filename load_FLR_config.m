%% Constants

n = 4;
m = 1;

%% Obstacle info

obs = struct('n_obs',0);
obs_mpc = struct('n_obs',0);

%% Setup Metric

load 'metric_FLR.mat';

W_fnc = @(x) W_mat(wrapToPi(x(1)));
dW_fnc = @(x) {dW_x1_mat(wrapToPi(x(1)))};

sigma_ThBw = 0.0475;
lambda =  2.5;
ctrl_bound = 5.91;
n_W = 1;

%% Dynamics

mass = 1;
l = 1;
I = (1/3)*mass*(2*l)^2;
sigma = 100;
J = 1;
b = 1;
g = 9.81;

%Dynamics for synthesis
f  = @(x) [x(2);
          (mass*g*l/I)*sin(x(1)) - (sigma/I)*(x(1)-x(3));
           x(4);
           (sigma/J)*(x(1)-x(3)) - (b/J)*x(4)];
       
B = [zeros(3,1);1/J];

df = @(x) [0, 1, 0, 0;
           (mass*g*l/I)*cos(x(1))-(sigma/I), 0, (sigma/I), 0;
           zeros(1,3), 1;
           (sigma/J), 0, -(sigma/J), -(b/J)];

B_w = [0, (1/I), 0, 0;
       0, 0, 0, (1/J)]';

%True Dynamics   
mass_true = 0.5 + 1*rand(1);
I_true = (1/3)*mass_true*(2*l)^2;
sigma_true = 80 + (120-80)*rand(1);
J_true = J;

f_true = @(x) [x(2);
          (mass_true*g*l/I_true)*sin(x(1)) - (sigma_true/I_true)*(x(1)-x(3));
           x(4);
           (sigma_true/J_true)*(x(1)-x(3)) - (b/J)*x(4)];
B_true =  [zeros(3,1);1/J_true];
B_w_true = [0,(1/I_true), 0, 0;
            0, 0, 0, (1/J_true)]';
% f_true = f;
% B_true = B;
% B_w_true = B_w;

   
%% Bounds

w_max = 1;

M_ccm = W_upper\eye(n);
d_bar = (w_max*sigma_ThBw/lambda);
ctrl_bound = ctrl_bound*w_max;
euc_bound = d_bar*sqrt(diag(W_upper));

In = eye(n);
M_ccm_ang = (1/d_bar^2)*((In([1;3],:)*W_upper*In([1;3],:)')\eye(2));
M_ccm_omg = (1/d_bar^2)*((In([2;4],:)*W_upper*In([2;4],:)')\eye(2));

P = 2.5*eye(n);
alpha = 1e-3;

%% Simulation constraints

state_constr_low = -[pi;5;pi;5]+euc_bound;
ctrl_constr = [-35+ctrl_bound, 35-ctrl_bound];
           
q_eq = 0*(pi/180);
th_eq = q_eq - (mass*g*l/sigma)*sin(q_eq);
x_eq = [q_eq; 0; th_eq; 0];
u_eq = -sigma*(q_eq-th_eq); 

link_ang = -pi+15*(pi/180);
mot_ang = link_ang + 5*(pi/180);
test_state = [link_ang;
              -30*(pi/180);
              mot_ang;
              0];
          
%% Design FL controller

phi = @(x) [x(1); 
            x(2);
            (mass*g*l/I)*sin(x(1)) - (sigma/I)*(x(1)-x(3));
           ((mass*g*l/I)*cos(x(1)) - (sigma/I))*x(2) + (sigma/I)*x(4)];
       
f_tilde_xi = @(xi) -(mass*g*l/I)*sin(xi(1))*(xi(2)^2) + (mass*g*l/I)*cos(xi(1))*xi(3) - (sigma/I)*xi(3);
f_tilde_x = @(x) (sigma/I)*((sigma/J)*(x(1)-x(3)) - (b/J)*x(4));

f_tilde = @(x) f_tilde_x(x) + f_tilde_xi(phi(x));
B_tilde = (sigma/I)*(1/J);

% Compute beforehand
Ac = [zeros(3,1), eye(3); zeros(1,4)];
Bc = [zeros(3,1);1];
% [K_FL, M_FL, E_FL] = lqr(Ac,Bc,diag(10*[2000;1000;2000;1000]),1);
K_FL = place(Ac,Bc,-[4;5;6;7]);





