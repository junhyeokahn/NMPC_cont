%% Constants
n = 2;
m = 1;

%% Obstacle info

obs = struct('n_obs',0);
obs_mpc = struct('n_obs',0);

%% Setup Metric

load 'metric_Allgower.mat';

% W_fnc = @(x) W_mat(wrapToPi(x(1)));
% dW_fnc = @(x) {dW_x1_mat(wrapToPi(x(1)))};

w_poly = @(x) w_poly_fnc(x);
W_fnc = struct('W_eval',W_eval,'w_poly_fnc',w_poly);
dW_fnc = @(x) {};

sigma_ThBw = 15.6182;
lambda =  1.7429;
ctrl_bound = 1.908;
n_W = []; %states that W is a function of

%% Dynamics & Cost

f = @(x) [-1*x(1) + 2*x(2); 
          -3*x(1) + 4*x(2) - 0.25*(x(2)^3)];
df = @(x) [-1, 2;
           -3, 4-0.75*(x(2)^2)];
B = [0.5;-2];

B_w = [0;1];

Q = diag([0.5;0.5]); R = 1;

%% Bounds

w_max = 0.1;

M_ccm = W_fnc.W_eval(1)\eye(n);
d_bar = (w_max*sigma_ThBw/lambda);
ctrl_bound = ctrl_bound*w_max;
euc_bound = d_bar*sqrt(diag(W_fnc.W_eval(1)));

alpha = 10;
P = [7.9997, -12.2019;
    -12.2019, 27.0777];

%% Simulation constraints

state_constr_low = -[5;5]+euc_bound;
state_constr = [state_constr_low, -state_constr_low];
ctrl_constr = [-2+ctrl_bound,2-ctrl_bound];

x_eq = [0;0];
u_eq = 0;
test_state = [3.4;-2.4];

