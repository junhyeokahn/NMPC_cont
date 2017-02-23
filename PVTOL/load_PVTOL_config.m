%% Constants

n = 6;
m = 2;

%% Obstacle info

obs_loc = [[3;-4],...
           [0.7;-3],...
           [-1;-0.5],...
           [2.5;-0.5],...
           [-4;1],...
           [1.0;1.7],...
           [2.5;3.8],...
           [-2.5;-4.5],...
           [-2;4]];
obs_loc_mpc = [[0.7;-3],...
           [-1;-0.5],...
           [2.5;-0.5],...
           [1.0;1.7],...
           [2.5;3.8],...
           [-2.5;-4.5]];
obs_rad = [1,0.9,0.8,1.2,1,0.9,0.5,1,0.6];
obs_rad_mpc = [0.9,0.8,1.2,0.9,0.5,1];
obs = struct('n_obs',length(obs_rad),'pos',obs_loc,'r',obs_rad);
obs_mpc = struct('n_obs',length(obs_rad_mpc),'pos',obs_loc_mpc,'r',obs_rad_mpc);

%% Setup Metric

load 'metric_PVTOL_vectorized.mat';

% W_fnc = @(x) W_mat(x);
% dW_fnc = @(x) {dW_p_mat(x),dW_vy_mat(x)};

W_fnc = struct('W_eval',W_eval,'w_poly_fnc',w_poly_fnc);
dW_fnc = @(x) {dw_poly_p_fnc(x), dw_poly_vy_fnc(x)};

% sigma_ThBw = 0.3296;
sigma_ThBw = 0.3185;
lambda =  0.8283;
ctrl_bound = 6.00;
n_W = [3,4];

%% Dynamics

mass = 0.486;
J = 0.00383;
g = 9.81;
len = 0.25;

f  = @(x) [x(4)*cos(x(3)) - x(5)*sin(x(3));
           x(4)*sin(x(3)) + x(5)*cos(x(3));
           x(6);
           x(6)*x(5)-g*sin(x(3));
           -x(6)*x(4)-g*cos(x(3));
           0];
       
B = [zeros(1,4),1/mass, len/J;
     zeros(1,4),1/mass,-len/J]';

df = @(x) [0,0,-x(4)*sin(x(3))-x(5)*cos(x(3)),cos(x(3)),-sin(x(3)),0;
           0,0, x(4)*cos(x(3))-x(5)*sin(x(3)),sin(x(3)), cos(x(3)),0;
           zeros(1,5),1;
           0,0,-g*cos(x(3)), 0, x(6), x(5);
           0,0, g*sin(x(3)), -x(6), 0, -x(4);
           zeros(1,6)];

B_w = [zeros(1,3),1,0,0;
       zeros(1,3),0,1,0]';
   
f_true = f;
B_true = B;
B_w_true = B_w;
   
%% Bounds

w_max = 0.1;

M_ccm = W_upper\eye(n);
d_bar = (w_max*sigma_ThBw/lambda);
ctrl_bound = ctrl_bound*w_max;
euc_bound = d_bar*sqrt(diag(W_upper));

In = eye(n);
M_ccm_pos_unscaled = ((In(1:2,:)*W_upper*In(1:2,:)')\eye(2));
M_ccm_pos = (1/d_bar^2)*((In(1:2,:)*W_upper*In(1:2,:)')\eye(2));
[U_pos,S_pos,V_pos] = svd(M_ccm_pos);
    
%Rescale ellipsoids by obstacle + robot radius
M_obs = zeros(2,2,obs.n_obs);
for i = 1:obs.n_obs
    S_new = (sqrt(S_pos\eye(2)) + (obs_rad(i)+len)*eye(2))^2\eye(2);
    M_obs(:,:,i) = U_pos*S_new*V_pos';
end
obs.M_obs = M_obs;

M_obs_mpc = zeros(2,2,obs_mpc.n_obs);
for i = 1:obs_mpc.n_obs
    S_new = (sqrt(S_pos\eye(2)) + (obs_mpc.r(i)+len)*eye(2))^2\eye(2);
    M_obs_mpc(:,:,i) = U_pos*S_new*V_pos';
end
obs_mpc.M_obs = M_obs_mpc;

P = 2.5*eye(n);
alpha = 1e-3;

%% Simulation constraints

state_constr_low = -[5.5;5.5;pi/4;2;1;pi/3]+euc_bound;
state_constr = [state_constr_low, -state_constr_low];
ctrl_constr = [0.1*mass*g+ctrl_bound, 2*mass*g-ctrl_bound;
               0.1*mass*g+ctrl_bound, 2*mass*g-ctrl_bound];
           
x_eq = [4.5;4.5;0;0;0;0];
u_eq = [0.5*mass*g; 0.5*mass*g]; 

test_state = [-4.4;
              -5;
               0;
               1.3;
               0;
               0];
