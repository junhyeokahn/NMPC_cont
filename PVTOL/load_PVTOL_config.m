%% Constants

n = 6;
m = 2;

%% Obstacle info

%all obstacles
obs_loc = [[-4;-3.5],...
           [3;-4],...
           [0.7;-3],...
          [-1;-0.5],...
           [2.5;-0.5],...
           [-4;1],...
           [1.0;1.7],...
           [2.5;3.8],...
           [-2;4]];

obs_rad = [0.5,1,0.9,0.8,0.9,1,0.9,0.5,0.6];       

%obstacles considered during MPC iterations (can ignore some) 
obs_loc_mpc = obs_loc(:,[1,3,4,5,7,8]);  
obs_rad_mpc = obs_rad([1;3;4;5;7;8]);

obs = struct('n_obs',length(obs_rad),'pos',obs_loc,'r',obs_rad);
obs_mpc = struct('n_obs',length(obs_rad_mpc),'pos',obs_loc_mpc,'r',obs_rad_mpc);
obs_adapt = obs;

%% Setup Metric

load 'metric_PVTOL_vectorized.mat';

W_fnc = struct('W_eval',W_eval,'w_poly_fnc',w_poly_fnc);
dW_fnc = @(x) {dw_poly_p_fnc(x), dw_poly_vy_fnc(x)};

sigma_ThBw = 0.3127;
lambda =  0.8283;
ctrl_bound = 6.00;
n_W = [3,4];

%% Dynamics & Cost

mass = 0.486;
J = 0.00383;
g = 9.81;
len = 0.25;

f  = @f_pvtol;
       
B = [zeros(1,4),1/mass, len/J;
     zeros(1,4),1/mass,-len/J]';

df = @df_pvtol;

B_w = [zeros(1,3),1,0,0;
       zeros(1,3),0,1,0]';
   
f_true = f;
B_true = B;
B_w_true = B_w;

Q = 0*diag([1;1;0;0;0;0]); R = eye(m); Q_T = 0*diag([1;1;0;0;0;0]);
   
%% Bounds

w_max = 0.1;

M_ccm = W_upper\eye(n);
d_bar = (w_max*sigma_ThBw/lambda);
ctrl_bound = ctrl_bound*w_max;
euc_bound = d_bar*sqrt(diag(W_upper));

In = eye(n);
%Unscaled position tube
M_ccm_pos_unscaled = ((In(1:2,:)*W_upper*In(1:2,:)')\eye(2)); 
%Maximal position tube
M_ccm_pos = (1/d_bar^2)*((In(1:2,:)*W_upper*In(1:2,:)')\eye(2)); 
[U_pos,S_pos,V_pos] = svd(M_ccm_pos);
    
%Rescale maximal ellipsoid by obstacle + robot radius
M_obs = zeros(2,2,obs.n_obs);
for i = 1:obs.n_obs
    S_new = (sqrt(S_pos\eye(2)) + (obs_rad(i)+len)*eye(2))^2\eye(2);
    M_obs(:,:,i) = U_pos*S_new*V_pos';
end
obs.M_obs = M_obs;

%Setup unscaled ellipsoids
[U_u,S_u,V_u] = svd(M_ccm_pos_unscaled);
obs_mpc.U = U_u; obs_mpc.V = V_u; obs_mpc.S = S_u;
obs_mpc.r = obs_rad_mpc + len;

obs_adapt.U = U_u; obs_adapt.V = V_u; obs_adapt.S = S_u;
obs_adapt.r = obs_rad + len;
adapt_EPS = struct('lambda',lambda,'d_bar',d_bar,...
                   'r',1/1.7);

%final state constraint
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
               0.5;
               0;
               0];

