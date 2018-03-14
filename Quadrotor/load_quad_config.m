%% Constants

Jx = 0.003; Jy =  0.004; Jz =  0.04;
global mq Jq;

Jq = diag([Jx;Jy;Jz]);
mq = 0.9574;

g = 9.81;

n = 9;
m = 3;

%% Dynamics

%pos,vel, roll, pitch, yaw, om, thrust_sp,ang_vel_err_state
f = @(x) [x(4:6);
          [0;0;g;]-[sin(x(8)); -cos(x(8))*sin(x(7)); cos(x(8))*cos(x(7))]*x(13);
          R_eul(x(7:9))*x(10:12);
          -Jq\cross(x(10:12),Jq*x(10:12));
          0];
      
B =      [zeros(9,4);
          zeros(3,1), Jq\eye(3);
          1,zeros(1,3)];
  
f_ctrl = @ f_xc;      
B_ctrl = [zeros(6,3);
          eye(3)];

B_w = [zeros(3);
       eye(3);
       zeros(7,3)];
   
w_max = 0.1;   
   
%% Generate desired trajectory

%Generate trajectory
poly_file = 'soln_traj.mat';

load(poly_file);

%% Setup metric 

pullback = 0;

if (~pullback)
    %SOS METRIC
    
    %Metric Controller state-space: xc: x,y,z, vx,vy,vz, th, r,p
    %Metric Controller control-space: uc: rd, pd, th_dot
    
    load ('metric_QUAD_vectorized.mat');
    
    W_fnc = struct('W_eval',W_eval,'w_poly_fnc',w_poly_fnc);
    dW_fnc = @(x) {dw_poly_th_fnc(x), dw_poly_r_fnc(x), dw_poly_p_fnc(x)};
    n_W = [7,8,9];
    
    sigma_ThBw = 0.3522;
    lambda = 0.84;
    
    %Bounds
    
    W_upper = W_upper_mat;
    M_ccm = W_upper\eye(n);
    d_bar = (w_max*sigma_ThBw/lambda);
    
    In = eye(n);
    %Maximal position tube
    M_ccm_pos = (1/d_bar^2)*((In(1:3,:)*W_upper*In(1:3,:)')\eye(3));
    [U_pos,S_pos,V_pos] = svd(M_ccm_pos);
    
else
    % PULLBACK METRIC
    % Auxiliary Metric state-space: xic = [x,y,z,vx,vy,vz,ax,ay,az] := phi(xc)
    load 'metric_QUAD_pullback.mat';
    In = eye(n);
    P = In(1:3,:);
    Ap = (P*(M_xi\P'))\eye(3);
    M_ccm_pos = (1/d_bar^2)*Ap;    
end

%% Setup lower-level controller

%Angular rate PI controller gains
global kp_om ki_om;
kp_om = 3*lambda;
ki_om = 0.1;



