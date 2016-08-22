clear all; close all; clc;
warning('off','MATLAB:nargchk:deprecated');
%% Constants

n = 12;
m = 4;

Jx = 0.003; Jy =  0.004; Jz =  0.04;
Jq = diag([Jx;Jy;Jz]);
a = (Jy-Jz)/Jx; b = (Jz-Jx)/Jy; c = (Jx-Jy)/Jz;
mq = 0.9574;
% Jx = 0.082; Jy =  0.0845; Jz =  0.1377;
% Jq = diag([Jx;Jy;Jz]);
% a = (Jy-Jz)/Jx; b = (Jz-Jx)/Jy; c = (Jx-Jy)/Jz;
% mq = 4.34;
g = 9.81;

%% Generate desired trajectory

T = 20;
dt = 0.008;
t = (0:dt:T)';

% [state_nom, ctrl_nom, ang_a, accel] = generate_quad_traj(t,Jq,mq,g);
state_nom = [zeros(length(t),6), mq*g*ones(length(t),1),zeros(length(t),7)];
ctrl_nom = zeros(length(t),4);
ang_a = zeros(length(t),3);
accel = zeros(length(t),3);

% Reduce to 12 dimensional
ctrl_nom = [state_nom(:,7), ctrl_nom(:,2:4)];
state_nom = [state_nom(:,1:6), state_nom(:,9:14)];

%Plotting
figure()
subplot(3,1,1)
plot(t,state_nom(:,7:9)*(180/pi),'linewidth',2);
grid on
xlabel('Time [s]'); ylabel('Attitude [deg]');
legend('Roll','Pitch','Yaw');

subplot(3,1,2)
plot(t,state_nom(:,10:12)*(180/pi),'linewidth',2);
grid on
xlabel('Time [s]'); ylabel('Angular rate [deg/s]');
legend('p','q','r');

subplot(3,1,3)
plot(t,ang_a(:,1:3)*(180/pi),'linewidth',2);
grid on
xlabel('Time [s]'); ylabel('Angular rate [deg/s^2]');
legend('pd','qd','rd');

%% Controller setup

eps_u = 0;
aux_prob = setup_opt_aux(3,eps_u);
geo_Ke = 1;

lambda = 2.0;
cvx_begin sdp
variable Wq_11(3,3) symmetric
variable Wq_12(3,3)
variable Wq_22(3,3) symmetric
variables wq_upper wq_lower
minimize (wq_upper - wq_lower)
subject to
Wq = [Wq_11, Wq_12;
     Wq_12', Wq_22];
Wq >= wq_lower*eye(6);
Wq <= wq_upper*eye(6);

wq_lower >= 1;

Wq_12 + Wq_12' <= -2*lambda*Wq_11;
cvx_end
% M_q = [0.87868*eye(3), 0.292893*eye(3);
%     0.292893*eye(3), 0.292893*eye(3)];

M_q = clean(Wq\eye(6),1e-3);
Mq_11 = M_q(1:3,1:3);
Mq_12 = M_q(1:3,4:6);
Mq_22 = M_q(4:6,4:6);

pause;

% kx = 16*mq;
% kv = 5.6*mq;
kx = 1*mq;
kv = 1*mq;

R_om = @(x) [1, 0, -sin(x(2));
    0, cos(x(1)), sin(x(1))*cos(x(2));
    0, -sin(x(1)), cos(x(1))*cos(x(2))]; %euler_rates -> body_rates

R_om_inv = @(x) [1, sin(x(1))*tan(x(2)), cos(x(1))*tan(x(2));
    0, cos(x(1)), -sin(x(1));
    0, sin(x(1))*sec(x(2)), cos(x(1))*sec(x(2))];

%Jacobian of transform from xi_q -> q
Psi_jac = @(xi_q) [eye(3), zeros(3);
    [(xi_q(5)*cos(xi_q(1))-xi_q(6)*sin(xi_q(1)))*tan(xi_q(2)), (xi_q(5)*sin(xi_q(1))+xi_q(6)*cos(xi_q(1)))*(sec(xi_q(2)))^2, 0;
    -xi_q(5)*sin(xi_q(1))-xi_q(6)*cos(xi_q(1)), 0, 0;
    (xi_q(5)*cos(xi_q(1))-xi_q(6)*sin(xi_q(1)))*sec(xi_q(2)), (xi_q(5)*sin(xi_q(1))+xi_q(6)*cos(xi_q(1)))*sec(xi_q(2))*tan(xi_q(2)),0], R_om_inv(xi_q)];

M = @(xi) M_q*(Psi_jac(xi));

%% Dynamics

%x:[x,y,z,vx,vy,vz,phi,th,psi,p,q,r]

roll_d = @(x) x(10) + x(11)*tan(x(8))*sin(x(7)) + x(12)*tan(x(8))*cos(x(7));
pitch_d = @(x) x(11)*cos(x(7)) - x(12)*sin(x(7));
yaw_d = @(x) x(11)*sec(x(8))*sin(x(7)) + x(12)*sec(x(8))*cos(x(7));

bq_1 = @(x) sin(x(7))*sin(x(9)) + cos(x(7))*sin(x(8))*cos(x(9));
bq_2 = @(x) -sin(x(7))*cos(x(9)) + cos(x(7))*sin(x(8))*sin(x(9));
bq_3 = @(x) cos(x(7))*cos(x(8));

f = @(x) [x(4);
    x(5);
    x(6);
    0;
    0;
    g;
    roll_d(x);
    pitch_d(x);
    yaw_d(x);
    a*x(11)*x(12);
    b*x(10)*x(12);
    c*x(10)*x(11)];

B = @(x)[zeros(3,4);
    (-1/mq)*bq_1(x), zeros(1,3);
    (-1/mq)*bq_2(x), zeros(1,3);
    (-1/mq)*bq_3(x), zeros(1,3);
    zeros(3,4);
    zeros(3,1), (Jq\eye(3))];

B_w = B;

f_rot = @(x)[x(4) + x(5)*tan(x(2))*sin(x(1)) + x(6)*tan(x(2))*cos(x(1));
    x(5)*cos(x(1)) - x(6)*sin(x(1));
    x(5)*sec(x(2))*sin(x(1)) + x(6)*sec(x(2))*cos(x(1));
    a*x(5)*x(6);
    b*x(4)*x(6);
    c*x(4)*x(5)];
B_rot = [zeros(3,3);
    (Jq\eye(3))];

dist_ratio = 0;
w_max = max(norms(dist_ratio*ctrl_nom(:,2:4),2,2));
sigma_Bw = sqrt(max(eig(B_rot'*B_rot)));
d_bar = w_max*sigma_Bw*sqrt(max(eig(M_q)))/lambda;

%% Continuous Simulation

% ode_options = odeset('RelTol', 1e-4, 'AbsTol', 1e-7);
start_p = state_nom(1,:) + [0.1*randn(6,1);
                            5*(pi/180)*randn(6,1)]';
% start_p = [state_nom(1,1:3)+0.1*randn(1,3), 0,0,0,0,0,40*(pi/180),0,0,0];

%% Discrete Simulation

T_steps = length(t)-1;

t_opt = cell(T_steps,1);
state = cell(T_steps,1);

x_act = zeros(T_steps+1,n);
x_act(1,:) = start_p';

ctrl = zeros(T_steps,m);
aux_torque = zeros(T_steps,3);

solved = ones(T_steps,1);

E = zeros(T_steps,1);
E_rate = zeros(T_steps,1);

u_prev = zeros(3,1);
att_des_prev = zeros(3,1);
om_des_prev = zeros(3,1);

for i = 1:T_steps
    
    %     fprintf('%d/%d \n',i, T_steps);
    
    x_nom = state_nom(i,:)';
    u_nom = ctrl_nom(i,:)';
    
    x = x_act(i,:)';
    
    ex = x(1:3) - x_nom(1:3);
    ev = x(4:6) - x_nom(4:6);
    
    a_des = accel(i,:)';
    thrust_des = -kx*ex - kv*ev - mq*g*[0;0;1] + mq*a_des;
    
    R = rot_matrix(x(9),x(8),x(7));
    thrust = thrust_des'*(R*[0;0;-1]);
    
    zb_des = -thrust_des/norm(thrust_des);
    yaw_des = x_nom(9);
    pitch_des = atan2(zb_des(1)*cos(yaw_des) + zb_des(2)*sin(yaw_des),...
        zb_des(3));
    roll_des = atan2(zb_des(1)*sin(yaw_des) - zb_des(2)*cos(yaw_des),...
        zb_des(3)/cos(pitch_des));
    
    att_des = [roll_des;pitch_des;yaw_des];
    
    if i>1
        rate_eul_des = (att_des - att_des_prev)/dt;
        om_des = R_om(att_des)*rate_eul_des;
        R_prev = rot_matrix(att_des_prev(3),att_des_prev(2),att_des_prev(1));
        R_des = rot_matrix(att_des(3),att_des(2),att_des(1));
        
        ang_accel_des = (om_des - R_des'*R_prev*om_des_prev)/dt;
        torque_nom = Jq*ang_accel_des + cross(om_des,Jq*om_des);
        
        att_des_prev = att_des;
        om_des_prev = om_des;
    else
        a_pred = [0;0;g] + (-1/mq)*[bq_1(x);bq_2(x);bq_3(x)]*thrust;
        xv_pred = [x(1:3); x(4:6)] + dt*[x(4:6);a_pred];
        
        ex_pred = xv_pred(1:3) - state_nom(2,1:3)';
        ev_pred = xv_pred(4:6) - state_nom(2,4:6)';
        a_des = accel(2,:)';
        thrust_des_pred = -kx*ex_pred - kv*ev_pred - mq*g*[0;0;1] + mq*a_des;
        zb_des_pred = -thrust_des_pred/norm(thrust_des_pred);
        yaw_des_pred = state_nom(2,9);
        pitch_des_pred = atan2(zb_des_pred(1)*cos(yaw_des_pred) + zb_des_pred(2)*sin(yaw_des_pred),...
            zb_des_pred(3));
        roll_des_pred = atan2(zb_des_pred(1)*sin(yaw_des_pred) - zb_des_pred(2)*cos(yaw_des_pred),...
            zb_des_pred(3)/cos(pitch_des_pred));
        
        att_des_pred = [roll_des_pred;pitch_des_pred;yaw_des_pred];
        
        rate_eul_des = (att_des_pred-att_des)/dt;
        om_des = R_om(att_des)*rate_eul_des;
        
        R_nom = rot_matrix(x_nom(9),x_nom(8),x_nom(7));
        R_des = rot_matrix(att_des(3),att_des(2),att_des(1));
        ang_accel_des = R_des'*R_nom*ang_a(1,:)';
%         ang_accel_des = zeros(3,1);
        torque_nom = Jq*ang_accel_des + cross(om_des,Jq*om_des);
        
        att_des_prev = att_des;
        om_des_prev = om_des;
    end
    
    rate_eul_act = R_om(x(7:9))\x(10:12);
    
    
    q_des = [att_des;rate_eul_des];
    q_act = [x(7:9);rate_eul_act];
    
    xi_q_des = [att_des;om_des];
    xi_q_act = [x(7:9);x(10:12)];
    
    att_act = x(7:9);
    om_act = x(10:12);
    
%     Xq_dot = kron(ones(1,2),q_act-q_des);
%     X_xi = [xi_q_des, xi_q_act];
%     X_q = [q_des, q_act];
%     
%     E(i) = (q_act - q_des)'*M_q*(q_act - q_des);
%     
%     [aux, solved(i), E_rate(i)] = compute_quad_aux(aux_prob,...
%         X_q,Xq_dot,X_xi,E(i),M,f_rot,B_rot,torque_nom,u_prev,eps_u,lambda);
%     aux_torque(i,:) = aux';
    
    q_err = att_act - att_des;
    qdot_err = rate_eul_act - rate_eul_des;
    
    Mq_Xq_dot_1 = Mq_11*q_err + Mq_12*qdot_err;
    Mq_Xq_dot_2 = Mq_12*q_err + Mq_22*qdot_err;
    E(i) = q_err'*Mq_Xq_dot_1 + qdot_err'*Mq_Xq_dot_2;
    
    % Setup control matrices
    R_om_inv_des = R_om_inv(att_des);
    R_om_inv_act = R_om_inv(x(7:9));
    
    % Normalize by inertia
    R_om_inv_act(1,1) = R_om_inv_act(1,1)/Jx;
    R_om_inv_act(1,2) = R_om_inv_act(1,2)/Jy;
    R_om_inv_act(2,2) = R_om_inv_act(2,2)/Jy;
    R_om_inv_act(3,2) = R_om_inv_act(3,2)/Jy;
    R_om_inv_act(1,3) = R_om_inv_act(1,3)/Jz;
    R_om_inv_act(2,3) = R_om_inv_act(2,3)/Jz;
    R_om_inv_act(3,3) = R_om_inv_act(3,3)/Jz;
        
    b_vec = R_om_inv_act'*Mq_Xq_dot_2;
    
    % Setup dynamics matrices
    Psi_jac_q_des = [(om_des(2)*cos(att_des(1))-om_des(3)*sin(att_des(1)))*tan(att_des(2)), (om_des(2)*sin(att_des(1))+om_des(3)*cos(att_des(1)))/(cos(att_des(2))^2.0), 0.0;
        -om_des(2)*sin(att_des(1))-om_des(3)*cos(att_des(1)),0.0,0.0;
        ( om_des(2)*cos(att_des(1))-om_des(3)*sin(att_des(1)))/cos(att_des(2)), (om_des(2)*sin(att_des(1))+om_des(3)*cos(att_des(1)))*tan(att_des(2))/cos(att_des(2)), 0.0];
    
    Psi_jac_q_act = [(om_act(2)*cos(att_act(1))-om_act(3)*sin(att_act(1)))*tan(att_act(2)), (om_act(2)*sin(att_act(1))+om_act(3)*cos(att_act(1)))/(cos(att_act(2))^2.0), 0.0;
        -om_act(2)*sin(att_act(1))-om_act(3)*cos(att_act(1)),0.0,0.0;
        ( om_act(2)*cos(att_act(1))-om_act(3)*sin(att_act(1)))/cos(att_act(2)), (om_act(2)*sin(att_act(1))+om_act(3)*cos(att_act(1)))*tan(att_act(2))/cos(att_act(2)), 0.0];
    
    f_des_2 = Psi_jac_q_des*rate_eul_des + R_om_inv_des*ang_accel_des;
    
    f_act_2 = [(Jy-Jz)*om_act(2)*om_act(3),(Jz-Jx)*om_act(1)*om_act(3),(Jx-Jy)*om_act(1)*om_act(2)]';
    
    f_act_2 = Psi_jac_q_act*rate_eul_act + R_om_inv_act*f_act_2;
    
    u_b = -E(i)*1.0 + (Mq_Xq_dot_1'*(-qdot_err)) + ...
                          Mq_Xq_dot_2'*(f_des_2 - f_act_2) -...
                         b_vec'*torque_nom;
    
    a = -(2.0)*u_b;
    b_vec = 2*b_vec;
    
    if (abs(a)<= 1e-5)
        a = 0;
    end
    for j = 1:3
        if (abs(b_vec(j)) <= 1e-5)
            b_vec(j) = 0.0;
        end
    end
    
    b_norm = norm(b_vec);
    if (a <= 0.0) || (b_norm <= 1e-4)
        aux = zeros(3,1);
    else
        rho = -(a / (b_norm*b_norm));
        aux = rho*b_vec;
    end    
    
    ctrl(i,:) = [thrust, torque_nom'+aux'];
    u_prev = ctrl(i,2:4)';
    
    w_dist = u_nom.*[0;dist_ratio*ones(3,1)];
    
    x_act(i+1,:) = x'+(f(x) + B(x)*ctrl(i,:)' + B_w(x)*w_dist)'*dt;
    
end


%% Plot

close all

%Trajectory plot
figure()
plot3(x_act(:,1),x_act(:,2),x_act(:,3),'b-','linewidth',2); hold on
plot3(state_nom(:,1),state_nom(:,2),state_nom(:,3),'r-','linewidth',2);
grid on; axis tight;
xlabel('x'); ylabel('y'); zlabel('h');
set(gca,'ZDir','Reverse');
set(gca,'YDir','Reverse');
set(findall(gcf,'type','text'),'FontSize',32);set(gca,'FontSize',32)

%Trajectory errors
figure()
subplot(2,1,1)
plot(t,state_nom(:,1:3)-x_act(:,1:3),'linewidth',2);
grid on
xlabel('Time [s]');
ylabel('[m]');
legend('e_x','e_y','e_z');
set(findall(gcf,'type','text'),'FontSize',32);set(gca,'FontSize',32)

subplot(2,1,2)
plot(t,(state_nom(:,9)-x_act(:,9))*(180/pi),'linewidth',2);
grid on
xlabel('Time [s]'); ylabel('[deg]');
legend('e_\psi');
set(findall(gcf,'type','text'),'FontSize',32);set(gca,'FontSize',32)

%Control effort
figure()
plot(t,ctrl_nom(:,1),'r-','linewidth',2); hold on
plot(t(1:end-1),ctrl(:,1),'b-','linewidth',2);
xlabel('Time [s]');
ylabel('Thrust [N]');
grid on
legend('nominal','net');
set(findall(gcf,'type','text'),'FontSize',32);set(gca,'FontSize',32)

figure()
subplot(3,1,1)
plot(t, ctrl_nom(:,2),'--','linewidth',2); hold on
plot(t(1:end-1),ctrl(:,2),'-','linewidth',2);
xlabel('Time [s]');
ylabel('Torque [Nm]');
grid on
set(findall(gcf,'type','text'),'FontSize',32);set(gca,'FontSize',32)

subplot(3,1,2)
plot(t, ctrl_nom(:,3),'--','linewidth',2); hold on
plot(t(1:end-1),ctrl(:,3),'-','linewidth',2);
xlabel('Time [s]');
ylabel('Torque [Nm]');
grid on
set(findall(gcf,'type','text'),'FontSize',32);set(gca,'FontSize',32)

subplot(3,1,3)
plot(t, ctrl_nom(:,4),'--','linewidth',2); hold on
plot(t(1:end-1),ctrl(:,4),'-','linewidth',2);
xlabel('Time [s]');
ylabel('Torque [Nm]');
grid on
set(findall(gcf,'type','text'),'FontSize',32);set(gca,'FontSize',32)

% Geodesic Energy
figure()
plot(t(1:end-1),E,'b-','linewidth',2); hold on
plot(t(1:end-1),(d_bar^2)*ones(T_steps,1),'k-','linewidth',2);
plot(t(1:end-1),E(1)*exp(-2*lambda*t(1:end-1)),'r-','linewidth',2);
grid on
xlabel('Time [s]'); ylabel('Energy');
set(findall(gcf,'type','text'),'FontSize',32);set(gca,'FontSize',32)

% Solve success
figure()
plot(t(1:end-1),solved(:,1),'go','markersize',10,'markerfacecolor','g');
grid on
xlabel('Time [s]');
set(findall(gcf,'type','text'),'FontSize',32);set(gca,'FontSize',32)



%% Animate

plot_quad_movie(state_nom(:,1),state_nom(:,2),state_nom(:,3),t,x_act,round(0.02/dt),n)

