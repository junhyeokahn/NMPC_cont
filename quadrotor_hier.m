clear all; close all; clc;
warning('off','MATLAB:nargchk:deprecated');
%% Constants

n = 12;
m = 4;

Jx = 0.082; Jy =  0.0845; Jz =  0.1377;
Jq = diag([Jx;Jy;Jz]);
a = (Jy-Jz)/Jx; b = (Jz-Jx)/Jy; c = (Jx-Jy)/Jz;
mq = 4.34;
g = 9.81;

%% Generate desired trajectory

T = 20;
dt = 0.005;
t = (0:dt:T)';

[state_nom, ctrl_nom, ang_a, accel] = generate_quad_traj(t,Jq,mq,g);
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

lambda = 1.0;
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

wq_lower >= 0.01;

Wq_12 + Wq_12' <= -2*lambda*Wq_11;
cvx_end

M_q = Wq\eye(6);

pause;

kx = 16*mq;
kv = 5.6*mq;

R_om = @(x) [1, 0, -sin(x(2));
             0, cos(x(1)), sin(x(1))*cos(x(2));
             0, -sin(x(1)), cos(x(1))*cos(x(2))]; %euler_rates -> body_rates

%Jacobian of transform from q->xi_q
phi_d = @(x) [eye(3), zeros(3);
            [0, -x(6)*cos(x(2)), 0;
           -x(5)*sin(x(1))+x(6)*cos(x(1))*cos(x(2)), -x(6)*sin(x(1))*sin(x(2)), 0;
           -x(5)*cos(x(1))-x(6)*sin(x(1))*cos(x(2)), -x(6)*cos(x(1))*sin(x(2)), 0], R_om(x)];
       
M = @(x) M_q*(phi_d(x)\eye(6));

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

%% Continuous Simulation

% ode_options = odeset('RelTol', 1e-4, 'AbsTol', 1e-7);
start_p = state_nom(1,:) + [0.05*randn(6,1);
                            3*(pi/180)*randn(6,1)]';

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

u_prev = zeros(3,1);

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
        ang_accel_des = (om_des - om_des_prev)/dt;
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
        torque_nom = Jq*ang_accel_des + cross(om_des,Jq*om_des);
        
        att_des_prev = att_des;
        om_des_prev = om_des;
    end
    
    rate_eul_act = R_om(x(7:9))\x(10:12);
    
    q_des = [att_des;rate_eul_des];
    q_act = [x(7:9);rate_eul_act];
    
    xi_q_des = [att_des;om_des];
    xi_q_act = [x(7:9);x(10:12)];
    
    Xq_dot = kron(ones(1,2),q_act-q_des);
    X_xi = [xi_q_des, xi_q_act];
    X_q = [q_des, q_act];
    
    E(i) = (q_act - q_des)'*M_q*(q_act - q_des);
    
    [aux, solved(i)] = compute_quad_aux(aux_prob,...
     X_q,Xq_dot,X_xi,E(i),M,f_rot,B_rot,torque_nom,u_prev,eps_u,lambda);
    aux_torque(i,:) = aux';
    
    ctrl(i,:) = [thrust, torque_nom'+aux'];
    u_prev = ctrl(i,2:4)';
    
    w_dist = u_nom.*[0;0.5*ones(3,1)];
    
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
plot(t,state_nom(:,1:3)-x_act(:,1:3),'linewidth',2);
grid on
xlabel('Time [s]');
ylabel('Traj errors');
legend('e_x','e_y','e_z');
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
% plot(t(1:end-1),(d_bar^2)*ones(T_steps,1),'r-','linewidth',2);
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

plot_quad_movie(state_nom(:,1),state_nom(:,2),state_nom(:,3),t,x_act,5,n)

