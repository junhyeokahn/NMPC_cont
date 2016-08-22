clear all; close all; clc;
warning('off','MATLAB:nargchk:deprecated');
%% Constants

n = 12;
m = 4;

Jx = 0.003; Jy =  0.004; Jz =  0.04;
Jq = diag([Jx;Jy;Jz]);
a = (Jy-Jz)/Jx; b = (Jz-Jx)/Jy; c = (Jx-Jy)/Jz;
mq = 0.93;
% Jx = 0.082; Jy =  0.0845; Jz =  0.1377;
% Jq = diag([Jx;Jy;Jz]);
% a = (Jy-Jz)/Jx; b = (Jz-Jx)/Jy; c = (Jx-Jy)/Jz;
% mq = 4.34;
g = 9.81;

%% Generate desired trajectory

T = 20;
dt = 0.002;
t = (0:dt:T)';

[state_nom, ctrl_nom, ang_a, accel] = generate_quad_traj(t,Jq,mq,g);
% state_nom = [zeros(length(t),6), mq*g*ones(length(t),1),zeros(length(t),7)];
% ctrl_nom = zeros(length(t),4);
% ang_a = zeros(length(t),3);
% accel = zeros(length(t),3);

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

% kx = 16*mq;
% kv = 5.6*mq
kx = 1*mq;
kv = 1*mq;
kR = 4;
k_om = 0.8;

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
     
dist_ratio = 0.0;

%% Continuous Simulation

% ode_options = odeset('RelTol', 1e-4, 'AbsTol', 1e-7);
start_p = state_nom(1,:) + [0.05*randn(6,1);
                            3*(pi/180)*randn(6,1)]';
% start_p = [state_nom(1,1:3)+0.1*randn(1,3), 0,0,0,0,0,40*(pi/180),0,0,0];
                        
% [t_vec,x_act] = ode113(@(t_vec,x_act)quad_sim_geom(t_vec,x_act,...
%     t,state_nom,ctrl_nom,ang_a,f,B,B_w,dist_ratio),...
%         t,start_p,ode_options);

%% Discrete Simulation

T_steps = length(t)-1;

t_opt = cell(T_steps,1);
state = cell(T_steps,1);

x_act = zeros(T_steps+1,n);
x_act(1,:) = start_p';

ctrl = zeros(T_steps,m);


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
    xb_des = [cos(pitch_des)*cos(yaw_des);
              cos(pitch_des)*sin(yaw_des);
             -sin(pitch_des)];
    yb_des = cross(zb_des,xb_des);

    R_des = [xb_des,yb_des,zb_des]; 
%     om_des = x_nom(10:12);
%     ang_accel_des = ang_a(i,1:3)';
%     
    if i>1
        R_des_rate = (R_des-R_des_prev)/dt;
        om_des = vee(R_des_prev'*R_des_rate);
        ang_accel_des = (om_des - om_des_prev)/dt;
        
        R_des_prev = R_des;
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
        xb_des_pred = [cos(pitch_des_pred)*cos(yaw_des_pred);
                       cos(pitch_des_pred)*sin(yaw_des_pred);
                      -sin(pitch_des_pred)];
        yb_des_pred = cross(zb_des_pred,xb_des_pred);

        R_des_pred = [xb_des_pred,yb_des_pred,zb_des_pred];
        
        R_des_rate = (R_des_pred-R_des)/dt;
        om_des = vee(R_des'*R_des_rate);
        
        R_nom = rot_matrix(x_nom(9),x_nom(8),x_nom(7));
        ang_accel_des = R_des'*R_nom*ang_a(1,:)';
        
        R_des_prev = R_des;
        om_des_prev = om_des;
    end
    
    eR = vee(0.5*(R_des'*R - R'*R_des));
    om = x(10:12);
    
    e_om = om - R'*R_des*om_des;
    
    torque = -kR*eR - k_om*e_om + Skew(om)*Jq*om -...
         Jq*(Skew(om)*R'*R_des*om_des - R'*R_des*ang_accel_des);
 
    
    ctrl(i,:) = [thrust, torque'];
    
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
ylabel('Traj errors');
legend('e_x','e_y','e_z');
set(findall(gcf,'type','text'),'FontSize',32);set(gca,'FontSize',32)

subplot(2,1,2)
plot(t,(state_nom(:,9)-x_act(:,9))*(180/pi),'linewidth',2);
grid on
xlabel('Time [s]');
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


%% Animate

plot_quad_movie(state_nom(:,1),state_nom(:,2),state_nom(:,3),t,x_act,round(0.02/dt),n)

