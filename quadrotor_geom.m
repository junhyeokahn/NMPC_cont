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
dt = 0.001;
t = (0:dt:T)';

[state_nom, ctrl_nom, ang_a, accel] = generate_quad_traj(t,Jq,mq,g);

ctrl_nom = [state_nom(:,7), ctrl_nom(:,2:4)];
state_nom = [state_nom(:,1:6), state_nom(:,9:14)];

%Plotting
figure()
subplot(3,1,1)
plot(t,state_nom(:,9:11)*(180/pi),'linewidth',2);
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
w_dist = [0;0.01*ones(3,1)];

%% Continuous Simulation

ode_options = odeset('RelTol', 1e-4, 'AbsTol', 1e-7);
start_p = state_nom(1,:);

[t_vec,x_act] = ode113(@(t_vec,x_act)quad_sim_geom(t_vec,x_act,...
    t,state_nom,ctrl_nom,ang_a,f,B,B_w,w_dist),...
        t,start_p,ode_options);

%% Plot

close all

%Trajectory plot
figure()
plot3(x_act(:,1),x_act(:,2),x_act(:,3),'b-','linewidth',2); hold on
plot3(state_nom(:,1),state_nom(:,2),state_nom(:,3),'r-','linewidth',2);
grid on
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

%% Animate

plot_quad_movie(state_nom(:,1),state_nom(:,2),state_nom(:,3),t,x_act,20,n)

