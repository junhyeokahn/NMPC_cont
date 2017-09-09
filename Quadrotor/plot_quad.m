
%% Trajectory plot
figure()
plot3(x_act(:,1),x_act(:,2),x_act(:,3),'b-','linewidth',2); hold on
plot3(MP_state(:,1),MP_state(:,2),MP_state(:,3),'r-','linewidth',2);
plot3(x_act(1,1),x_act(1,2),x_act(1,3),'go','markersize',15,'markerfacecolor','g');
grid on; axis tight;
xlabel('x'); ylabel('y'); zlabel('h');
set(gca,'ZDir','Reverse');
set(gca,'YDir','Reverse');
set(findall(gcf,'type','text'),'FontSize',32);set(gca,'FontSize',32)

% Obstacles
close all;
figure()
goal.plot('color','green','alpha',0.3)
hold on
for i = 1:length(obs)
    obs{i}.plot('color','black','alpha',0.3);
end

%% Trajectory errors
figure()
subplot(2,1,1)
plot(t,MP_state(1:(dt_sim/dt):end,1:3)-x_act(1:end,1:3),'linewidth',2);
grid on
xlabel('Time [s]');
ylabel('[m]');
legend('e_x','e_y','e_z');
set(findall(gcf,'type','text'),'FontSize',32);set(gca,'FontSize',32)

subplot(2,1,2)
plot(t,(MP_state(1:(dt_sim/dt):end,9)-x_act(:,9))*(180/pi),'linewidth',2);
grid on
xlabel('Time [s]'); ylabel('[deg]');
legend('e_\psi');
set(findall(gcf,'type','text'),'FontSize',32);set(gca,'FontSize',32)

%% Control effort
figure()
hold on
% plot(t,thrust_nom,'r-','linewidth',1.5);
plot(t,MP_state(:,7)*mq,'b-','linewidth',2);
plot(t,x_act(:,13)*mq,'k-','linewidth',2);
xlabel('Time [s]');
ylabel('Thrust [N]');
grid on
legend('net','commanded');
set(findall(gcf,'type','text'),'FontSize',32);set(gca,'FontSize',32)


figure()
subplot(3,1,1)
plot(t, Nom_ctrl(:,2),'--','linewidth',2); hold on
plot(t ,True_ctrl(:,2),'-','linewidth',2);
xlabel('Time [s]');
ylabel('$\dot{\phi}$ [rad/s]','interpreter','latex');
grid on
set(findall(gcf,'type','text'),'FontSize',32);set(gca,'FontSize',32)

subplot(3,1,2)
plot(t, Nom_ctrl(:,3),'--','linewidth',2); hold on
plot(t, True_ctrl(:,3),'-','linewidth',2);
xlabel('Time [s]');
ylabel('$\dot{\theta}$ [rad/s]','interpreter','latex');
grid on
set(findall(gcf,'type','text'),'FontSize',32);set(gca,'FontSize',32)

subplot(3,1,3)
plot(t, Nom_ctrl(:,4),'--','linewidth',2); hold on
plot(t, True_ctrl(:,4),'-','linewidth',2);
xlabel('Time [s]');
ylabel('$\dot{\psi}$ [rad/s]','interpreter','latex');
grid on
set(findall(gcf,'type','text'),'FontSize',32);set(gca,'FontSize',32)

%% Rotation angles
figure()
subplot(2,1,1)
plot(t, MP_state(:,8)*(180/pi),'--','linewidth',2); hold on
plot(t, x_act(:,7)*(180/pi),'-','linewidth',2); hold on
xlabel('Time [s]');
ylabel('$\phi$ [deg]','interpreter','latex');
grid on
set(findall(gcf,'type','text'),'FontSize',32);set(gca,'FontSize',32)

subplot(2,1,2)
plot(t, MP_state(:,9)*(180/pi),'--','linewidth',2); hold on
plot(t, x_act(:,8)*(180/pi),'-','linewidth',2); hold on
xlabel('Time [s]');
ylabel('$\theta$ [deg]','interpreter','latex');
grid on
set(findall(gcf,'type','text'),'FontSize',32);set(gca,'FontSize',32)

%% Geodesic Energy
figure()
plot(solve_t(1:end-1),geo_energy,'b-','linewidth',2); hold on
plot(solve_t(1:end-1),(d_bar^2)*ones(length(solve_t)-1,1),'k-','linewidth',2);
% plot(t(1:end-1),E(1)*exp(-2*lambda*t(1:end-1)),'r-','linewidth',2);
grid on
xlabel('Time [s]'); ylabel('Energy');
set(findall(gcf,'type','text'),'FontSize',32);set(gca,'FontSize',32)

%% Computation time

figure()
hold on
plot(solve_t(1:end-1),ctrl_solve_time(:,1),'bs','markersize',10,'markerfacecolor','m');
plot(solve_t(1:end-1),ctrl_solve_time(:,2),'rd','markersize',10,'markerfacecolor','k');
grid on
legend('Geo','Aux');
xlabel('Time [s]');
ylabel('Solve time [s]'); title('Solve time');
set(findall(gcf,'type','text'),'FontSize',32);set(gca,'FontSize',32)

%% Solve success

figure()
hold on
plot(solve_t(1:end-1),opt_solved(:,1),'bs','markersize',15,'markerfacecolor','m');
plot(solve_t(1:end-1),opt_solved(:,2),'rd','markersize',10,'markerfacecolor','k');
grid on
legend('Geodesic (0,1,6)','Aux (0)');
xlabel('Time [s]');
title('Convergence');
set(findall(gcf,'type','text'),'FontSize',32);set(gca,'FontSize',32)
 

%%

keyboard;
%% Animate
close all;
plot_quad_movie(MP_state(:,1),MP_state(:,2),MP_state(:,3),solve_t,x_act,round(0.1/dt_sim),M_pos,goal,obs)

%% Save 

% save('quad_sim.mat');