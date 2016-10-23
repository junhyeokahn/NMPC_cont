%% Plot
close all

% State Trajectory
figure()
hold on
plot(solve_t, x_act(:,3:6),'linewidth',2);
grid on
xlabel('Time [s]'); ylabel('States'); 
h_leg = legend('$\phi$','$v_x$','$v_z$','$\dot{\phi}$');
set(h_leg,'interpreter','latex');
set(findall(gcf,'type','text'),'FontSize',32);set(gca,'FontSize',32)

%Control Trajectory
figure()
subplot(2,1,1)
hold on
plot(0:dt:t_end,Nom_ctrl(:,1),'b-','linewidth',2);
plot(0:dt:t_end,True_ctrl(:,1),'r-','linewidth',2);
xlabel('Time [s]');
ylabel('u(t)'); 
set(findall(gcf,'type','text'),'FontSize',32);set(gca,'FontSize',32)
grid on

subplot(2,1,2)
hold on
plot(0:dt:t_end,Nom_ctrl(:,2),'b-','linewidth',2);
plot(0:dt:t_end,True_ctrl(:,2),'r-','linewidth',2);
xlabel('Time [s]');
ylabel('u(t)'); 
set(findall(gcf,'type','text'),'FontSize',32);set(gca,'FontSize',32)
grid on

figure()
plot(solve_t(1:end-1),Aux_ctrl,'linewidth',2); 
xlabel('Time [s]');
ylabel('k(x^{*},x)');
set(findall(gcf,'type','text'),'FontSize',32);set(gca,'FontSize',32)
grid on

%2D State Plot
figure()
hold on
plot(MP_state(:,1),MP_state(:,2),'g--','linewidth',1);
for i_mpc = 1:T_steps_MPC
    plot(MPC_state{i_mpc}(1:round(T_mpc/dt)+1,1),MPC_state{i_mpc}(1:round(T_mpc/dt)+1,2),'r--','linewidth',1.5);
    Ellipse_plot(M_ccm_pos,MPC_state{i_mpc}(1,1:2)',30,'k');
    Ellipse_plot(M_ccm_pos,MPC_state{i_mpc}(round(delta/dt)+1,1:2)',30,'k');
%     for j = 1:round((T_mpc/2)/dt):round(T_mpc/dt)+1
%         Ellipse_plot(M_ccm_pos,MPC_state{i_mpc}(j,1:2)',30,'k');
%     end
%     pause;
end
plot(x_act(:,1),x_act(:,2),'k-','linewidth',2);
plot(x_act(1:round(delta/dt_sim):end,1),x_act(1:round(delta/dt_sim):end,2),'ko','markersize',7,'markerfacecolor','k');

Ellipse_plot(P(1:2,1:2)*(1/(alpha)), x_eq(1:2),30,'k');
xlabel('$X$','interpreter','latex'); 
ylabel('$Z$','interpreter','latex');
set(findall(gcf,'type','text'),'FontSize',32);set(gca,'FontSize',32)
grid on; 
% axis equal

%Solve Time
figure()
hold on
plot(solve_t(1:end-1),ctrl_solve_time(:,1),'ro','markersize',10,'markerfacecolor','g');
plot(solve_t(1:end-1),ctrl_solve_time(:,2),'bs','markersize',10,'markerfacecolor','m');
plot(solve_t(1:end-1),ctrl_solve_time(:,3),'rd','markersize',10,'markerfacecolor','k');
grid on
legend('MPC','Geo','Aux');
xlabel('Time [s]');
ylabel('Solve time [s]'); title('Solve time');
set(findall(gcf,'type','text'),'FontSize',32);set(gca,'FontSize',32)

%Solve success
figure()
hold on
plot(solve_t(1:end-1),opt_solved(:,1),'ro','markersize',15,'markerfacecolor','g');
plot(solve_t(1:end-1),opt_solved(:,2),'bs','markersize',15,'markerfacecolor','m');
plot(solve_t(1:end-1),opt_solved(:,3),'rd','markersize',10,'markerfacecolor','k');
grid on
legend('MPC (1,2,3)','Geodesic (0,1,6)','Aux (0)');
xlabel('Time [s]');
title('Convergence');
set(findall(gcf,'type','text'),'FontSize',32);set(gca,'FontSize',32)
 
figure()
plot(solve_t(1:end-1),sqrt(geo_energy(:,1)),'b-','linewidth',2);
hold on
plot(solve_t(1:end-1),d_bar*ones(T_steps,1),'r-','linewidth',2);
plot(solve_t(1:end-1),sqrt(geo_energy(:,2)),'bo','markersize',10,'markerfacecolor','b','linewidth',2);
grid on
legend('d(x^{*},x)','RCI bound');
xlabel('Time [s]');
set(findall(gcf,'type','text'),'FontSize',32);set(gca,'FontSize',32)

save('quad_sim.mat');