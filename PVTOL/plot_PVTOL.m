%% Plot

%% 2D State Plot

close all; 

figure(1)
hold on
%Nominal motion plan
plot(MP_state(:,1),MP_state(:,2),'g--','linewidth',1);

%Plot obstacles
for i_ob = 1:obs.n_obs
    Ellipse_plot(eye(2)*(1/(obs.r(i_ob)+len)^2),obs.pos(:,i_ob), 25,'r',0.7);
end

if (~track_traj)
    for i_mpc = 1:T_steps_MPC
        %MPC reference trajectory executed segement
        plot(MPC_state{i_mpc}(1:round(delta/dt)+1,1),MPC_state{i_mpc}(1:round(delta/dt)+1,2),'r--','linewidth',1.5);
        
%         if i_mpc<T_steps_MPC
%             num_pnts = 3;
%         else
%             num_pnts = 5;
%         end
        num_pnts = delta/(4*dt_sim);
        
        %Outer geodesic ball at start of MPC segment
        E_start = geo_energy(1+(i_mpc-1)*(delta/dt_sim),2);
        Ellipse_plot(M_ccm_pos_unscaled*(1/E_start),MPC_state{i_mpc}(1,1:2)',30,'g',0.2);
%         Ellipse_plot(M_ccm_pos,MPC_state{i_mpc}(1,1:2)',30,'k');
        
        t_mpc_span = linspace(0,delta,num_pnts);
        %Evolution of outer geodesic ball over MPC segment
        for j = 2:length(t_mpc_span)-1
            E_j = (sqrt(E_start)*exp(-lambda*t_mpc_span(j)) + d_bar*(1-exp(-lambda*t_mpc_span(j))))^2;
            Ellipse_plot(M_ccm_pos_unscaled*(1/E_j),MPC_state{i_mpc}(round(t_mpc_span(j)/dt)+1,1:2)',30,'b',0.2);
%             Ellipse_plot(M_ccm_pos,MPC_state{i_mpc}(round(t_mpc_span(j)/dt)+1,1:2)',30,'k');
        end
        %Final outer geodesic ball
        E_end = (sqrt(E_start)*exp(-lambda*delta) + d_bar*(1-exp(-lambda*delta)))^2;
        Ellipse_plot(M_ccm_pos_unscaled*(1/E_end),MPC_state{i_mpc}(round(delta/dt)+1,1:2)',30,'r',0.2);
%         Ellipse_plot(M_ccm_pos,MPC_state{i_mpc}(round(delta/dt)+1,1:2)',30,'k');
    end
else
    d0 = sqrt(geo_energy(1,1));
    for i = 1:round(delta/dt_sim):length(solve_t)
        time = solve_t(i);
        d_i = d0*exp(-lambda*time) + d_bar*(1-exp(-lambda*time));
        Ellipse_plot(M_ccm_pos_unscaled*(1/d_i^2),MP_state(1+(i-1)*(dt_sim/dt),1:2),30,'k');
    end
end

%Plot actual trajectory
plot(x_act(:,1),x_act(:,2),'k-','linewidth',2);
plot(x_act(1:round(delta/dt_sim):end,1),x_act(1:round(delta/dt_sim):end,2),'ko','markersize',7,'markerfacecolor','k');

%Final condition
Ellipse_plot(P(1:2,1:2)*(1/(alpha)), x_eq(1:2),30,'k');
xlabel('$X$','interpreter','latex'); 
ylabel('$Z$','interpreter','latex');
set(findall(gcf,'type','text'),'FontSize',32);set(gca,'FontSize',32)
grid on; 
% axis equal

%% State Trajectory
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
 
%Geod distances
figure()
plot(solve_t(1:end-1),sqrt(geo_energy(:,1)),'b-','linewidth',2);
hold on
if (~track_traj)
%     plot(solve_t(1:end-1),d_bar*ones(T_steps,1),'r-','linewidth',2);
    for i_mpc = 1:T_steps_MPC
        E_start = geo_energy(1+(i_mpc-1)*(delta/dt_sim),2);
        t_start = solve_t(1+(i_mpc-1)*(delta/dt_sim));
        t_mpc_span = 0:dt_sim:delta;
        d_mpc = sqrt(E_start)*exp(-lambda*t_mpc_span) + d_bar*(1-exp(-lambda*t_mpc_span));
        plot(t_mpc_span+t_start,d_mpc,'r-','linewidth',2);
    end
    plot(solve_t(1:end-1),sqrt(geo_energy(:,2)),'bo','markersize',10,'markerfacecolor','b','linewidth',2);
else
    plot(solve_t(1:end-1),d0*exp(-lambda*solve_t(1:end-1)) + d_bar*(1-exp(-lambda*solve_t(1:end-1))),'r-','linewidth',2);
end
grid on
legend('d(x^{*},x)','RCI bound');
xlabel('Time [s]');
set(findall(gcf,'type','text'),'FontSize',32);set(gca,'FontSize',32)

%MPC final time
figure()
plot(MPC_time(:,1),MPC_time(:,2),'ro','markersize',10,'markerfacecolor','r');
hold on
plot(MPC_time(:,1),min([MPC_time(:,1)+delta,Tp*ones(T_steps_MPC,1)],[],2),'bo','markersize',10,'markerfacecolor','b');
grid on
xlabel('Time [s]');ylabel('t_i + T');
legend('Optimal re-join time','Minimum re-join time');
set(findall(gcf,'type','text'),'FontSize',32);set(gca,'FontSize',32)


save('quad_sim.mat');