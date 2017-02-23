%2D State Plot
% P_inv = diag([39.0251,486.0402]);
figure()
for i_mpc = 1:T_steps_MPC
    plot(MPC_state{i_mpc}(:,1),...
         MPC_state{i_mpc}(:,2),'r-','linewidth',2);
      hold on;
    Ellipse_plot(M_ccm*(1/d_bar^2),MPC_state{i_mpc}(1,1:2),25,'k');
end
plot(x_act(:,1),x_act(:,2),'b-','linewidth',2);    

Ellipse_plot((1/alpha)*P,[0;0],25,'r');
xlabel('x_1'); ylabel('x_2');
set(findall(gcf,'type','text'),'FontSize',32);set(gca,'FontSize',32)
grid on; 