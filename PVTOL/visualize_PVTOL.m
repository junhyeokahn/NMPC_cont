
% close all

if strcmp(select_path,'ad')
    MP_state = MP_ad_state;
    MP_ctrl = MP_ad_ctrl;
end

figure(1); 
hold on
%Obstacles
for i_ob = 1:obs.n_obs
    Ellipse_plot(eye(2)*(1/(obs.r(i_ob)+len)^2),obs.pos(:,i_ob), 25,'r',1);
end
%RCI set
for i = 1:(1/dt):length(MP_state)
    Ellipse_plot(M_ccm_pos,MP_state(i,1:2)',25,'k');
end
%Terminal set
Ellipse_plot(P(1:2,1:2)*(1/alpha),x_eq(1:2),25,'r');

%Full MP traj
plot(MP_state(:,1),MP_state(:,2),'b-','linewidth',1);
quiver(MP_state(1:(0.5/dt):end-1,1),MP_state(1:(0.5/dt):end-1,2),...
       (MP_ctrl(1:(0.5/dt):end-1,1)+MP_ctrl(1:(0.5/dt):end-1,2)).*sin(MP_state(1:(0.5/dt):end-1,3)),...
       (MP_ctrl(1:(0.5/dt):end-1,1)+MP_ctrl(1:(0.5/dt):end-1,2)).*cos(MP_state(1:(0.5/dt):end-1,3)));

plot(test_state(1),test_state(2),'go','markersize',10,'markerfacecolor','g');
grid on
axis equal
xlabel('X'); ylabel('Z','interpreter','latex');
set(findall(gcf,'type','text'),'FontSize',32);set(gca,'FontSize',32)