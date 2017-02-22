% function plot_quad_movie(t,x_act,MPC_state,M_ccm_pos,x_eq,obs,t_step)
close all;
clc;

t_step = (0.004/dt_sim);

[U_pos,S_pos,V_pos] = svd(M_ccm_pos);
S_new = (sqrt(S_pos\eye(2)) + 0*eye(2))^2\eye(2);

M_ccm_pos_infl = U_pos*S_new*V_pos';
%% Setup geometry
len = 0.25;

quad_bound = [[0.8*len;0.05],...
             [len;0.1],...
             [len;-0.1],...
             [0.8*len;-0.05],...
             [-0.8*len;-0.05],...
             [-len;-0.1],...
             [-len;0.1],...
             [-0.8*len;0.05],...
             [-0.05; 0.05],...
             [0;0.08],...
             [0.05;0.05]];

J = size(quad_bound,2);

%% Plot paths

figure();
%actual
% plot(x_act(:,1),x_act(:,2),'r-','linewidth',2);
plot(MP_state(:,1),MP_state(:,2),'r-','linewidth',2);
hold on;
%obstacles
for i = 1:obs.n_obs
    Ellipse_plot(eye(2)*(1/(obs.r(i))^2),obs.pos(:,i),25,'k',1);
end
% for i = 1:length(MPC_state) %each cell is new MPC segment
%     plot(MPC_state{i}(1:round(delta/dt_sim),1),MPC_state{i}(1:round(delta/dt_sim),2),'b-','linewidth',2); %resolve after delta
%     for j = 1:(delta/dt_sim):round(delta/dt_sim)+1
%         Ellipse_plot(M_ccm_pos_infl,MPC_state{i}(j,1:2)',25,'b',0.05);
%     end
% end
for i = 1:0.5/dt:length(MP_state)
    Ellipse_plot(M_ccm_pos_infl,MP_state(i,1:2)',25,'b',0.05);
end
% Ellipse_plot(M_ccm_pos,x_eq(1:2),25,'r');
line([-5 -5],[-5, 5],'color','k','linewidth',2);
line([-5  5],[ 5, 5],'color','k','linewidth',2);
line([ 5  5],[ 5,-5],'color','k','linewidth',2);
line([ 5 -5],[-5,-5],'color','k','linewidth',2);

i = 1;
quad_p = quad_bound;
R = [cos(x_act(i,3)), sin(x_act(i,3));
    -sin(x_act(i,3)), cos(x_act(i,3))]';

for j = 1:J
    quad_p(:,j) = x_act(i,1:2)' + R*quad_bound(:,j);
end
hp = patch(quad_p(1,:),quad_p(2,:),'k','FaceAlpha',0.8,'linewidth',2);
hold off

start_patch = [-5,-4,-4,-5 ;
               -5,-5,-4.5,-4.5];
patch(start_patch(1,:),start_patch(2,:),'g','FaceAlpha',0.5,'linewidth',2); 

end_patch = [5,4,4,5;
             5,5,4,4];
patch(end_patch(1,:),end_patch(2,:),'r','FaceAlpha',0.5,'linewidth',2); 

% xl = get(gca,'Xlim');
% yl = get(gca,'Ylim');
% xl = 1.1*xl;
% yl = 1.1*yl;
% set(gca,'Xlim',xl);
% set(gca,'Ylim',yl);
xlim(1.1*[-5,5]); ylim(1.1*[-5,5]);

% grid on;
grid off

axis manual;
% pause;
for i = t_step:t_step:length(solve_t)
    
    quad_p = quad_bound;
    R = [cos(x_act(i,3)), sin(x_act(i,3));
        -sin(x_act(i,3)), cos(x_act(i,3))]';
    
    for j = 1:J
        quad_p(:,j) = x_act(i,1:2)' + R*quad_bound(:,j);
    end
    hold on;
    set(hp,'XData',quad_p(1,:),'YData',quad_p(2,:));
    hold off;
    
    title(sprintf('t = %.2f',solve_t(i)));
    drawnow;
    pause(0.001);
       
    
end
% end