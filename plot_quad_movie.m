function plot_quad_movie(x_des,y_des,z_des,t,x_act,t_step,n)

if n==14
    yaw_idx = 11;
    pitch_idx = 10;
    roll_idx = 9;
else
    yaw_idx = 9;
    pitch_idx = 8;
    roll_idx = 7;
end
bound = [0.2,0.2, 0;
           0.2, -0.2, 0;
           -0.2,-0.2,0;
           -0.2,0.2,0]';

ax = 0.3*eye(3);

figure();
plot3(x_des,y_des,z_des,'r-','linewidth',1);
set(gca,'ZDir','Reverse');
set(gca,'YDir','Reverse');
hold on;

xl = get(gca,'Xlim');
yl = get(gca,'Ylim');
zl = get(gca,'Zlim');

xl = 1.5*xl;
yl = 1.5*yl;
zl = 1.5*zl;

i = 1;
q_bound = bound;
R = rot_matrix(x_act(i,yaw_idx),x_act(i,pitch_idx),x_act(i,roll_idx));
q_ax = R*ax + kron(ones(1,3),x_act(i,1:3)');

for j = 1:4
    q_bound(:,j) = x_act(i,1:3)' + R*bound(:,j);
end
h_bound = patch(q_bound(1,:),q_bound(2,:),q_bound(3,:),'k','FaceAlpha',0.1,'linewidth',2);
h_line = zeros(3,1);
col_ax = {'r','b','g'};
for j = 1:3
h_line(j) = line([x_act(i,1);q_ax(1,j)],[x_act(i,2),q_ax(2,j)],[x_act(i,3),q_ax(3,j)],...
                 'color',col_ax{j},'linewidth',2);
end
h_corn = plot3(q_bound(1,:)',q_bound(2,:)',q_bound(3,:)',...
               'bo','markersize',10,'markerfacecolor','g');
hold off

set(gca,'Xlim',xl);
set(gca,'Ylim',yl);
set(gca,'Zlim',zl);
view([39 42]);
xlabel('X'); ylabel('Y');
grid on;

axis manual;

for i = t_step:t_step:length(t)
    
    q_bound = bound;
    R = rot_matrix(x_act(i,yaw_idx),x_act(i,pitch_idx),x_act(i,roll_idx));
    q_ax = R*ax + kron(ones(1,3),x_act(i,1:3)');
    for j = 1:4
        q_bound(:,j) = x_act(i,1:3)' + R*bound(:,j);
    end
    hold on;
    set(h_bound,'XData',q_bound(1,:),'YData',q_bound(2,:),'ZData',q_bound(3,:));
    for j = 1:3
        set(h_line(j),'XData',[x_act(i,1);q_ax(1,j)],...
                      'YData',[x_act(i,2),q_ax(2,j)],...
                      'ZData',[x_act(i,3),q_ax(3,j)]);
    end
    set(h_corn,'XData',q_bound(1,:)',...
               'YData',q_bound(2,:)',...
               'ZData',q_bound(3,:)');
    hold off;
    
    title(sprintf('t = %.2f',t(i)));
    drawnow;
    
end
end