function plot_quad_movie(x_act,t,t_step)

len = 0.25;

quad_bound = [[0.9*len;0.05],...
             [len;0.1],...
             [len;-0.1],...
             [0.9*len;-0.05],...
             [-0.9*len;-0.05],...
             [-len;-0.1],...
             [-len;0.1],...
             [-0.9*len;0.05],...
             [-0.05; 0.05],...
             [0;0.08],...
             [0.05;0.05]];

J = size(quad_bound,2);

figure();
plot(x_act(:,1),x_act(:,2),'r-','linewidth',2);
hold on;

i = 1;
quad_p = quad_bound;
R = [cos(x_act(i,3)), sin(x_act(i,3));
    -sin(x_act(i,3)), cos(x_act(i,3))]';

for j = 1:J
    quad_p(:,j) = x_act(i,1:2)' + R*quad_bound(:,j);
end
hp = patch(quad_p(1,:),quad_p(2,:),'k','FaceAlpha',0.1,'linewidth',2);
hold off

set(gca,'Xlim',[min(x_act(:,1))-2*len,max(x_act(:,1))+2*len]);
set(gca,'Ylim',[min(x_act(:,2))-2*len,max(x_act(:,2))+2*len]);

grid on;

axis manual;
pause;
for i = t_step:t_step:length(t)
    
    quad_p = quad_bound;
    R = [cos(x_act(i,3)), sin(x_act(i,3));
        -sin(x_act(i,3)), cos(x_act(i,3))]';
    
    for j = 1:J
        quad_p(:,j) = x_act(i,1:2)' + R*quad_bound(:,j);
    end
    hold on;
    set(hp,'XData',quad_p(1,:),'YData',quad_p(2,:));
    hold off;
    
    title(sprintf('t = %.2f',t(i)));
    drawnow;
    pause(0.001);
       
    
end
end