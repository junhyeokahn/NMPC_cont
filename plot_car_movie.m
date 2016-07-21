function plot_car_movie(x_des,y_des,t,x_act,t_step)


tri_bound = [[0.2;0],0.1*[cosd(-120);sind(-120)],0.1*[cosd(120);sind(120)]];

figure();
plot(x_des,y_des,'r-','linewidth',1);
hold on;

xl = get(gca,'Xlim');
yl = get(gca,'Ylim');

xl = 1.1*xl;
yl = 1.1*yl;

i = 1;
tri_p = tri_bound;
R = [cos(x_act(i,4)), sin(x_act(i,4));
    -sin(x_act(i,4)), cos(x_act(i,4))]';

for j = 1:3
    tri_p(:,j) = x_act(i,1:2)' + R*tri_bound(:,j);
end
hp = patch(tri_p(1,:),tri_p(2,:),'k','FaceAlpha',0.1,'linewidth',2);
hold off

set(gca,'Xlim',xl);
set(gca,'Ylim',yl);

grid on;

axis manual;

for i = t_step:t_step:length(t)
    
    tri_p = tri_bound;
    R = [cos(x_act(i,4)), sin(x_act(i,4));
        -sin(x_act(i,4)), cos(x_act(i,4))]';
    
    for j = 1:3
        tri_p(:,j) = x_act(i,1:2)' + R*tri_bound(:,j);
    end
    hold on;
    set(hp,'XData',tri_p(1,:),'YData',tri_p(2,:));
    hold off;
    
    title(sprintf('t = %.2f',t(i)));
    drawnow;
       
    
end
end