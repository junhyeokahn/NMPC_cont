function plot_quad_movie(x_des,y_des,z_des,t,x_act,t_step,M,goal,obs)

yaw_idx = 9;
pitch_idx = 8;
roll_idx = 7;

%quadrotor body 
bound = [0.1,0.1, 0;
           0.1, -0.1, 0;
           -0.1,-0.1,0;
           -0.1,0.1,0]';

%rotor base       
N_circ = 20;       
base_circle = circle([0;0;0],0.05,N_circ); 

%axes base
ax = 0.3*eye(3);

%workspace
figure();
plot3(x_des,y_des,z_des,'r-','linewidth',2);
hold on
plot3(x_act(:,1),x_act(:,2),x_act(:,3),'b-','linewidth',2);
for i = 1:t_step:length(t)
    Ellipse_plot(M,[x_des(i);y_des(i);z_des(i)],30,'r',0.4,3);
end
set(gca,'ZDir','Reverse');
set(gca,'YDir','Reverse');
hold on;

goal.plot('color','green','alpha',0.3)
hold on
for i = 1:length(obs)
    obs{i}.plot('color','black','alpha',0.3);
end

xl = get(gca,'Xlim');
yl = get(gca,'Ylim');
zl = get(gca,'Zlim');

xl = 0.75*xl;
yl = 0.75*yl;
zl = 0.75*zl;

i = 1;

%define points
q_bound = bound;
q_rotors = cell(4,1);

R = rot_matrix(x_act(i,roll_idx),x_act(i,pitch_idx),x_act(i,yaw_idx));
q_ax = R*ax + kron(ones(1,3),x_act(i,1:3)');

c_bound = R*base_circle;
for j = 1:4
    q_bound(:,j) = x_act(i,1:3)' + R*bound(:,j);
    q_rotors{j} = kron(ones(1,N_circ),q_bound(:,j)) + ...
                   c_bound;
end

%draw objects
h_rotors = zeros(4,1);
for j = 1:4
    h_rotors(j) = patch(q_rotors{j}(1,:),q_rotors{j}(2,:),q_rotors{j}(3,:),'k','FaceAlpha',0.1,'linewidth',2);
end

h_ax = zeros(3,1);
col_ax = {'r','b','g'};
for j = 1:3
    h_ax(j) = line([x_act(i,1);q_ax(1,j)],[x_act(i,2),q_ax(2,j)],[x_act(i,3),q_ax(3,j)],...
                     'color',col_ax{j},'linewidth',2);
end

h_arms = zeros(2,1);
for j = 1:2
    h_arms(j) = line([q_bound(1,j);q_bound(1,j+2)],...
                     [q_bound(2,j);q_bound(2,j+2)],...
                     [q_bound(3,j);q_bound(3,j+2)],...
                     'color','k','linewidth',2);
end

hold off

set(gca,'Xlim',xl);
set(gca,'Ylim',yl);
set(gca,'Zlim',zl);
view([43 25]);
xlabel('X'); ylabel('Y');
grid on;

axis equal;

view_center = [0;0;0];
rel_view_az = atan2(x_act(i,1)-view_center(1),x_act(i,2)-view_center(2));
rel_view_el = atan2(view_center(3)-x_act(i,3),norm(x_act(i,1:2)-view_center(1:2)));
az = rel_view_az*(180/pi) - 10;
el = rel_view_el*(180/pi) + 5;
view([az,el]);

%% Setup Movie
record_vid = 0;

if (record_vid)
    writerObj = VideoWriter('Quad_sim_poly.mp4');
    writerObj.FrameRate = 1/(t_step*(t(2)-t(1)));
    writerObj.Quality = 100;
    open(writerObj);
    set(gcf, 'renderer', 'zbuffer')
end

%% Record

keyboard;
for i = 1:t_step:length(t)
    tic
%     set(gcf, 'Units', 'Normalized', 'Outerposition', [0, 0, 1, 1]);
    R = rot_matrix(x_act(i,roll_idx),x_act(i,pitch_idx),x_act(i,yaw_idx));
    q_ax = R*ax + kron(ones(1,3),x_act(i,1:3)');
    
    c_bound = R*base_circle;
    for j = 1:4
        q_bound(:,j) = x_act(i,1:3)' + R*bound(:,j);
        q_rotors{j} = kron(ones(1,N_circ),q_bound(:,j)) + ...
                        c_bound;
    end
    hold on;
    
    for j = 1:4
        set(h_rotors(j),'XData',q_rotors{j}(1,:),'YData',q_rotors{j}(2,:),'ZData',q_rotors{j}(3,:));
    end
    for j = 1:3
        set(h_ax(j),'XData',[x_act(i,1);q_ax(1,j)],...
                      'YData',[x_act(i,2),q_ax(2,j)],...
                      'ZData',[x_act(i,3),q_ax(3,j)]);
    end
    for j = 1:2
        set(h_arms(j),'XData',[q_bound(1,j);q_bound(1,j+2)],...
                      'YData',[q_bound(2,j),q_bound(2,j+2)],...
                      'ZData',[q_bound(3,j),q_bound(3,j+2)]);
    end
    hold off;
    
    rel_view_az = atan2(x_act(i,1)-view_center(1),x_act(i,2)-view_center(2));
    rel_view_el = atan2(1,norm(x_act(i,1:2)-view_center(1:2)));
    az = rel_view_az*(180/pi) - 10;
    el = rel_view_el*(180/pi) + 5;
    view([az,el]);
    
    title(sprintf('t = %.2f',t(i)));
    drawnow;
    pause(t_step*(t(2)-t(1))-toc);
    
    if (record_vid)
        thisFrame = getframe(gcf);
        %Write this frame out to a new video file.
        writeVideo(writerObj, thisFrame);
    end
        
end

if (record_vid)
    close(writerObj);
end
end