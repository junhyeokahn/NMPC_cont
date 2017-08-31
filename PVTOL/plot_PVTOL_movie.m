% function plot_quad_movie(t,x_act,MPC_state,M_ccm_pos,x_eq,obs,t_step)

t_step = 0.05/dt;

plot_time_var = 0;

[U_pos,S_pos,V_pos] = svd(M_ccm_pos);
S_new = (sqrt(S_pos\eye(2)) + len*eye(2))^2\eye(2);

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

%Initial Motion plan
plot(MP_state(:,1),MP_state(:,2),'r-','linewidth',2);
hold on;

%Maximal RCI tube on MP
for i = 1:0.1/dt:length(MP_state)
    Ellipse_plot(M_ccm_pos,MP_state(i,1:2)',25,'r',0.05);
end

%obstacles
for i = 1:obs.n_obs
    Ellipse_plot(eye(2)*(1/(obs.r(i)+len)^2),obs.pos(:,i),25,'k',1);
end

%MPC
h_mpc = plot(MPC_state{1}(:,1),MPC_state{1}(:,2),'b-','linewidth',2);
t_vec = 0:0.1:T_mpc;
E_start = geo_energy(1,2);
E_bnd = (sqrt(E_start)*exp(-lambda*t_vec) + d_bar*(1-exp(-lambda*t_vec))).^2;

%Online RCI tube on MPC traj
S_new = (sqrt(E_bnd(1)*(obs_mpc.S\eye(2))) + 0*eye(2))^2\eye(2);
if plot_time_var
    h_rci = Ellipse_plot(obs_mpc.U*S_new*obs_mpc.V',MPC_state{1}(1,1:2)',25,'b',0.05);
else
    h_rci = Ellipse_plot(M_ccm_pos,MPC_state{1}(1,1:2)',25,'b',0.05);
end
h_rci = repmat(h_rci,length(t_vec),1);
for j = 2:length(t_vec)
    if plot_time_var
        S_new = (sqrt(E_bnd(j)*(obs_mpc.S\eye(2))) + 0*eye(2))^2\eye(2);
        h_rci(j) = Ellipse_plot(obs_mpc.U*S_new*obs_mpc.V',MPC_state{1}(1+(j-1)*(0.1/dt),1:2)',25,'b',0.1);
    else
        h_rci(j) = Ellipse_plot(M_ccm_pos,MPC_state{1}(1+(j-1)*(0.1/dt),1:2)',25,'b',0.1);
    end
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

xlim(1.1*[-5,5]); ylim(1.1*[-5,5]);

% grid on;
grid off

axis manual;

%% Setup Movie
record_vid = 0;

if (record_vid)
    writerObj = VideoWriter('PVTOL_sim.mp4');
    writerObj.FrameRate = 1/(t_step*(solve_t(2)-solve_t(1)));
    writerObj.Quality = 100;
    open(writerObj);
    set(gcf, 'renderer', 'zbuffer')
end
%%

keyboard;
i_mpc =  1;
for i = t_step:t_step:length(solve_t)
    tic
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
%     pause(0.001);
    
    if (floor(solve_t(i)/delta)+1 > i_mpc)
        i_mpc = i_mpc + 1;
        set(h_mpc,'XData',MPC_state{i_mpc}(:,1),'YData',MPC_state{i_mpc}(:,2));
        E_start = geo_energy(1+(i_mpc-1)*(delta/dt_sim),2);
        E_bnd = (sqrt(E_start)*exp(-lambda*t_vec) + d_bar*(1-exp(-lambda*t_vec))).^2;
        for j = 1:length(t_vec)
            delete(h_rci(j));
            if plot_time_var
                S_new = (sqrt(E_bnd(j)*(obs_mpc.S\eye(2))) + 0*eye(2))^2\eye(2);
                h_rci(j) = Ellipse_plot(obs_mpc.U*S_new*obs_mpc.V',MPC_state{i_mpc}(1+(j-1)*(0.1/dt),1:2)',25,'b',0.05);
            else
                h_rci(j) = Ellipse_plot(M_ccm_pos,MPC_state{i_mpc}(1+(j-1)*(0.1/dt),1:2)',25,'b',0.05);
            end
        end
%         pause(0.25);
    end      
    pause(t_step*(solve_t(2)-solve_t(1))-toc);
    
    if (record_vid)
        thisFrame = getframe(gcf);
        %Write this frame out to a new video file.
        writeVideo(writerObj, thisFrame);
    end
    
end

if (record_vid)
    close(writerObj);
end
% end