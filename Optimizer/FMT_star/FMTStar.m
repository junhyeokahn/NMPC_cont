function path = FMTStar(V,EPS,start_pos,d_start,goal_pos,goal_r,obs)

%% Add Start node to set of nodes
q_start.coord = start_pos';
q_start.cost = d_start;
q_start.parent = 0;
q_start.idx = 1;
q_start.N = [];

V = [q_start;V];

V(1).N = Near(V,q_start,EPS);
start_N = V(1).N;
for j = 1:length(start_N)
   V(start_N(j)).N = [1;V(start_N(j)).N]; % include node 1 
end

%% Initialize
V_unvisited = ones(length(V),1); V_unvisited(1) = 0;
V_open = Inf(length(V),1); V_open(1) = 1;
C_open = Inf(length(V),1); C_open(1) = 0;
z = V(1);
z_in_goal = 0;

%% Search

tic
while ~z_in_goal
    
    N_z = z.N;
    
    idx_open_new = [];
    C_open_new = [];

    for i = 1:length(N_z)
        if V_unvisited(N_z(i)) %N_z(i) in V_unvisited
            x = V(N_z(i)); 
            N_x = x.N;
            c_min = Inf; y_min = Inf;
            for j = 1:length(N_x)
                if V_open(N_x(j))< Inf
                    y = V(N_x(j));
%                     c = y.cost + dist_3d(y.coord,x.coord); 
                      c = dist_geod(y,x,EPS);
                    if c < c_min
                        c_min = c;
                        y_min = N_x(j); %store min node index
                    end
                end
            end
            if (y_min < Inf)  && (checkCollision(V(y_min),x,obs,EPS))
                %add x to graph
                V(N_z(i)).parent = y_min; V(N_z(i)).cost = c_min; 
                %add x to V_open_new
                idx_open_new = [idx_open_new;N_z(i)]; C_open_new = [C_open_new;c_min];
                %remove x from V_unvisited
                V_unvisited(N_z(i)) = 0;
            end
        end
    end
    %add V_open_new to V_open
    V_open(idx_open_new) = 1;
    C_open(idx_open_new) = C_open_new;
    
    %remove z from V_open 
    V_open(z.idx) = Inf;
    
    %find next min cost node
    [~,idx] = min(C_open.*V_open);
    
    z = V(idx);
    z_in_goal = (norm(goal_pos-z.coord') <= goal_r);
    
end
FMT_time = toc

%% Plot tree
% figure(2)
% for i = 2:length(V)
%     p = V(i).parent;
%     if (p~=0)
%         line([V(i).coord(1),V(p).coord(1)],[V(i).coord(2),V(p).coord(2)],[V(i).coord(3),V(p).coord(3)],'Color','k');
%     end
% end

%% Recover optimal path

% Search backwards from goal to start to find the optimal least cost path
figure(1)
%Obstacles
for i_ob = 1:obs.n_obs
    Ellipse_plot(eye(2)*(1/(obs.r(i_ob))^2),obs.pos(:,i_ob), 25,'r',1);
end
q_end = z;
path(1) = q_end;
while q_end.parent ~= 0
    start = q_end.parent;
    line([q_end.coord(1), V(start).coord(1)], [q_end.coord(2), V(start).coord(2)], 'Color', 'b', 'LineWidth', 2);
    hold on
    q_end = V(start);
    path = [path,q_end];
end

%flip to give start to end
path = fliplr(path);

end

