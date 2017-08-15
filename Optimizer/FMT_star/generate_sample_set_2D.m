function V = generate_sample_set_2D(x_lim,y_lim,goal_pos,goal_r,N,N_goal,EPS)

%% Generate sample set
   
freeCoords = [-x_lim+2*x_lim*rand(N,1), -y_lim+2*y_lim*rand(N,1)];

rand_r = goal_r*rand(N_goal,1);
theta_r = 2*pi*rand(N_goal,1);

goalCoords = kron(ones(N_goal,1),goal_pos')+[rand_r.*cos(theta_r),rand_r.*sin(theta_r)];

randCoords = [freeCoords;
              goalCoords];
maxNodes = N+N_goal;          

%% Generate node struct

q_empty.coord = [0 0]; %coords
q_empty.cost = 0; %accumulated cost 
q_empty.parent = 0; %parent
q_empty.idx = 0; %index in node set
q_empty.N = []; %nearest neighbor node set

V = repmat(q_empty,maxNodes,1);
for n = 1:maxNodes
   V(n).coord = randCoords(n,:); 
   V(n).idx = n+1; %to accommodate start node 
end

%% Compute nearest neighbor set

for n = 1:maxNodes
   V(n).N = Near(V,V(n),EPS); 
end

end