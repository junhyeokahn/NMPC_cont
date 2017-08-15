function isValid = checkCollision(q_start,q_end,obs,EPS)

isValid = 1;

Tp = norm(q_end.coord-q_start.coord)/1.7; %estimated total time
tau = 0:0.1:Tp;
q_traj = [linspace(q_start.coord(1),q_end.coord(1),length(tau));
          linspace(q_start.coord(2),q_end.coord(2),length(tau))];

E_time_bound = (q_start.cost*exp(-EPS.lambda*tau) + EPS.d_bar*(1-exp(-EPS.lambda*tau))).^2;
for t = 1:length(tau)
    for i = 1:obs.n_obs
        S_new = (sqrt(E_time_bound(t)*(obs.S\eye(2))) + (obs.r(i))*eye(2))^2\eye(2);
        M = obs.U*S_new*obs.V';
        if ((q_traj(:,t)-obs.pos(:,i))'*M*(q_traj(:,t)-obs.pos(:,i)) <= 1)
            isValid = 0;
            return;
        end
    end
end

end