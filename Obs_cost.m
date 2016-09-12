function J_ob = Obs_cost(xu,n,N,obs)

global US_B;

J_obs = zeros(N+1,obs.n_obs);

for i = 1:obs.n_obs
    o_pos = obs.pos(:,i);
    o_r = obs.r(i) + obs.infl;
    decay = obs.decay(i);
    
    for k = 1:N+1
        r_k = norm(xu(1+(k-1)*n:2+(k-1)*n)-o_pos);
        J_obs(k,i) = 2*exp(-decay*(r_k/o_r)^2);
    end
end

US_B = J_obs;

J_ob = sum(sum(J_obs));


end